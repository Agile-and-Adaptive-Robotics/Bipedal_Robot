"""Hunt C: fit unclamped four-cell joint PF activity with frozen A/B inputs.

Inactive cells may hyperpolarize. The loss compares rectified PF activity,
rather than forcing inactive membrane voltage to zero and compensating its
reciprocal inhibition. No motor output or physical curriculum is certified here.
"""
import os,sys,json,time,hashlib,shutil
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import numpy as np
from numba import njit
from scipy.optimize import least_squares
import controller as C
from pf_targets import mono_command_targets

H=Path(__file__).resolve().parent
source=H/'runs/feasible_hunt_20261030_v1'
starter=H/'runs/pf_inhibition_015_20261030'
O=H/'runs'/os.environ.get('HUNT_NATIVE_PF_CAMPAIGN','native_joint_pf_20261030_v1');O.mkdir(exist_ok=True);W=O/'weights';W.mkdir(exist_ok=True)
slow_weight=float(os.environ.get('HUNT_NATIVE_SLOW_WEIGHT','1'))
tie_inhibition=os.environ.get('HUNT_NATIVE_PF_TIE_INHIBITION','0')=='1'
target_mode=os.environ.get('HUNT_NATIVE_TARGET','net_active_torque')
if (O/'NATIVE_PF_COMPLETE.json').exists():raise SystemExit('Existing native PF fit retained')
for name in ('pf_weights.json','timing_weights.json','circuit_options.json'):shutil.copyfile(starter/'weights'/name,W/name)
acts=np.load(H/'runs/targets/targets.npz',allow_pickle=True)['act_names'].tolist()
n=C.build(acts,options=dict(weights_directory=str(W),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True));n.apply_weights(1.);c=n.compiled

@njit(cache=True,nogil=True)
def simulate(x,pre,post,g,rev,lo,ran,tf,gm,rest,bias,out):
    state=rest.copy();following=np.empty(4)
    for k in range(len(x)):
        for j in range(4):following[j]=-gm[j]*(state[j]-rest[j])+bias[j]
        for e in range(len(g)):
            value=state[pre[e]] if pre[e]<4 else x[k,pre[e]-4]
            saturation=min(1.,max(0.,(value-lo[e])/ran[e]))
            following[post[e]]+=g[e]*saturation*(rev[e]-state[post[e]])
        for j in range(4):
            state[j]+=tf[j]*following[j]
            if not np.isfinite(state[j]) or abs(state[j])>100.:
                out[:]=np.nan;return
            out[k,j]=state[j]

start=time.perf_counter();weights=json.loads((W/'pf_weights.json').read_text());results={}
for side in ('r','l'):
    for joint in ('HIP','KNEE'):
        cells=[f'PF_{joint}-E_{side}',f'PF_{joint}-F_{side}',f'PF_IN_{joint}-E_{side}',f'PF_IN_{joint}-F_{side}']
        ids=np.array([n.idx[name] for name in cells]);owned=set(ids.tolist())
        edges=np.array([e for e in range(len(c.pre)) if c.post[e] in owned and (c.gmax[e]>0. or e in n.pfbank_edges)],dtype=int)
        external=sorted(set(c.pre[edges])-owned);mapping={j:i for i,j in enumerate(ids.tolist()+external)}
        pre=np.array([mapping[int(c.pre[e])] for e in edges]);post=np.array([mapping[int(c.post[e])] for e in edges]);g=c.gmax[edges].copy()
        variable=np.array([j for j,e in enumerate(edges) if e in n.pfbank_edges or e in n.pf_internal_edges]);assert len(variable)==50,len(variable)
        params=(c.reversal[edges],c.elo[edges],c.erange[edges],c.c.time_factor_membrane[ids],c.c.g_m[ids],c.c.V_rest[ids],c.c.i_b[ids])
        datasets=[];parity=[]
        for mode,T in [('afferented',1.233),('deafferented',2.466)]:
            teacher=np.load(source/('pf_'+mode)/'pf_teacher.npz');names=teacher['neuron_names'].tolist();cols=teacher['pf_cells'].tolist()
            x=np.vstack([c.c.V_rest[external][None,:],teacher['voltage'][:-1,[names.index(n.names[i]) for i in external]]])
            target=np.maximum(teacher['pf_target'][:,[cols.index(cell) for cell in cells[:2]]],0.)
            if target_mode=='mono_command':
                primitives=mono_command_targets(teacher['desired_mn'],acts,n.hunt_pairs)
                target=np.column_stack([primitives[cell] for cell in cells[:2]])
            reference=np.load(starter/'weights'/('pf_only_'+mode+'.npz'))
            expected=reference['pf'][:,[cols.index(cell) for cell in cells[:2]]]
            out=np.empty((len(x),4));simulate(x,pre,post,g,*params,out)
            error=float(np.max(abs(out[:,:2]-expected)));parity.append(error);assert error<1e-10,(side,joint,mode,error)
            keep=teacher['t']>=4*T;sample=np.flatnonzero(keep)[::10]
            datasets.append((x,target,out,sample,keep))
        free=np.arange(len(variable));inhibitory=np.array([j for j,e in enumerate(edges[variable]) if e in n.pf_internal_edges],dtype=int)
        if tie_inhibition:
            assert len(inhibitory)==2
            free=free[free!=inhibitory[1]]
        def expand(parameters):
            if not tie_inhibition:return parameters
            values=np.empty(len(variable));values[free]=parameters;values[inhibitory[1]]=values[inhibitory[0]];return values
        def residual(parameters):
            values=expand(parameters)
            candidate=g.copy();candidate[variable]=values;rows=[]
            for condition,(x,target,out,sample,keep) in enumerate(datasets):
                simulate(x,pre,post,candidate,*params,out)
                rows.append(np.sqrt(1. if condition==0 else slow_weight)*np.nan_to_num((np.maximum(out[sample,:2],0.)-target[sample]).ravel(),nan=1000.)/np.sqrt(len(sample)))
            rows.append(.003*values/20.)
            return np.concatenate(rows)
        seed=g[variable].copy()
        if tie_inhibition:seed[inhibitory]=seed[inhibitory].mean()
        seed=seed[free]
        baseline=float(np.dot(residual(seed),residual(seed)))
        fit=least_squares(residual,seed,bounds=(np.full(len(free),.001),np.full(len(free),20.)),
            max_nfev=100,ftol=1e-5,xtol=1e-5,gtol=1e-5,x_scale='jac') if np.all(seed>=.001) else None
        # Existing NNLS zeros are legal; move only the initial search point into bounds.
        if fit is None:fit=least_squares(residual,np.maximum(seed,.001001),bounds=(.001,20.),max_nfev=100,ftol=1e-5,xtol=1e-5,gtol=1e-5,x_scale='jac')
        selected=fit.x if np.dot(residual(fit.x),residual(fit.x))<baseline else seed
        values=expand(selected)
        if tie_inhibition:assert values[inhibitory[0]]==values[inhibitory[1]]
        candidate=g.copy();candidate[variable]=values;waveform={}
        for mode,(x,target,out,sample,keep) in zip(('afferented','deafferented'),datasets):
            simulate(x,pre,post,candidate,*params,out);rows={}
            for j,cell in enumerate(cells[:2]):
                actual=np.maximum(out[keep,j],0.);desired=target[keep,j]
                corr=float(np.corrcoef(actual,desired)[0,1]) if np.std(actual)>1e-8 else 0.
                modulation=float(np.ptp(actual)/max(float(np.ptp(desired)),1e-8))
                off=float(actual[desired<=np.percentile(desired,20)].mean()/max(float(desired.max()),1e-8))
                rows[cell]=dict(correlation=corr,modulation_ratio=modulation,relative_off_activity=off,
                    accepted=bool(corr>.5 and modulation>.5 and off<.5))
            waveform[mode]=rows;np.savez_compressed(O/(side+'_'+joint+'_'+mode+'.npz'),target=target,voltage=out)
        for j,value in zip(variable,values):
            e=edges[j];weights['edges'][n.names[c.pre[e]]+' -> '+n.names[c.post[e]]]=float(value)
        result=dict(native_reference_max_error=max(parity),baseline_loss=baseline,fitted_loss=float(np.dot(residual(selected),residual(selected))),
            evaluations=int(fit.nfev),free_parameters=len(free),tied_opposite_inhibition=tie_inhibition,waveform=waveform)
        results[side+'_'+joint]=result;print(side,joint,json.dumps(result),flush=True)
weights['scope']='Free four-cell joint PF activity fit across both frozen A/B conditions; includes learned reciprocal inhibitory conductances'
(W/'pf_weights.json').write_text(json.dumps(weights,indent=2))
report=dict(results=results,wall_s=time.perf_counter()-start,target_mode=target_mode,slow_target_loss_weight=slow_weight,tied_opposite_inhibition=tie_inhibition,accepted=all(row['accepted'] for r in results.values() for rows in r['waveform'].values() for row in rows.values()),
    scope='Native PF fitting check only; no motor, reflex, or physical curriculum pass',upstream=str(source),starter=str(starter))
(O/'NATIVE_PF_COMPLETE.json').write_text(json.dumps(report,indent=2))
(O/'sources').mkdir(exist_ok=True)
for path in (Path(__file__),H/'controller.py'):shutil.copyfile(path,O/'sources'/path.name)
print('NATIVE_PF_COMPLETE',json.dumps(dict(accepted=report['accepted'],wall_s=report['wall_s'])),flush=True)
