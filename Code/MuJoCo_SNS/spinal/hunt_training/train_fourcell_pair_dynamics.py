"""Hunt C: fit native paired-PF dynamics before MN and reflex refinement."""
import os,sys,json,time,shutil,hashlib
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
H=Path(__file__).resolve().parent;prior=H/'runs/fourcell_multicondition_20261029_v1';O=H/'runs/fourcell_dynamic_20261029_v1';O.mkdir(exist_ok=True);fits=O/'weights';fits.mkdir(exist_ok=True)
if (O/'TRAINING_COMPLETE.json').exists():print('Existing completed training retained:',O);raise SystemExit(0)
import numpy as np,optuna
from numba import njit
import controller as C
start=time.perf_counter()
if not (fits/'pf_weights.json').exists():shutil.copyfile(prior/'weights/pf_weights.json',fits/'pf_weights.json')
if not (O/'DYNAMIC_TRAINING.json').exists():
    shutil.copyfile(prior/'weights/bi_weights.json',fits/'bi_weights.json')
    acts=np.load(H/'runs/targets/targets.npz',allow_pickle=True)['act_names'].tolist()
    net=C.build(acts,options=dict(weights_directory=str(fits),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True));net.apply_weights(1.);c=net.compiled
    @njit(cache=True)
    def simulate(x,pre,post,g,rev,lo,ran,tf,gm,rest,bias,out):
        state=rest.copy();following=np.empty(4)
        for k in range(len(x)):
            for j in range(4):following[j]=-gm[j]*(state[j]-rest[j])+bias[j]
            for e in range(len(g)):
                value=state[pre[e]] if pre[e]<4 else x[k,pre[e]-4]
                sat=min(1.,max(0.,(value-lo[e])/ran[e]));following[post[e]]+=g[e]*sat*(rev[e]-state[post[e]])
            for j in range(4):state[j]+=tf[j]*following[j];out[k,j]=state[j]
    weights={};report={};optuna.logging.set_verbosity(optuna.logging.WARNING)
    for side in ('r','l'):
        cells=['PF_BI_JOINT_HAM_'+side,'PF_BI_JOINT_RECT_'+side,'PF_IN_BI_PAIR_HAM_'+side,'PF_IN_BI_PAIR_RECT_'+side]
        inputs=['PF_HIP-E_'+side,'PF_KNEE-F_'+side,'PF_HIP-F_'+side,'PF_KNEE-E_'+side]
        ids=np.array([net.idx[name] for name in cells]);mapping={net.idx[name]:i for i,name in enumerate(cells+inputs)}
        edges=np.flatnonzero(np.isin(c.post,ids));assert len(edges)==8,[(net.names[c.pre[e]],net.names[c.post[e]]) for e in edges]
        pre=np.array([mapping[c.pre[e]] for e in edges]);post=np.array([mapping[c.post[e]] for e in edges]);g=c.gmax[edges].copy()
        variable=np.array([j for j,e in enumerate(edges) if e in net.bi_edges]);assert len(variable)==6
        params=(c.reversal[edges],c.elo[edges],c.erange[edges],c.c.time_factor_membrane[ids],c.c.g_m[ids],c.c.V_rest[ids],c.c.i_b[ids])
        datasets=[];parity=[]
        for condition,T in [('afferented',1.233),('deafferented',2.466)]:
            z=np.load(prior/('mn_'+condition)/'mn_teacher.npz');names=z['neuron_names'].tolist();v=z['voltage'] if 'voltage' in z else z['v']
            x=np.vstack([np.zeros((1,4)),v[:-1,[names.index(name) for name in inputs]]]);expected=v[:,[names.index(name) for name in cells]]
            teacher=np.load(prior/('bi_'+condition)/'bi_teacher.npz');columns=teacher['pf_cells'].tolist();target=teacher['pf_target'][:,[columns.index(name) for name in cells[:2]]]
            out=np.empty((len(x),4));keep=z['t']>=4*T;simulate(x,pre,post,g,*params,out);error=float(np.max(abs(out-expected)));parity.append(error)
            assert error<1e-10,(side,condition,error)
            datasets.append((x,target,out,keep))
        def evaluate(values):
            trialg=g.copy();trialg[variable]=values;loss=0.
            for x,target,out,keep in datasets:
                simulate(x,pre,post,trialg,*params,out);loss+=float(np.mean((out[keep,:2]-target[keep])**2))/2
            return loss
        baseline=evaluate(g[variable]);study=optuna.create_study(direction='minimize',sampler=optuna.samplers.TPESampler(seed=291029))
        seed={str(j):max(float(g[e]),.001) for j,e in enumerate(variable)};study.enqueue_trial(seed)
        def objective(trial):
            values=np.array([trial.suggest_float(str(j),.001 if pre[e]>=4 else .02,5. if pre[e]>=4 else 1.,log=True) for j,e in enumerate(variable)])
            return evaluate(values)
        study.optimize(objective,n_trials=384)
        best=np.array([study.best_params[str(j)] for j in range(len(variable))]);best_loss=evaluate(best)
        if best_loss>=baseline:best=g[variable];best_loss=baseline
        for e,value in zip(variable,best):weights[net.names[c.pre[edges[e]]]+' -> '+net.names[c.post[edges[e]]]]=float(value)
        report[side]=dict(baseline_voltage_mse=baseline,fitted_voltage_mse=best_loss,max_reference_error=max(parity),trials=len(study.trials),parameters=best.tolist())
        print(side,json.dumps(report[side]),flush=True)
    (fits/'bi_weights.json').write_text(json.dumps(dict(edges=weights,scope='Native four-cell dynamics fitted across both frozen timing conditions; joint inputs and reciprocal inhibitory conductances learned.'),indent=2))
    (O/'DYNAMIC_TRAINING.json').write_text(json.dumps(dict(results=report,wall_s=time.perf_counter()-start,upstream=str(prior)),indent=2))
sources=O/'sources';sources.mkdir(exist_ok=True);shutil.copyfile(Path(__file__),sources/Path(__file__).name)
os.environ['HUNT_CAMPAIGN']=O.name
import train_fourcell_multicondition
