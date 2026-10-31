"""Hunt C/D: sequential conductance fits with timing fixed.

Teacher-forced PF/MN voltages are used only for offline conductance regression.
Validation and physical runs use native neural outputs with no target clamp.
"""
import os,sys,json,time
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;S=H.parent;sys.path[:0]=[str(H),str(S),str(S/'functional_gait2392')]
import numpy as np,mujoco
from scipy.optimize import lsq_linear
from scipy.ndimage import gaussian_filter1d
import controller as C,bsolve_ik as B,params as P
from muscle_map import EXTENSOR_STANCE_GROUPS
DT=.002;O=H/'runs/fitted';O.mkdir(parents=True,exist_ok=True)

def setup(cycles=8):
    z=np.load(H/'runs/targets/targets.npz',allow_pickle=True)
    m,d=B.build_air_model(False);acts=z['act_names'].tolist();net=C.build(acts);net.bind_mujoco(m,d);net.training_gate=1.
    t0,t1=.617,1.85;T=t1-t0;t=np.arange(int(cycles*T/DT))*DT;query=t0+t%T
    def interpolate(name):
        value=z[name];return np.stack([np.interp(query,z['t'],value.reshape(len(z['t']),-1)[:,k]) for k in range(value[0].size)],axis=1).reshape((len(t),)+value.shape[1:])
    q=interpolate('qpos');vel=interpolate('qvel');a=interpolate('activation')
    strength=z['rho']*(z['unit_force']*z['activation'])[:,None,:]
    base_pf={}
    for j,pair in net.hunt_pairs.items():
        row=C.JOINTS.index(j)
        for sign,cell in zip((1.,-1.),pair):
            if net.options.get('pf_target_mode','gross_muscle_torque')=='net_active_torque':
                positive=np.maximum(sign*strength[:,row,:].sum(axis=1),0.)
            else:positive=np.maximum(sign*strength[:,row,:],0.).sum(axis=1)
            target=np.interp(query,z['t'],positive)/max(float(np.percentile(positive,95)),1.)
            base_pf[cell]=gaussian_filter1d(np.clip(target,0.,1.)*2.75,10.,mode='wrap')
    return z,m,d,net,t,q,vel,a,base_pf

def afferents(net,u,d):
    c=net.compiled
    if not hasattr(net,'fit_encoder'):
        acts=net.targets_data['act_names'].tolist()
        net.fit_encoder=dict(ext=np.array([net.muscles[a].groups[0] in EXTENSOR_STANCE_GROUPS for a in acts]),
            rg=np.array([net.idx['RG_E_'+net.muscles[a].side] for a in acts]),
            ports={prefix:np.array([net.ports[prefix+'_'+a] for a in acts]) for prefix in ('Ia','II','Ib')})
    cache=net.fit_encoder;stance=np.clip(c.V[cache['rg']]/5.,0,1);ext=cache['ext']
    ranges=net.live[0].actuator_lengthrange;length=np.clip((d.actuator_length-ranges.mean(axis=1))/np.maximum(np.ptp(ranges,axis=1)/2,.001),-1.,1.)
    velocity=np.clip(d.actuator_velocity/P.AFF['ia_vel_ref'],-1.,1.);force=np.clip(abs(d.actuator_force)/net.Fmax,0,1.5)
    drive_fraction=float(np.clip(c.V[net.idx['DRIVE']]/P.E_HI,0.,1.))
    ia=P.AFF['ia_gain']*(1+P.MOD['ia']*drive_fraction)*np.where(ext,.5+.5*stance,1.)
    ii=P.AFF['ii_gain']*(1+P.MOD['ii']*drive_fraction)*np.where(ext,.3+.7*stance,.5)
    ib=P.AFF['ib_gain']*np.where(ext,1-.7*stance,1.)
    u[cache['ports']['Ia']]=np.maximum(ia*velocity,0.)
    u[cache['ports']['II']]=np.maximum(P.AFF['i0_ii']*(.3+.7*stance)+ii*length,0.)
    u[cache['ports']['Ib']]=.5*ib*force

def replay(kind,period_s=1.233):
    z,m,d,n,t,q,vel,a,pftarget=setup();c=n.compiled;N=len(t);u=n.make_inputs();u[n.ports['DRIVE']]=2.5
    u[n.ports['POSTURE']]=0.;keep=t>=4*period_s
    rows={};matrix={};rhs={};lookup={}
    if kind=='BI':
        for family in ('HAM','RECT'):
            for side in ('r','l'):
                indices=[i for i,act in enumerate(z['act_names']) if act.endswith('_'+side) and (act.startswith(('semimem','semiten','bifemlh')) if family=='HAM' else act.startswith('rect_fem'))]
                wave=np.mean(a[:,indices],axis=1)
                pftarget['PF_BI_JOINT_'+family+'_'+side]=gaussian_filter1d(wave/max(float(np.percentile(wave,95)),.05)*2.75,10.,mode='wrap')
    edges=n.pfbank_edges if kind=='PF' else n.bi_input_edges if kind=='BI' else n.pf_edges if kind=='MN' else n.reflex_edges
    if kind in ('PF','BI','MN'):n.reflex_gain[:]=0.;n.options['reciprocal_gain']=0.
    bypost={post:edges[c.post[edges]==post] for post in np.unique(c.post[edges])}
    others={post:np.flatnonzero((c.post==post)&~np.isin(np.arange(len(c.post)),ids)) for post,ids in bypost.items()}
    mn_columns={n.idx['MN_'+act]:i for i,act in enumerate(z['act_names'])}
    for post,ids in bypost.items():
        matrix[post]=np.zeros((int(keep.sum()),len(ids)));rhs[post]=np.zeros(int(keep.sum()));lookup[post]=ids
    trace=np.empty((N,len(c.V)));desired_mn=a*5.*(.8 if kind in ('PF','MN') else 1.)
    index=0
    for k in range(N):
        d.qpos[:]=q[k];n.arm_model.project(d);d.qvel[:]=vel[k];d.act[:]=n.expected_activation[k] if hasattr(n,'expected_activation') else a[k];mujoco.mj_forward(m,d)
        afferents(n,u,d)
        for side in n.sides:
            u[n.ports['HUNT_INIT_E_'+side]]=2. if side=='r' and k<50 else 0.
            u[n.ports['HUNT_INIT_F_'+side]]=2. if side=='l' and k<50 else 0.
        if kind in ('PF','BI'):
            for cell,y in pftarget.items():
                if kind=='PF' or cell.startswith('PF_BI_JOINT_'):c.c.V[n.idx[cell]]=y[k]
        else:
            for i,act in enumerate(z['act_names']):c.c.V[n.idx['MN_'+act]]=desired_mn[k,i]
        # All feature/state values are from the same pre-update instant.
        if keep[k] and k+1<N:
            for post,ids in bypost.items():
                current=float(c.V[post]);pre=c.V[c.pre[ids]]
                saturation=np.clip((pre-c.elo[ids])/c.erange[ids],0.,1.)
                factor=n.geometry_factor[ids] if kind=='MN' else np.ones(len(ids))
                matrix[post][index]=saturation*(c.reversal[ids]-current)*factor
                other=others[post]
                sat=np.clip((c.V[c.pre[other]]-c.elo[other])/c.erange[other],0.,1.)
                existing=np.sum(c.gmax[other]*sat*(c.reversal[other]-current))
                cell=n.names[post]
                nextvalue=pftarget[cell][k+1] if kind in ('PF','BI') else desired_mn[k+1,mn_columns[post]]
                rhs[post][index]=(nextvalue-current)/c.c.time_factor_membrane[post]+c.c.g_m[post]*(current-c.c.V_rest[post])-c.c.i_b[post]-existing
            index+=1
        n.step(u);trace[k]=c.V
    weights={};fits={}
    for post,X in (matrix.items() if os.environ.get('HUNT_BILATERAL_TIED')!='1' else ()):
        y=rhs[post];X=X[:index:5];y=y[:index:5]
        regularizer=.03*np.sqrt(len(y));aug=np.vstack([X,regularizer*np.eye(X.shape[1])]);target=np.r_[y,np.zeros(X.shape[1])]
        solution=lsq_linear(aug,target,bounds=(0.,20. if kind in ('PF','BI') else 10. if kind=='MN' else 2.),method='bvls',tol=1e-8,max_iter=200)
        for i,value in zip(lookup[post],solution.x):weights[n.names[c.pre[i]]+' -> '+n.names[post]]=float(value)
        fits[n.names[post]]=dict(current_rmse=float(np.sqrt(np.mean((X@solution.x-y)**2))),parameters=len(solution.x))
    if os.environ.get('HUNT_BILATERAL_TIED')=='1':
        from bilateral import tied_solutions,assert_tied
        weights,fits=tied_solutions({p:X[:index] for p,X in matrix.items()},{p:y[:index] for p,y in rhs.items()},lookup,n.names,c.pre,c.post,lsq_linear,20. if kind in ('PF','BI') else 10. if kind=='MN' else 2.)
        assert_tied(weights)
    name={'PF':'pf_weights.json','BI':'bi_weights.json','MN':'mn_weights.json','REFLEX':'reflex_weights.json'}[kind]
    (O/name).write_text(json.dumps(dict(block=kind,edges=weights,regression=fits,
        scope='Offline teacher-forced conductance identification; runtime uses native voltages',timing_frozen=True),indent=2))
    np.savez_compressed(O/(kind.lower()+'_teacher.npz'),t=t,voltage=trace,desired_mn=desired_mn,
        pf_cells=list(pftarget),pf_target=np.stack(list(pftarget.values()),axis=1),neuron_names=n.names[:len(c.V)])
    print(kind,'fitted',len(weights),'edges; median current RMSE',np.median([v['current_rmse'] for v in fits.values()]),flush=True)

def validate(target_scale=1.,prefix='neural_validation'):
    z,m,d,n,t,q,vel,a,pftarget=setup();c=n.compiled;u=n.make_inputs();u[n.ports['DRIVE']]=2.5
    mnidx=np.array([n.idx['MN_'+act] for act in z['act_names']]);output=[];pf=[]
    for k in range(len(t)):
        d.qpos[:]=q[k];n.arm_model.project(d);d.qvel[:]=vel[k];d.act[:]=n.expected_activation[k] if hasattr(n,'expected_activation') else a[k];mujoco.mj_forward(m,d);afferents(n,u,d)
        u[n.ports['HUNT_INIT_E_r']]=2. if k<50 else 0.;u[n.ports['HUNT_INIT_F_l']]=2. if k<50 else 0.
        n.step(u);output.append(np.clip(c.V[mnidx]/5.,0.,1.));pf.append([c.V[n.idx[cell]] for cell in pftarget])
    output=np.asarray(output);keep=t>4*1.233
    target=a[keep]*target_scale;error=output[keep]-target
    result=dict(activation_rmse=float(np.sqrt(np.mean(error**2))),peak_target=float(target.max()),
        finite=bool(np.all(np.isfinite(output))),source='Native neural replay with expected physical afferents; no PF or MN voltage clamping',
        target_scale=target_scale,closed_loop_body_test=False,accepted_as_walking=False)
    active=np.max(target,axis=0)>.05
    result['command_rmse']=result['activation_rmse']
    result['comparison']='MN voltage / 5 versus desired excitation command; activation_rmse is a legacy field name, not native activation-state error'
    variance=np.var(target,axis=0);r2=1.-np.mean(error**2,axis=0)/np.maximum(variance,1e-12)
    result['active_muscle_count']=int(active.sum())
    result['active_muscle_rmse']=float(np.sqrt(np.mean(error[:,active]**2))) if active.any() else None
    result['median_active_muscle_r_squared']=float(np.median(r2[active])) if active.any() else None
    result['per_muscle']={str(act):dict(r_squared=float(r2[i]),target_peak=float(target[:,i].max()),predicted_range=float(np.ptp(output[keep,i]))) for i,act in enumerate(z['act_names']) if active[i]}
    (O/(prefix+'.json')).write_text(json.dumps(result,indent=2));np.savez_compressed(O/(prefix+'.npz'),t=t,activation=output,target=a*target_scale,pf=np.asarray(pf))
    print(json.dumps(result),flush=True)
if __name__=='__main__':
    for kind in ('PF','MN','REFLEX'):replay(kind)
    validate()
