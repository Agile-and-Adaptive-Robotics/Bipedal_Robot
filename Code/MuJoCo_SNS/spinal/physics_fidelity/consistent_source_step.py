"""Opt-in whole-step source-equilibrium Newton adapter.

The A403 source law is imported unchanged. Each Newton trial restores the same
native state and old fibers. Candidate tension is constant within the native
implicitfast trial; the converged tension equals the nonlinear source law at
the actual end length. Source velocity dependence is solved implicitly, not
removed. Native passive damping and contact equations are retained.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import numpy as np
import mujoco
from pathlib import Path
import hashlib,json,argparse,sys
from time import perf_counter
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent/'functional_gait2392')]
from scipy.linalg import cho_factor,cho_solve,solve
from thelen_dynamic import DynamicThelenBank,fiber_components


def sparse_dense(values,rownnz,rowadr,colind,nrow,ncol):
    result=np.zeros((nrow,ncol))
    for row in range(nrow):
        start=int(rowadr[row]);count=int(rownnz[row])
        result[row,colind[start:start+count]]=values[start:start+count]
    return result


def tendon_integral(strain,strain0):
    edge=99.*np.exp(3.)/(166.*np.exp(3.)-67.)*strain0
    x=np.maximum(0.,np.minimum(strain,edge))
    value=.33/np.expm1(3.)*(edge/3.*np.expm1(3.*x/edge)-x)
    above=np.maximum(0.,strain-edge)
    return value+.33*above+.5*.67/(strain0-edge)*above**2


def passive_integral(normalized):
    x=np.maximum(0.,normalized-1.)
    return (.6/4.*np.expm1(4.*x/.6)-x)/np.expm1(4.)


def constitutive_branch(bank,length,activation,old=None):
    p=bank.parameters;state=bank.proposed;l=np.hypot(state,p[:,3]);old=bank.state if old is None else old;v=(state-old)*state/(l*bank.dt*p[:,1]);c=.25+.75*np.maximum(.01,activation)
    k=(2.+2./p[:,6])/(p[:,7]-1.);edge=.95*p[:,7];threshold=c*(edge-1.)/(k*(p[:,7]-edge))
    fv=np.select([v/p[:,8]<=-c,v<=0.,v/p[:,8]<=threshold],[0,1,2],default=3)
    strain=(length-state)/p[:,2]-1.;toe=99.*np.exp(3.)/(166.*np.exp(3.)-67.)*p[:,4]
    return fv+4*(l>p[:,1])+8*np.select([strain<=0.,strain<toe],[0,1],default=2)


class ConsistentSourceBackend:
    """Scalar-coordinate,pyramidal-contact prototype; no production hooks."""
    def __init__(self,model,data,parameters,max_iterations=12,tolerance_N=1e-7):
        if model.nq!=model.nv or int(model.opt.integrator)!=int(mujoco.mjtIntegrator.mjINT_IMPLICITFAST):
            raise ValueError('Scalar coordinates and existing implicitfast required')
        if int(model.opt.cone)!=int(mujoco.mjtCone.mjCONE_PYRAMIDAL):
            raise ValueError('Prototype constrained tangent currently covers pyramidal contacts')
        if np.any(model.actuator_forcelimited) or np.any(getattr(model,'jnt_actfrclimited',False)):raise ValueError('Prototype force Jacobian requires unclipped source actuators')
        self.model,self.data=model,data;self.parameters=np.asarray(parameters);self.dt=float(model.opt.timestep)
        if np.any(self.parameters[:,5]!=0.):raise ValueError('This source-preserving prototype expects the original zero explicit fiber-damping parameter')
        self.bank=DynamicThelenBank(parameters,self.dt);self.bank.initialize(data.actuator_length,data.act)
        self.max_iterations=max_iterations;self.tolerance_N=tolerance_N
        self.mask=mujoco.mjtState.mjSTATE_INTEGRATION;self.saved=np.empty(mujoco.mj_stateSize(model,self.mask))
        self.old_state=np.zeros(model.nu);self.activation=np.zeros(model.nu);self.length=np.zeros(model.nu)
        self.predicted_force=np.zeros(model.nu);self.length_tangent=np.zeros(model.nu);self.components=np.zeros((3,model.nu))
        self.M=np.empty((model.nv,model.nv));self.steps=0;self.native_trials=0;self.iterations=0
        self.max_residual=0.;self.max_native_balance=0.;self.max_forward_balance=0.;self.max_energy_identity=0.;self.max_native_work_identity=0.;self.minimum_convex_dissipation=np.inf
        self.max_prediction_error=0.;self.max_prediction_error_F0=0.;self.bound_steps=0
        ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'knee_angle_'+s) for s in ('r','l')]
        self.knees=np.array([int(model.jnt_qposadr[j]) for j in ids if j>=0],int);self.knee_max=data.qpos[self.knees].copy()
        model.actuator_gaintype[:]=mujoco.mjtGain.mjGAIN_FIXED;model.actuator_gainprm[:]=0.
        model.actuator_biastype[:]=mujoco.mjtBias.mjBIAS_AFFINE;model.actuator_biasprm[:]=0.

    def _check(self):
        self.max_residual=max(self.max_residual,float(max(abs(self.bank.error))))
        if np.any(self.bank.bounds):self.bound_steps+=1;raise RuntimeError('Fiber bound gate')
        if max(abs(self.bank.error))>1e-5:raise RuntimeError('Fiber equilibrium gate')
        if not np.isfinite(self.bank.force).all():raise RuntimeError('Nonfinite source tension')
        if np.any(self.bank.force>20*self.parameters[:,0]):raise RuntimeError('Fiber force guard')

    def _trial(self,tension,jacobian=True):
        m,d=self.model,self.data
        mujoco.mj_setState(m,d,self.saved,self.mask)
        m.actuator_biasprm[:]=0.;m.actuator_biasprm[:,0]=-tension
        mujoco.mj_step1(m,d);moment0=d.actuator_moment.copy();mujoco.mj_fullM(m,self.M,d.qM)
        velocity0=d.qvel.copy()
        mujoco.mj_step2(m,d);self.native_trials+=1
        force=d.qfrc_smooth+d.qfrc_constraint
        deriv=sparse_dense(d.qDeriv,d.D_rownnz,d.D_rowadr,d.D_colind,m.nv,m.nv)
        # implicitfast consumes the symmetric lower triangle of qDeriv.
        deriv=np.tril(deriv)+np.tril(deriv,-1).T;B=self.M-self.dt*deriv
        native_error=B@((d.qvel-velocity0)/self.dt)-force
        self.max_native_balance=max(self.max_native_balance,float(max(abs(native_error))))
        self.max_forward_balance=max(self.max_forward_balance,float(max(abs(self.M@d.qacc-force))))
        dv=d.qvel-velocity0;kinetic_change=.5*(d.qvel@self.M@d.qvel-velocity0@self.M@velocity0)
        work=self.dt*d.qvel@(force+deriv@dv);implicit_loss=.5*dv@self.M@dv
        self.max_native_work_identity=max(self.max_native_work_identity,float(abs(kinetic_change-work+implicit_loss)))
        self.last_native_matrix=B.copy();self.last_native_derivative=deriv.copy();self.last_native_force=force.copy();self.last_moment0=moment0.copy()
        if jacobian:
            J=sparse_dense(d.efc_J,d.efc_J_rownnz,d.efc_J_rowadr,d.efc_J_colind,d.nefc,m.nv) if mujoco.mj_isSparse(m) else d.efc_J.reshape(d.nefc,m.nv).copy()
            states=np.asarray(d.efc_state[:d.nefc]);weights=np.where(states==int(mujoco.mjtConstraintState.mjCNSTRSTATE_QUADRATIC),d.efc_D[:d.nefc],0.)
            if np.any(states==int(mujoco.mjtConstraintState.mjCNSTRSTATE_CONE)):raise RuntimeError('Unsupported elliptic constraint Hessian')
            A=self.M+J.T@(weights[:,None]*J)
            response=cho_solve(cho_factor(B,lower=True,check_finite=False),self.M@cho_solve(cho_factor(A,lower=True,check_finite=False),moment0.T,check_finite=False),check_finite=False)
        mujoco.mj_forward(m,d);self.bank.state[:]=self.old_state
        self.bank.evaluate(d.actuator_length,self.activation);self._check()
        residual=tension-self.bank.force
        matrix=np.eye(m.nu)+self.dt**2*self.bank.tangent[:,None]*(d.actuator_moment@response) if jacobian else None
        return residual,matrix

    def step(self):
        old_fiber=self.bank.state.copy();old_force=self.bank.force.copy();old_bias=self.model.actuator_biasprm.copy()
        mujoco.mj_getState(self.model,self.data,self.saved,self.mask)
        try:self._step_impl()
        except Exception:
            mujoco.mj_setState(self.model,self.data,self.saved,self.mask);self.model.actuator_biasprm[:]=old_bias
            self.bank.state[:]=old_fiber;self.bank.force[:]=old_force;mujoco.mj_forward(self.model,self.data)
            raise

    def _step_impl(self):
        m,d=self.model,self.data;self.activation[:]=d.act;self.old_state[:]=self.bank.state
        old_length=d.actuator_length.copy();mujoco.mj_getState(m,d,self.saved,self.mask)
        self.bank.evaluate(old_length+self.dt*d.actuator_velocity,self.activation);self._check();tension=self.bank.force.copy()
        iterations=0;converged=False
        for iterations in range(self.max_iterations):
            residual,matrix=self._trial(tension);error=float(max(abs(residual)))
            if error<=self.tolerance_N:converged=True;break
            change=solve(matrix,-residual,assume_a='gen',check_finite=False);accepted=False
            for power in range(8):
                candidate=np.maximum(0.,tension+(2.**-power)*change)
                proposed,_=self._trial(candidate,False)
                if max(abs(proposed))<error:
                    tension=candidate;accepted=True
                    if max(abs(proposed))<=self.tolerance_N:
                        residual=proposed;error=float(max(abs(proposed)));converged=True
                    break
            if not accepted:raise RuntimeError('Whole-step source Newton line search failed')
            if converged:break
        else:raise RuntimeError('Whole-step source Newton did not converge')
        # Last candidate is one native physical step; no repeated fiber advance.
        self.max_prediction_error=max(self.max_prediction_error,error);self.max_prediction_error_F0=max(self.max_prediction_error_F0,float(max(abs(residual)/self.parameters[:,0])))
        self.bank.state[:]=self.bank.proposed;self.length[:]=d.actuator_length;self.predicted_force[:]=tension;self.length_tangent[:]=self.bank.tangent
        fiber_components(self.bank.state,self.old_state,self.activation,self.parameters,self.dt,self.components)
        self.energy_check(old_length)
        m.actuator_biasprm[:]=0.;m.actuator_biasprm[:,0]=-self.bank.force;mujoco.mj_forward(m,d)
        self.steps+=1;self.iterations+=iterations+1;self.knee_max=np.maximum(self.knee_max,d.qpos[self.knees])
        if np.any(self.knee_max>0.):raise RuntimeError('Strict zero-hyperextension gate')
        if not np.isfinite(d.qpos).all() or max(abs(d.qvel))>300.:raise RuntimeError('Nonfinite/velocity guard')

    def energy_check(self,old_length):
        p=self.parameters;s0=self.old_state;s1=self.bank.state;f0=p[:,0];l0=p[:,1];slack=p[:,2]
        lf0=np.hypot(s0,p[:,3]);lf1=np.hypot(s1,p[:,3]);lt0=old_length-s0;lt1=self.data.actuator_length-s1
        du_t=f0*slack*(tendon_integral(lt1/slack-1.,p[:,4])-tendon_integral(lt0/slack-1.,p[:,4]))
        du_p=f0*l0*(passive_integral(lf1/l0)-passive_integral(lf0/l0))
        co=s1/lf1;active=self.components[0]/co;passive=self.components[1]/co
        work=self.bank.force*(self.data.actuator_length-old_length);active_work=active*(lf1-lf0)
        diss_t=self.bank.force*(lt1-lt0)-du_t;diss_p=passive*(lf1-lf0)-du_p;diss_geometry=(active+passive)*(co*(s1-s0)-(lf1-lf0))
        defect=work-du_t-du_p-active_work-diss_t-diss_p-diss_geometry
        self.max_energy_identity=max(self.max_energy_identity,float(max(abs(defect))));self.minimum_convex_dissipation=min(self.minimum_convex_dissipation,float(min(diss_t.min(),diss_p.min(),diss_geometry.min())))

    def diagnostics(self):
        return dict(variant='whole_step_source_equilibrium_newton_v1',steps=self.steps,native_trials=self.native_trials,newton_iterations=self.iterations,
            max_equilibrium_residual_N=self.max_residual,max_whole_step_force_residual_N=self.max_prediction_error,max_native_implicit_balance_residual=self.max_native_balance,max_native_forward_balance_residual=self.max_forward_balance,
            max_constitutive_energy_identity_error_J=self.max_energy_identity,max_frozen_mass_native_work_identity_error_J=self.max_native_work_identity,minimum_convex_discrete_dissipation_J=self.minimum_convex_dissipation,maximum_knee_extension_deg=np.rad2deg(self.knee_max).tolist(),
            scope='Original source curves,activation timing,contacts,native damping,and implicitfast retained. Newton source tension is self-consistent with actual whole-body end length; original backend/default unchanged.')


def capture_crossings(model_folder,start_folder,target_folder,out):
    from predictive_source_recruitment import AugmentedTransition
    from source_force_comparison import load_recruitment
    if (out/'crossing_snapshots.npz').exists():raise ValueError('Preserve captured states')
    a=load_recruitment(start_folder);target=load_recruitment(target_folder)
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));m.opt.timestep=.001
    r=AugmentedTransition(m,a['parameters'],.05);m,d,b=r.model,r.data,r.backend
    d.qpos[:]=a['qpos'];d.act[:]=a['activation'];d.ctrl[:]=a['activation'];d.xfrc_applied[m.body('pelvis').id,2]=.6*float(m.body_mass.sum()*9.81)
    mujoco.mj_forward(m,d);b.bank.initialize(d.actuator_length,d.act);m.actuator_biasprm[:,0]=-b.bank.force;mujoco.mj_forward(m,d);r.commit_snapshot();r.reset()
    selected={7275,7281,7303,7307,7308};rows=[];checksum=hashlib.sha256();started=perf_counter()
    for step in range(1,7313):
        phase=np.clip((d.time-.5)/1.812,0.,1.);blend=phase**3*(10.+phase*(-15.+6.*phase));command=a['activation']+blend*(target['activation']-a['activation']);r.excitation+=r.beta*(command-r.excitation);d.ctrl[:]=r.excitation
        if step in selected:
            state=np.empty(mujoco.mj_stateSize(m,r.mask));mujoco.mj_getState(m,d,state,r.mask)
            row=dict(integration_state=state,bias=m.actuator_biasprm.copy(),old_fiber=b.bank.state.copy(),old_force=b.bank.force.copy(),moment0=d.actuator_moment.copy(),velocity0=d.qvel.copy(),old_length=d.actuator_length.copy())
        b.step();checksum.update(np.r_[d.qpos,d.qvel,d.act,b.bank.state,b.bank.force,d.ctrl].tobytes())
        if step in selected:
            derivative=sparse_dense(d.qDeriv,d.D_rownnz,d.D_rowadr,d.D_colind,m.nv,m.nv);derivative=np.tril(derivative)+np.tril(derivative,-1).T
            row.update(end_time=d.time,end_qpos=d.qpos.copy(),end_qvel=d.qvel.copy(),end_activation=d.act.copy(),end_fiber=b.bank.state.copy(),end_force=b.bank.force.copy(),original_predicted_force=b.predicted_force.copy(),original_derivative=derivative);rows.append(row)
    digest=checksum.hexdigest()
    expected='1682c36c4f6c9581ef88642c4a4db90a4544fb3666627993ec874ea0f1c3aba1'
    if digest!=expected:raise RuntimeError('Frozen original input/trajectory replay changed')
    out.mkdir(parents=True,exist_ok=True);np.savez_compressed(out/'crossing_snapshots.npz',state_spec=int(r.mask),parameters=a['parameters'],muscles=a['muscles'],**{key:np.array([row[key] for row in rows]) for key in rows[0]})
    report=dict(original_state_force_sha256=digest,expected_sha256=expected,exact_original_replay=True,steps=7312,wall_s=perf_counter()-started,source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model=str(model_folder),start=str(start_folder),target=str(target_folder),scope='Five pre-step snapshots with identical old physical/fiber/input state; original A403 trajectory retained.')
    (out/'CAPTURE.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


def load_context(model_folder,archive,index,consistent=True):
    from thelen_dynamic import ConstitutiveBackend
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));m.opt.timestep=.001;d=mujoco.MjData(m)
    mujoco.mj_setState(m,d,archive['integration_state'][index],int(archive['state_spec']));mujoco.mj_forward(m,d)
    b=(ConsistentSourceBackend if consistent else ConstitutiveBackend)(m,d,archive['parameters'])
    mujoco.mj_setState(m,d,archive['integration_state'][index],int(archive['state_spec']));m.actuator_biasprm[:]=archive['bias'][index]
    b.bank.state[:]=archive['old_fiber'][index];b.bank.force[:]=archive['old_force'][index];b.bank.bounds[:]=0;mujoco.mj_forward(m,d)
    return m,d,b


def validate_contexts(model_folder,folder,out):
    if (out/'CONTEXT_VALIDATION.json').exists():raise ValueError('Preserve existing frozen-context validation')
    z=np.load(folder/'crossing_snapshots.npz');rows=[];started=perf_counter()
    for index,time in enumerate(z['end_time']):
        mo,do,original=load_context(model_folder,z,index,False);original.step()
        parity=max(float(np.max(abs(value-z[key][index]))) for key,value in dict(end_qpos=do.qpos,end_qvel=do.qvel,end_fiber=original.bank.state,end_force=original.bank.force).items())
        m,d,b=load_context(model_folder,z,index,True);mujoco.mj_getState(m,d,b.saved,b.mask);b.activation[:]=d.act;b.old_state[:]=b.bank.state
        b.bank.evaluate(d.actuator_length+b.dt*d.actuator_velocity,b.activation);guess=b.bank.force.copy();base,J=b._trial(guess)
        direction=np.sin(np.arange(m.nu)+1.);direction/=np.linalg.norm(direction);checks=[]
        for epsilon in (.001,.0001,.00001):
            rp,_=b._trial(guess+epsilon*direction,False);rm,_=b._trial(guess-epsilon*direction,False);measured=(rp-rm)/(2*epsilon);predicted=J@direction
            checks.append(dict(force_step_N=epsilon,relative_error=float(np.linalg.norm(measured-predicted)/np.linalg.norm(measured)),absolute_error=float(max(abs(measured-predicted)))))
        # Each law derivative is queried from unchanged old fibers and input.
        lengths=do.actuator_length.copy();alpha=b.activation.copy();bank=DynamicThelenBank(z['parameters'],.001);bank.state[:]=z['old_fiber'][index];bank.evaluate(lengths,alpha);gain=bank.tangent.copy();base_force=bank.force.copy();base_branch=constitutive_branch(bank,lengths,alpha);central=[]
        for epsilon in (1e-8,1e-9):
            bank.evaluate(lengths+epsilon,alpha);plus=bank.force.copy();plus_branch=constitutive_branch(bank,lengths+epsilon,alpha)
            bank.evaluate(lengths-epsilon,alpha);minus=bank.force.copy();minus_branch=constitutive_branch(bank,lengths-epsilon,alpha);fd=(plus-minus)/(2*epsilon)
            same_plus=plus_branch==base_branch;same_minus=minus_branch==base_branch;one_sided=np.where(same_plus&same_minus,fd,np.where(same_plus,(plus-base_force)/epsilon,np.where(same_minus,(base_force-minus)/epsilon,np.nan)))
            central.append(dict(length_step_m=epsilon,maximum_central_relative_error=float(max(abs(fd-gain)/np.maximum(1.,abs(gain)))),maximum_branch_consistent_relative_error=float(np.nanmax(abs(one_sided-gain)/np.maximum(1.,abs(gain)))),unresolved_channels=int(np.isnan(one_sided).sum())))
        # Solve from exactly the same frozen old state,not a reset end pose.
        m,d,b=load_context(model_folder,z,index,True);b.step()
        source_defect=z['moment0'][index].T@(z['end_force'][index]-z['original_predicted_force'][index])+(z['original_derivative'][index]-b.last_native_derivative)@(z['end_qvel'][index]-z['velocity0'][index])
        full_source_D=-b.dt*z['moment0'][index].T@(original.length_tangent[:,None]*z['moment0'][index]);native_source_D=z['original_derivative'][index]-b.last_native_derivative
        rows.append(dict(time_s=float(time),original_replay_error_max=parity,native_full_force_Jacobian_FD=checks,source_force_length_tangent_FD=central,consistent=b.diagnostics(),native_source_D_relative_omission=float(np.linalg.norm(full_source_D-native_source_D)/np.linalg.norm(full_source_D)),original_whole_step_source_balance_defect_max=float(max(abs(source_defect))),activation_timing_error_max=float(max(abs(d.act-z['end_activation'][index]))),end_qpos_change_max=float(max(abs(d.qpos-z['end_qpos'][index]))),mtp_acceleration_rad_s2=[float(d.qacc[m.joint('mtp_angle_'+s).dofadr[0]]) for s in ('r','l')]))
    report=dict(scope='Frozen identical native state,fibers,and inputs; original versus whole-step Newton integration. No parameter/gain search.',rows=rows,wall_s=perf_counter()-started,source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.mkdir(parents=True,exist_ok=True);(out/'CONTEXT_VALIDATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


def replay_reference(model_folder,start_folder,target_folder,spec_path,out,backend='consistent',dt=.001):
    """One frozen C2 input,not a gain/horizon search or controller promotion."""
    from thelen_dynamic import ConstitutiveBackend
    from source_force_comparison import load_recruitment
    from predictive_source_recruitment import NativeContactKinematics,AugmentedTransition
    from mechanics import CoupledArms
    from contact_standing_solver import _external
    if (out/'REPLAY.json').exists():raise ValueError('Preserve existing integration replay')
    a=load_recruitment(start_folder);target=load_recruitment(target_folder);p=a['parameters']
    if not np.array_equal(p,target['parameters']) or not a['force_channel_enabled'].all() or not target['force_channel_enabled'].all():raise ValueError('Unchanged original92 source policy required')
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));m.opt.timestep=dt
    runner=AugmentedTransition(m,p,.05,backend_class=ConsistentSourceBackend if backend=='consistent' else None);m,d,b=runner.model,runner.data,runner.backend
    d.qpos[:]=a['qpos'];d.act[:]=a['activation'];d.ctrl[:]=a['activation'];bw=float(m.body_mass.sum()*9.81);d.xfrc_applied[m.body('pelvis').id,2]=.6*bw;mujoco.mj_forward(m,d)
    b.bank.initialize(d.actuator_length,d.act);m.actuator_biasprm[:,0]=-b.bank.force;mujoco.mj_forward(m,d);runner.commit_snapshot();runner.reset()
    filtered=runner.excitation;beta=runner.beta;ar=CoupledArms(m);contact=NativeContactKinematics(m,d);mtp=np.array([int(m.joint('mtp_angle_'+s).qposadr[0]) for s in ('r','l')])
    spec=json.loads(spec_path.read_text());address=np.array([v['qpos_address'] for v in spec['joints'].values()]);lo=np.array([v['native_bounds_rad'][0] for v in spec['joints'].values()]);hi=np.array([v['native_bounds_rad'][1] for v in spec['joints'].values()])
    ref=np.load(spec['source_target']);speed=np.max(abs(ref['qvel'][:,address]),axis=0);acceleration=np.max(abs(ref['qacc'][:,address]),axis=0);geometry=[]
    for g in range(m.ngeom):
        name=m.body(int(m.geom_bodyid[g])).name
        if m.geom_type[g]==mujoco.mjtGeom.mjGEOM_MESH and name in ('calcn_r','toes_r','calcn_l','toes_l'):
            mesh=int(m.geom_dataid[g]);first=int(m.mesh_vertadr[mesh]);number=int(m.mesh_vertnum[mesh]);geometry.append((0 if name.endswith('_r') else 1,g,m.mesh_vert[first:first+number].copy()))
    n=round(7.312/dt);phase_trace=[];branch_trace=[];sparse=[];late=[];checksum=hashlib.sha256();guard='';maxslip=0.;maxroot=0.;maxrot=0.;mingap=np.full(2,np.inf);started=perf_counter()
    out.mkdir(parents=True,exist_ok=True)
    for i in range(n):
        phase=np.clip((d.time-.5)/1.812,0.,1.);blend=phase**3*(10.+phase*(-15.+6.*phase));u=a['activation']+blend*(target['activation']-a['activation']);filtered+=beta*(u-filtered);d.ctrl[:]=filtered
        try:b.step()
        except RuntimeError as error:guard=str(error);break
        measure=contact.measure(True);normal=measure['normal_N']/bw;gap=np.full(2,np.inf)
        for side,g,vertices in geometry:gap[side]=min(gap[side],float(min(vertices@d.geom_xmat[g].reshape(3,3)[2]+d.geom_xpos[g,2])))
        mingap=np.minimum(mingap,gap);maxslip=max(maxslip,float(max(measure['slip_mps'])));maxroot=max(maxroot,float(max(abs(d.qpos[:3]-a['qpos'][:3]))));maxrot=max(maxrot,float(max(abs(d.qpos[3:6]-a['qpos'][3:6]))))
        phase_trace.append(np.r_[d.time,d.qpos[mtp],d.qvel[mtp],d.qacc[mtp],measure['normal_by_body_N'],measure['signature_u32']]);branch_trace.append(constitutive_branch(b.bank,d.actuator_length,b.activation,b.old_state))
        checksum.update(np.r_[d.qpos,d.qvel,d.act,b.bank.state,b.bank.force,d.ctrl].tobytes())
        if i%max(1,round(.01/dt))==0:sparse.append(np.r_[d.time,d.qpos,d.qvel,d.act,b.bank.state,d.ctrl,normal,gap])
        if d.time>=6.812-1e-10:
            length=np.hypot(b.bank.state,p[:,3]);velocity=(b.bank.state-b.old_state)*b.bank.state/(length*dt*p[:,1]);late.append((d.qpos.copy(),d.qvel.copy(),b.bank.force.copy(),velocity,normal.copy()))
        if np.any(d.qpos[address]<lo-1e-10) or np.any(d.qpos[address]>hi+1e-10):guard='Frozen source ROM'
        elif max(abs(d.qvel[address])/speed)>1.+1e-10 or max(abs(d.qacc[address])/acceleration)>1.+1e-10:guard='Frozen source rate envelope'
        elif maxroot>.005 or maxrot>np.deg2rad(2.):guard='Root envelope'
        elif maxslip>.005 or mingap.min()<-.002:guard='Contact slip/mesh envelope'
        elif measure['maximum_friction_cone_violation_N']>1e-6 or measure['minimum_normal_N']< -1e-8:guard='Native contact cone'
        other=[int(m.joint(int(j)).qposadr[0]) for j in ar.jids if m.joint(int(j)).name not in spec['joints']]
        if max(abs(d.qpos[other]))>np.deg2rad(10.)+1e-10:guard=guard or 'Existing nonsagittal bounds'
        if guard:break
    probe=mujoco.MjData(m);probe.qpos[:]=d.qpos;probe.act[:]=d.act;probe.xfrc_applied[:]=d.xfrc_applied;mujoco.mj_inverse(m,probe);required=probe.qfrc_inverse-_external(m,probe)
    static=DynamicThelenBank(p,dt);static.initialize(probe.actuator_length,d.act);raw_static=-probe.actuator_moment.T@static.force-required;raw_dynamic=-probe.actuator_moment.T@b.bank.force-required
    pt=np.array(phase_trace);window=pt[:,0]>=6.812-1e-10;targetnorm=np.array(json.loads((target_folder/'MANIFEST.json').read_text())['stationary_load_split_reference']['normal_target_BW']);digest=checksum.hexdigest()
    metrics=dict(duration_s=float(d.time),guard=guard,completed_without_live_guard=bool(not guard and i==n-1),wall_s=perf_counter()-started,maximum_loaded_slip_mps=maxslip,maximum_root_translation_m=maxroot,maximum_root_rotation_deg=float(np.rad2deg(maxrot)),minimum_mesh_gap_m=mingap.tolist(),maximum_knee_extension_deg=np.rad2deg(b.knee_max).tolist(),
        late_MTP_acceleration_abs_max_rad_s2=np.max(abs(pt[window,5:7]),axis=0).tolist() if window.any() else None,late_MTP_velocity_sign_changes=[int(np.sum(pt[window,3+j][1:]*pt[window,3+j][:-1]<0.)) for j in range(2)] if window.any() else None,
        late_joint_speed_max_deg_s=max((float(np.rad2deg(max(abs(v[1][ar.rows])))) for v in late),default=None),late_fiber_speed_l0_s=max((float(max(abs(v[3]))) for v in late),default=None),late_normal_error_BW=max((float(max(abs(v[4]-targetnorm))) for v in late),default=None),
        final_full_static_residual_max=float(max(abs(raw_static))),final_full_dynamic_force_static_residual_max=float(max(abs(raw_dynamic))),final_joint_pose_error_deg=float(np.rad2deg(max(abs(d.qpos[ar.addresses]-target['qpos'][ar.addresses])))),backend=b.diagnostics())
    report=dict(accepted_motion=False,backend=backend,timestep_s=dt,support_fraction=.6,model=str(model_folder),start=str(start_folder),target=str(target_folder),source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),physical_state_force_sha256=digest,original_1ms_parity=bool(digest=='1682c36c4f6c9581ef88642c4a4db90a4544fb3666627993ec874ea0f1c3aba1') if backend=='original' and dt==.001 else None,metrics=metrics,scope='Same frozen C2 input,source parameters,and native physics; opt-in integration comparison only.')
    (out/'REPLAY.json').write_text(json.dumps(report,indent=2));np.savez_compressed(out/'replay_sparse.npz',phase=pt,corrected_branch=np.array(branch_trace,dtype=np.uint8),samples=np.array(sparse),muscles=a['muscles'],nq=m.nq,nu=m.nu)
    print(json.dumps(report,indent=2))


def decompose_crossings(model_folder,folder,out):
    """Read-only full-M acceleration decomposition,without follower projection."""
    from contact_standing_solver import _external
    z=np.load(folder/'crossing_snapshots.npz');rows=[]
    for index in range(len(z['end_time'])):
        m,d,b=load_context(model_folder,z,index,False)
        for when in ('before','after'):
            if when=='after':b.step()
            # Before-step fiber velocity is not stored; exact total muscle tension
            # is used for both states rather than inventing active/passive splits.
            moment=d.actuator_moment.copy();M=np.empty((m.nv,m.nv));mujoco.mj_fullM(m,M,d.qM)
            probe=mujoco.MjData(m);probe.qpos[:]=d.qpos;probe.qvel[:]=0.;mujoco.mj_forward(m,probe)
            J=sparse_dense(d.efc_J,d.efc_J_rownnz,d.efc_J_rowadr,d.efc_J_colind,d.nefc,m.nv) if mujoco.mj_isSparse(m) else d.efc_J.reshape(d.nefc,m.nv).copy()
            kinds=np.asarray(d.efc_type[:d.nefc]);force=np.asarray(d.efc_force[:d.nefc]);parts=dict(muscles=d.qfrc_actuator.copy(),native_position_passive=probe.qfrc_passive.copy(),native_velocity_passive=d.qfrc_passive-probe.qfrc_passive,gravity=-probe.qfrc_bias.copy(),coriolis=-(d.qfrc_bias-probe.qfrc_bias),external=_external(m,d),equality=J.T@(force*(kinds==0)),contact=J.T@(force*(kinds>=5)),friction=J.T@(force*((kinds==1)|(kinds==2))),limits=J.T@(force*((kinds==3)|(kinds==4))))
            keys=list(parts);acc=solve(M,np.column_stack([parts[k] for k in keys]),assume_a='pos',check_finite=False);residual=acc.sum(axis=1)-d.qacc
            muscle_acc=solve(M,-moment.T*b.bank.force[None,:],assume_a='pos',check_finite=False)
            sides={}
            for s in ('r','l'):
                joint=m.joint('mtp_angle_'+s);row=int(joint.dofadr[0]);order=np.argsort(abs(muscle_acc[row]))[::-1][:8]
                sides[s]=dict(actual_acceleration_rad_s2=float(d.qacc[row]),generalized_force_Nm={k:float(v[row]) for k,v in parts.items()},whole_body_acceleration_components_rad_s2={k:float(acc[row,j]) for j,k in enumerate(keys)},largest_individual_muscle_acceleration_components=[dict(muscle=str(z['muscles'][j]),contribution_rad_s2=float(muscle_acc[row,j]),tension_N=float(b.bank.force[j])) for j in order])
            rows.append(dict(time_s=float(d.time),position=when,full105_acceleration_reconstruction_error=float(max(abs(residual))),sides=sides))
    result=dict(scope='Linear decomposition of actual instantaneous native forward acceleration using the full mass matrix and actual constraint forces. Large terms can cancel; individual size is not causal attribution. No trajectory/physics change.',rows=rows,source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    out.mkdir(parents=True,exist_ok=True);(out/'MTP_ACCELERATION_COMPONENTS.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--model',required=True,type=Path);parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--capture-start',type=Path);parser.add_argument('--capture-target',type=Path);parser.add_argument('--validate-contexts',type=Path)
    parser.add_argument('--decompose-contexts',type=Path)
    parser.add_argument('--replay-start',type=Path);parser.add_argument('--replay-target',type=Path);parser.add_argument('--spec',type=Path);parser.add_argument('--backend',choices=('original','consistent'),default='original');parser.add_argument('--physics-dt',type=float,default=.001)
    args=parser.parse_args()
    if args.decompose_contexts:decompose_crossings(args.model,args.decompose_contexts,args.out)
    elif args.replay_start:replay_reference(args.model,args.replay_start,args.replay_target,args.spec,args.out,args.backend,args.physics_dt)
    elif args.capture_start:capture_crossings(args.model,args.capture_start,args.capture_target,args.out)
    elif args.validate_contexts:validate_contexts(args.model,args.validate_contexts,args.out)
    else:parser.error('Select capture or frozen-context validation')
