"""Offline W2L cycle references: geometry before dynamics before neural fitting.

No reference is injected into an accepted SNS run. The first gate registers
stance material points to the moving belt and enforces original joint bounds.
Toe gravity initialization is only an IK seed; passive dynamics must follow.
"""
import json,itertools,time
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
from import_virtual_walker import mujoco


def smooth(x):
    x=np.clip(x,0,1);return x*x*x*(10+x*(-15+6*x))


def register_measured_cycle(root):
    """Register whole measured strides; never confuse near closure with a new gait."""
    from native_kernel import NativeController,polynomial
    root=Path(root);ledger=json.loads((root/'INTERLOCK_GATES.json').read_text());row=next(c for c in ledger['cases'] if Path(c['path']).name=='pace_0p2_supported10_v1');base=Path(row['path']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());z=np.load(base/'trajectory.npz');qtrace=z['q'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));dt=rec['result']['physics_dt'];strikes=[v['touchdown_s'] for v in row['landings'] if v['side']=='L' and 2<=v['touchdown_s']<=8];legs=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n+'_'+s) for s in 'LR' for n in ('hip','knee','ankle','toe')];qa=m.jnt_qposadr[legs];vel=np.zeros((len(qtrace),m.nv))
    for k in range(1,len(qtrace)):mujoco.mj_differentiatePos(m,vel[k],dt,qtrace[k-1],qtrace[k])
    cases=[]
    for a,b in zip(strikes[:-1],strikes[1:]):
        i=int(round(a/dt));j=int(round(b/dt));q0=qtrace[i];q1=qtrace[j];delta=np.zeros(m.nv);mujoco.mj_differentiatePos(m,delta,1.,q0,q1);hinge=np.degrees(q1[qa]-q0[qa]);angle=np.degrees(np.linalg.norm(delta[3:6]));dyz=float(np.linalg.norm((q1-q0)[1:3]));score=float(np.mean(hinge**2)+angle**2+(1000*dyz)**2);ranges=[]
        for jid,adr in zip(legs,qa):
            if m.jnt_limited[jid]:ranges.append(float(max(0.,m.jnt_range[jid,0]-qtrace[i:j+1,adr].min(),qtrace[i:j+1,adr].max()-m.jnt_range[jid,1])))
        cases.append(dict(start_s=a,end_s=b,period_s=b-a,indices=[i,j],selection_score=score,world_translation_x_m=float(q1[0]-q0[0]),belt_relative_progress_m=float(q1[0]-q0[0]+.2*(b-a)),root_yz_closure_m=dyz,body_orientation_closure_deg=float(angle),joint_angle_closure_deg=hinge.tolist(),maximum_nominal_joint_limit_overshoot_deg=float(np.degrees(max(ranges))),root_velocity_closure_m_s=(vel[j,:3]-vel[i,:3]).tolist(),maximum_joint_velocity_closure_deg_s=float(np.max(np.abs(np.degrees(vel[j,m.jnt_dofadr[legs]]-vel[i,m.jnt_dofadr[legs]]))))))
    assert cases;chosen=min(cases,key=lambda c:c['selection_score']);i,j=chosen['indices'];bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];core=NativeController(rec['source_configuration']['source'],bank);outputs=[(pre,k,coeff) for pre,k,coeff in core.outputs if bank[k]['type']=='LinearHillMuscle'];vm=np.column_stack([polynomial(coeff,z['neuro'][i:j+1,pre]*.001) for pre,k,coeff in outputs]);folder=root/'periodic_reference';folder.mkdir(exist_ok=True);keep={k:z[k][i:j+1] for k in ('q','angles','foot_force_world','foot_cop_world','foot_clearance','com_position','com_velocity','rail_force_world','guide_force_world','harness_force')};np.savez_compressed(folder/'MEASURED_CYCLE_SEED.npz',t=z['t'][i:j+1]-z['t'][i],qvel_reconstructed=vel[i:j+1],motor_voltage_V=vm,motor_names=[core.names[p] for p,k,c in outputs],muscle_names=[bank[k]['name'] for p,k,c in outputs],muscle_tension_N=z['forces'][i:j+1][:,[k for p,k,c in outputs]],**keep)
    report=dict(machine='easteregg2',baseline=str(base),source_sha256=rec['source_sha256'],candidates=cases,selected=chosen,selection_rule='All adjacent actual left force-qualified strikes between2and8s. Minimize mean squared leg-angle closure(deg), body-angle closure(deg), and root yz closure(mm); no outcome cherry-pick or parameter fit.',accepted_periodic_reference=False,scope='Complete measured stride under the preserved real SNS/Hill/MuJoCo model, not a new result or a force-allocator success. Coordinates translate in worldX; relative-belt progression separately reported. Velocities reconstructed from exact2ms pose differences, not directly logged. Source soft joint-limit overshoot is explicitly measured and not re-labelled nominal ROM. Closure/acceleration/contact/Hill reachability of any fitted periodic or warped reference still must pass; no quarter-cycle promotion.');(folder/'MEASURED_CYCLE_REGISTRATION.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def validate_measured_dynamics(root):
    """Complete-cycle one-step reconstruction and inverse-dynamics control gate."""
    root=Path(root);folder=root/'periodic_reference';reg=json.loads((folder/'MEASURED_CYCLE_REGISTRATION.json').read_text());base=Path(reg['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];z=np.load(base/'trajectory.npz');q=z['q'];F=z['forces'];seed=np.load(folder/'MEASURED_CYCLE_SEED.npz');v=seed['qvel_reconstructed'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));act=[(i,mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,b['name'])) for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];motor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide');i,j=reg['selected']['indices'];d=mujoco.MjData(m);inv=mujoco.MjData(m);rows=[];started=time.perf_counter()
    accelerations=[]
    for k in range(i,j):
        mujoco.mj_resetData(m,d);d.qpos[:]=q[k];d.qvel[:]=v[k-i]
        for b,a in act:d.ctrl[a]=F[k+1,b]
        if motor>=0:d.ctrl[motor]=.2
        mujoco.mj_step1(m,d);mujoco.mj_step2(m,d);qacc=d.qacc.copy();accelerations.append(qacc);applied=d.qfrc_actuator.copy();err=np.zeros(m.nv);mujoco.mj_differentiatePos(m,err,1.,q[k+1],d.qpos);mujoco.mj_resetData(m,inv);inv.qpos[:]=q[k];inv.qvel[:]=v[k-i];inv.qacc[:]=qacc;mujoco.mj_inverse(m,inv);rows.append([float(abs(err).max()),float(abs(d.qvel-v[k-i+1]).max()),float(abs(inv.qfrc_inverse-applied).max())])
    np.savez_compressed(folder/'MEASURED_ACCELERATION_BASIS.npz',qacc_forward=np.array(accelerations));a=np.array(rows);report=dict(machine='easteregg2',source=str(base),source_period_s=reg['selected']['period_s'],physics_steps=len(rows),physics_dt=float(m.opt.timestep),integrator_option=int(m.opt.integrator),max_position_tangent_error=float(a[:,0].max()),max_velocity_error=float(a[:,1].max()),max_inverse_applied_force_error=float(a[:,2].max()),median_inverse_applied_force_error=float(np.median(a[:,2])),warmed_wall_s=time.perf_counter()-started,scope='Every2ms interval in the selected complete measured cycle. Replay actual next Hill tensions through original actuators, same split-step MuJoCo API, original contacts/passive toes and apparatus; reconstruct velocity from prior pose. Inverse-force control uses the acceleration actually computed by that forward step, not a smoothed finite-difference surrogate. No neural/physics parameter change or new gait success.');report['position_reconstruction_pass']=bool(a[:,0].max()<1e-7);report['inverse_control_pass']=bool(a[:,2].max()<1e-4);(folder/'MEASURED_DYNAMICS_CONTROL_GATE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def close_and_audit(root):
    """Construct closure corrections; derived acceleration never labelled measured."""
    from scipy.optimize import lsq_linear
    root=Path(root);folder=root/'periodic_reference';reg=json.loads((folder/'MEASURED_CYCLE_REGISTRATION.json').read_text());base=Path(reg['baseline']);seed=np.load(folder/'MEASURED_CYCLE_SEED.npz');q0=seed['q'];v0=seed['qvel_reconstructed'];acc=np.load(folder/'MEASURED_ACCELERATION_BASIS.npz')['qacc_forward'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);N=len(acc);dt=float(m.opt.timestep);u=np.arange(N+1)/N;q=q0.copy();drift=np.zeros(m.nq);drift[0]=q0[-1,0]-q0[0,0]
    for j in range(m.njnt):
        if m.jnt_type[j]==mujoco.mjtJoint.mjJNT_SLIDE:
            a=m.jnt_qposadr[j];drift[a]=q0[-1,a]-q0[0,a]
    error=q0[-1]-q0[0]-drift
    for a in list(range(3))+list(range(7,m.nq)):q[:,a]-=smooth(u)*error[a]
    rots=Rotation.from_quat(np.roll(q0[:,3:7],-1,axis=1));end=(rots[0].inv()*rots[-1]).as_rotvec()
    q[:,3:7]=np.roll((rots*Rotation.from_rotvec(-smooth(u)[:,None]*end)).as_quat(),1,axis=1)
    v=np.empty((N+1,m.nv));previous=q[N-1]-drift;mujoco.mj_differentiatePos(m,v[0],dt,previous,q[0])
    for k in range(1,N+1):mujoco.mj_differentiatePos(m,v[k],dt,q[k-1],q[k])
    derived=acc+(np.diff(v,axis=0)-np.diff(v0,axis=0))/dt;bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];tids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act];driven=[m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n+'_'+s)] for s in 'LR' for n in ('hip','knee','ankle')];toes=[m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'toe_'+s)] for s in 'LR'];required=[];fitres=[];feet=[]
    for k in range(N):
        d.qpos[:]=q[k];d.qvel[:]=v[k];d.qacc[:]=derived[k];mujoco.mj_inverse(m,d);required.append(d.qfrc_inverse.copy());A=-d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tids,driven)].T;fit=lsq_linear(A,d.qfrc_inverse[driven],bounds=(np.zeros(len(act)),p[act,5]),method='bvls',tol=1e-10);fitres.append(float(np.max(np.abs(A@fit.x-d.qfrc_inverse[driven]))));feet.append([float(d.xpos[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'foot_'+s),1]) for s in 'LR'])
    req=np.array(required);bw=float(m.body_subtreemass[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root')])*9.81;report=dict(machine='easteregg2',period_s=reg['selected']['period_s'],translation_per_cycle=drift.tolist(),position_closure_tangent_norm=float(np.linalg.norm(q[-1,:3]-q[0,:3]-drift[:3])),velocity_closure_max=float(abs(v[-1]-v[0]).max()),root_force_required_peak_BW=(abs(req[:,:3]).max(axis=0)/bw).tolist(),root_torque_required_peak_Nm=abs(req[:,3:6]).max(axis=0).tolist(),unactuated_toe_required_peak_Nm=abs(req[:,toes]).max(axis=0).tolist(),necessary_instantaneous_muscle_capacity_residual_peak_Nm=float(max(fitres)),dynamic_feasible=False,accepted_complete_cycle=False,scope='Constructed periodic pose/velocity closure, not measured motion. Acceleration=validated source-forward acceleration plus the discrete velocity-derivative correction; correction explicitly derived, not labelled measured. Actual MuJoCo inverse/contact/passive operator evaluated every2ms. Root and unactuated-toe required applied forces must vanish; muscle0..maximum capacity test is necessary only, not dynamic Hill-state reachability. No Root force/clamp or desired torque applied in a simulation, no optimizer or SNS fit. Failing this gate prevents forward-teacher promotion.');np.savez_compressed(folder/'CLOSED_REFERENCE.npz',qpos=q,qvel=v,qacc_derived=derived,required_applied_generalized_force=req);(folder/'CLOSED_REFERENCE_GATE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def project_free_acceleration(root):
    """Offline mixed inverse diagnostic, never a prescribed-body forward trial.

    Keep six source-derived joint accelerations, solve all unactuated accelerations
    for zero applied force, and retain the actual belt servo. This can diagnose a
    bad prescribed body/toe trajectory; it does not make that trajectory periodic
    or establish dynamically reachable Hill tensions.
    """
    from scipy.optimize import lsq_linear
    root=Path(root);folder=root/'periodic_reference';reg=json.loads((folder/'MEASURED_CYCLE_REGISTRATION.json').read_text());base=Path(reg['baseline']);seed=np.load(folder/'MEASURED_CYCLE_SEED.npz');closed=np.load(folder/'CLOSED_REFERENCE.npz');basis=np.load(folder/'MEASURED_ACCELERATION_BASIS.npz')['qacc_forward'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];tids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act];driven=np.array([m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n+'_'+s)] for s in 'LR' for n in ('hip','knee','ankle')]);free=np.array([i for i in range(m.nv) if i not in driven]);motor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide');records={};started=time.perf_counter()
    for name,q,v,a in [('measured_control',seed['q'][:-1],seed['qvel_reconstructed'][:-1],basis),('periodically_closed',closed['qpos'][:-1],closed['qvel'][:-1],closed['qacc_derived'])]:
        residuals=[];changes=[];fits=[];iterations=[];accs=[];forces=[];unactuated_J=[]
        for k in range(len(a)):
            mujoco.mj_resetData(m,d);d.qpos[:]=q[k];d.qvel[:]=v[k];d.ctrl[motor]=.2;mujoco.mj_forward(m,d);environment=d.qfrc_actuator.copy();A=-d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tids,driven)].T;unactuated_J.append(float(abs(d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tids,free)]).max()));current=a[k].copy()
            def evaluate(values):
                current[free]=values;d.qacc[:]=current;mujoco.mj_inverse(m,d);return (d.qfrc_inverse-environment)[free].copy()
            x=current[free].copy();r=evaluate(x)
            for it in range(12):
                if abs(r).max()<1e-5:break
                J=np.empty((len(free),len(free)))
                for j in range(len(free)):
                    step=1e-4*max(1.,abs(x[j]));xx=x.copy();xx[j]+=step;rr=evaluate(xx);J[:,j]=(rr-r)/step
                delta=np.linalg.lstsq(J,-r,rcond=1e-12)[0];best=(float(np.linalg.norm(r)),x,r)
                for scale in (1.,.5,.25,.125):
                    xx=x+scale*delta;rr=evaluate(xx);candidate=(float(np.linalg.norm(rr)),xx,rr)
                    if candidate[0]<best[0]:best=candidate
                if best[0]>=float(np.linalg.norm(r))*(1-1e-10):break
                _,x,r=best
            r=evaluate(x);tau=(d.qfrc_inverse-environment)[driven].copy();fit=lsq_linear(A,tau,bounds=(np.zeros(len(act)),p[act,5]),method='bvls',tol=1e-10);residuals.append(float(abs(r).max()));changes.append((current-a[k]).copy());fits.append(float(abs(A@fit.x-tau).max()));iterations.append(it);accs.append(current.copy());forces.append(tau)
        change=np.array(changes);report=dict(samples=len(a),inverse_free_force_residual_max=float(max(residuals)),zero_unactuated_wrench_samples=int(sum(x<1e-4 for x in residuals)),maximum_actuator_Jacobian_in_free_rows=float(max(unactuated_J)),free_acceleration_change_peak=abs(change[:,free]).max(axis=0).tolist(),necessary_capacity_residual_peak_Nm=float(max(fits)),necessary_capacity_pass_samples=int(sum(x<1e-4 for x in fits)),maximum_iterations=int(max(iterations)))
        records[name]=report
        if name=='periodically_closed':np.savez_compressed(folder/'FREE_ACCELERATION_PROJECTION.npz',qacc_projected=np.array(accs),qacc_change=change,required_driven_torque_Nm=np.array(forces),free_force_residual=np.array(residuals),necessary_capacity_residual_Nm=np.array(fits))
    report=dict(machine='easteregg2',free_dof_indices=free.tolist(),driven_dof_indices=driven.tolist(),comparisons=records,warmed_wall_s=time.perf_counter()-started,accepted_complete_cycle=False,scope='Offline mixed inverse/contact calculation only. Six source-derived joint accelerations held; all other accelerations solved for zero external applied wrench except actual belt servo. No body forces, joint torques, trajectory clamps, forward teacher, optimizer, or SNS fitting applied. Even a converged instantaneous inverse solve does not validate a periodic trajectory: projected free acceleration must integrate to the candidate body/toe trajectory and source-Hill states must dynamically reach the required muscle forces. Capacity0..sourceMaximum is necessary only. Large raw closed-reference residuals describe incompatible constructed contact accelerations, NOT forces observed in a simulation.')
    (folder/'FREE_ACCELERATION_PROJECTION_GATE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def contact_preserving_warps(root):
    """Full measured-cycle geometry warps with continuous Jacobian/RNE derivatives.

    These retain the measured near-cycle, not an artificially closed trajectory.
    Coordinate, velocity, and acceleration controls are explicitly distinguished.
    All force results remain offline inverse residuals, never applied commands.
    """
    from scipy.optimize import lsq_linear
    root=Path(root);folder=root/'periodic_reference';reg=json.loads((folder/'MEASURED_CYCLE_REGISTRATION.json').read_text());base=Path(reg['baseline']);seed=np.load(folder/'MEASURED_CYCLE_SEED.npz');qs=seed['q'][:-1];vs=seed['qvel_reconstructed'][:-1];acc=np.load(folder/'MEASURED_ACCELERATION_BASIS.npz')['qacc_forward'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));src=mujoco.MjData(m);d=mujoco.MjData(m);rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');gids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_foot_'+s) for s in 'LR'];jids=[np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n+'_'+s) for n in ('hip','knee','ankle')]) for s in 'LR'];qadr=[m.jnt_qposadr[j] for j in jids];dadr=[m.jnt_dofadr[j] for j in jids];driven=np.concatenate(dadr);free=np.array([i for i in range(m.nv) if i not in driven]);motor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide');bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];act=[i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle'];tids=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+bank[i]['name']) for i in act];bw=float(m.body_subtreemass[rid])*9.81;cases=[];started=time.perf_counter()
    def jac(data,g):
        jp=np.zeros((3,m.nv));jr=jp.copy();mujoco.mj_jacGeom(m,data,jp,jr,g);return np.vstack((jr,jp))
    def geomacc(data,g):
        out=np.zeros(6);mujoco.mj_objectAcceleration(m,data,mujoco.mjtObj.mjOBJ_GEOM,g,out,0);return out
    for dx,dy in ((0.,0.),(.01,0.),(.02,.002),(.04,.005)):
        rows=[];qkeep=[];vkeep=[];akeep=[]
        for k in range(len(acc)):
            mujoco.mj_resetData(m,src);src.qpos[:]=qs[k];src.qvel[:]=vs[k];mujoco.mj_forward(m,src);R=src.xmat[rid].reshape(3,3).copy();targets=[src.geom_xpos[g].copy() for g in gids];rots=[src.geom_xmat[g].reshape(3,3).copy() for g in gids];source_J=np.vstack([jac(src,g) for g in gids]);targetv=source_J@vs[k];src.qacc[:]=acc[k];mujoco.mj_rnePostConstraint(m,src);targeta=np.concatenate([geomacc(src,g) for g in gids]);q=qs[k].copy();v=vs[k].copy();a=acc[k].copy();shift=R@np.array([dx,dy,0.]);q[:3]+=shift;fiterr=[]
            for side,ids,qa,g,target,rot in zip('LR',jids,qadr,gids,targets,rots):
                if dx==dy==0:fiterr.append(0.);continue
                def residual(x):
                    d.qpos[:]=q;d.qpos[qa]=x;d.qvel[:]=0.;mujoco.mj_forward(m,d);local=R.T@(d.geom_xpos[g]-target);now=R.T@d.geom_xmat[g].reshape(3,3)[:,0];old=R.T@rot[:,0];angle=np.arctan2(now[1],now[0])-np.arctan2(old[1],old[0]);return np.r_[local[:2],.1*angle]
                lo=m.jnt_range[ids,0];hi=m.jnt_range[ids,1];fit=least_squares(residual,np.clip(q[qa],lo+1e-9,hi-1e-9),bounds=(lo,hi),gtol=1e-11,xtol=1e-11,ftol=1e-11,max_nfev=40);q[qa]=fit.x;fiterr.append(float(abs(fit.fun).max()))
            mujoco.mj_resetData(m,d);d.qpos[:]=q;mujoco.mj_forward(m,d);position_error=max(float(np.linalg.norm(d.geom_xpos[g]-target)) for g,target in zip(gids,targets));rotation_error=max(float(np.linalg.norm(Rotation.from_matrix(rot.T@d.geom_xmat[g].reshape(3,3)).as_rotvec())) for g,rot in zip(gids,rots));geometry_ok=max(fiterr)<1e-6 and position_error<1e-6 and rotation_error<1e-5
            if geometry_ok and (dx!=0 or dy!=0):
                jp=np.zeros((3,m.nv));jr=jp.copy();mujoco.mj_jacBody(m,src,jp,jr,rid);omega=jr@vs[k];alpha=jr@acc[k];v[:3]+=np.cross(omega,shift);a[:3]+=np.cross(alpha,shift)+np.cross(omega,np.cross(omega,shift));J=np.vstack([jac(d,g) for g in gids]);v[driven]=0.;v[driven]=np.linalg.lstsq(J[:,driven],targetv-J@v,rcond=1e-12)[0];d.qvel[:]=v;mujoco.mj_fwdVelocity(m,d);d.qacc[:]=a;d.qacc[driven]=0.;mujoco.mj_rnePostConstraint(m,d);bias=np.concatenate([geomacc(d,g) for g in gids]);a[driven]=np.linalg.lstsq(J[:,driven],targeta-bias,rcond=1e-12)[0]
            else:J=np.vstack([jac(d,g) for g in gids])
            d.qvel[:]=v;mujoco.mj_fwdVelocity(m,d);d.qacc[:]=a;mujoco.mj_rnePostConstraint(m,d);v_error=float(abs(J@v-targetv).max());a_error=float(abs(np.concatenate([geomacc(d,g) for g in gids])-targeta).max())
            if geometry_ok:
                d.ctrl[motor]=.2;mujoco.mj_fwdActuation(m,d);environment=d.qfrc_actuator.copy();d.qacc[:]=a;mujoco.mj_inverse(m,d);req=d.qfrc_inverse-environment;A=-d.ten_J.reshape(m.ntendon,m.nv)[np.ix_(tids,driven)].T;fit=lsq_linear(A,req[driven],bounds=(np.zeros(len(act)),p[act,5]),method='bvls',tol=1e-10);capacity=float(abs(A@fit.x-req[driven]).max());rootf=(req[:3]/bw).tolist();roott=req[3:6].tolist();freeforce=float(abs(req[free]).max())
            else:capacity=None;rootf=None;roott=None;freeforce=None
            knees=[float(np.degrees(q[qadr[j][1]])+(10 if j==0 else 55)) for j in range(2)];rows.append(dict(phase=k/len(acc),geometric_feasible=bool(geometry_ok),foot_position_error_m=position_error,foot_rotation_error_rad=rotation_error,foot_velocity_error=v_error,foot_acceleration_error=a_error,root_residual_force_BW=rootf,root_residual_torque_Nm=roott,free_force_residual=freeforce,necessary_capacity_residual_Nm=capacity,knee_flexion_deg=knees));qkeep.append(q);vkeep.append(v);akeep.append(a)
        valid=[r for r in rows if r['geometric_feasible']];report=dict(root_local_forward_offset_m=dx,root_local_superior_offset_m=dy,geometric_feasible_samples=len(valid),total_samples=len(rows),maximum_foot_position_error_m=max(r['foot_position_error_m'] for r in rows),maximum_foot_acceleration_error=max(r['foot_acceleration_error'] for r in valid) if valid else None,peak_inverse_root_force_BW=np.max(abs(np.array([r['root_residual_force_BW'] for r in valid])),axis=0).tolist() if valid else None,peak_inverse_root_torque_Nm=np.max(abs(np.array([r['root_residual_torque_Nm'] for r in valid])),axis=0).tolist() if valid else None,necessary_capacity_pass_samples=sum(r['necessary_capacity_residual_Nm']<1e-4 for r in valid),zero_unactuated_wrench_samples=sum(r['free_force_residual']<1e-4 for r in valid),minimum_knee_flexion_deg=np.min([r['knee_flexion_deg'] for r in rows],axis=0).tolist(),accepted_complete_cycle=False);cases.append(report)
        if dx==.02:np.savez_compressed(folder/'CONTACT_PRESERVING_WARP.npz',qpos=np.array(qkeep),qvel_transformed=np.array(vkeep),qacc_transformed=np.array(akeep))
    report=dict(machine='easteregg2',source=str(base),physics_samples=len(acc),period_s=reg['selected']['period_s'],cases=cases,warmed_wall_s=time.perf_counter()-started,accepted_complete_cycle=False,scope='Offline registered-measured-cycle body/landing geometry transformations. Foot poses retained where original nominal hip/knee/ankle ranges allow; original passive-toe coordinates retained. Velocities and accelerations transformed using actual Jacobians and MuJoCo RNE body acceleration, not finite-difference interpolated acceleration labelled measured. Actual source forward acceleration is the unmodified basis. Identity warp is a mandatory operator control. Root/carriage/rails remain actual geometry; no reaction-force cancellation. Any nonzero unactuated inverse wrench, original-ROM failure, or source-force-capacity failure rejects the full reference. Near-periodicity remains measured, not repaired by fiat. No forward teacher, applied Root wrench, trajectory clamp, optimizer, or neural fit.')
    (folder/'CONTACT_PRESERVING_WARP_GATE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def periodic_motor_shooting(root,cycles=6,phase_parameters=None,save_trace=True,artifact="PERIODIC_MOTOR_SHOOTING"):
    """Actual source-voltage/Hill forward seed for OFFLINE dynamic shooting.

    Reconstruct the original physical prefix before repeating the measured full
    MN-voltage cycle. No prescribed positions, body wrench, force allocation,
    changed muscle law, neural fitting, or claimed new SNS controller.
    """
    from native_kernel import NativeController,polynomial,hill_step
    from contact_guides import Measurements
    root=Path(root);folder=root/'periodic_reference';reg=json.loads((folder/'MEASURED_CYCLE_REGISTRATION.json').read_text());base=Path(reg['baseline']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];z=np.load(base/'trajectory.npz');sourceq=z['q'];sourceforces=z['forces'];sourcev=z['neuro'];names=z['neuron_names'].tolist();core=NativeController(rec['source_configuration']['source'],bank);nb=len(bank);vm=np.full((len(sourceq),nb),-.15)
    for pre,target,coeff in core.outputs:vm[:,target]=polynomial(coeff,sourcev[:,names.index(core.names[pre])]*.001)
    m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);dt=float(m.opt.timestep);i,j=reg['selected']['indices'];N=j-i;candidate_steps=N if phase_parameters is None else int(round(N*phase_parameters['period_scale']));steps=j+1+(cycles-1)*candidate_steps;assert cycles<=7,'Finite source belt permits at most seven registered periods, with explicit stroke guards';tids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+b['name']) for b in bank]);act=[(b,mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,row['name'])) for b,row in enumerate(bank) if row['type']=='LinearHillMuscle'];motor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide');rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');measure=Measurements(m,rid);geoms=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+n) for n in ('foot_L','toe_L','foot_R','toe_R')];belt=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');knees=[m.jnt_qposadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'knee_'+s)] for s in 'LR'];bw=float(m.body_subtreemass[rid])*9.81;d.qpos[:]=rec['source_configuration']['initial_qpos'];mujoco.mj_forward(m,d);last=d.ten_length[tids].copy();internal=np.zeros(nb);F=internal.copy();ia=internal.copy();ii=internal.copy();ib=internal.copy();hill_step(last,np.zeros(nb),vm[0],p,internal.copy(),F.copy(),ia.copy(),ii.copy(),ib.copy(),dt);parityq=parityF=0.;records=[];cycle_states=[];failure=None;started=time.perf_counter();shared={};slide_joints=[jj for jj in range(m.njnt) if m.jnt_type[jj]==mujoco.mjtJoint.mjJNT_SLIDE and m.jnt_limited[jj]]
    predicted_belt_end=float(sourceq[j,m.jnt_qposadr[m.actuator_trnid[motor,0]]]+(cycles-1)*candidate_steps*dt*.201)
    assert predicted_belt_end<1.95,'Reject predicted finite-belt boundary before launch'
    if phase_parameters is not None:
        assert set(phase_parameters)=={'period_scale','hip_lead_phase','knee_delay_phase'}
        assert .95<=phase_parameters['period_scale']<=1.05 and abs(phase_parameters['hip_lead_phase'])<=.12 and 0<=phase_parameters['knee_delay_phase']<=.15
        bindings={core.names[pre]:b for pre,b,coeff in core.outputs if bank[b]['type']=='LinearHillMuscle'};grid=np.arange(1,N+1)/N
        for name,b in bindings.items():
            if not name.startswith('L '):continue
            other=bindings['R '+name[2:]];left=1/(1+np.exp(np.clip(p[b,8]*(p[b,6]-vm[i+1:j+1,b]),-700,700)));right=1/(1+np.exp(np.clip(p[other,8]*(p[other,6]-vm[i+1:j+1,other]),-700,700)));pooled=.5*(left+np.interp(grid+.5,grid,right,period=1.));shared[name[2:]]=(b,other,pooled)
        assert len(shared)==6,bindings
    for k in range(steps):
        if any(d.qpos[m.jnt_qposadr[jj]]<=m.jnt_range[jj,0]+.05 or d.qpos[m.jnt_qposadr[jj]]>=m.jnt_range[jj,1]-.05 for jj in slide_joints):failure=dict(step=k,relative_t_s=(k-i)*dt,reason='Stopped before finite apparatus stroke:50mm safety margin');break
        index=k if k<=j else i+1+(k-i-1)%N
        actual_vm=vm[index].copy()
        if phase_parameters is not None and k>j:
            phi=(k-j)/candidate_steps
            for group,(left,right,pooled) in shared.items():
                shift=phase_parameters['hip_lead_phase'] if group.startswith('Hip') else -phase_parameters['knee_delay_phase'] if group.startswith('Knee') else 0.
                for b,offset in ((left,0.),(right,.5)):
                    recruitment=float(np.clip(np.interp(phi+shift+offset,grid,pooled,period=1.),1e-12,1-1e-12));actual_vm[b]=p[b,6]-np.log(1/recruitment-1)/p[b,8]
            assert min(actual_vm[b] for b,a in act)>=-.100001 and max(actual_vm[b] for b,a in act)<=.000001
        mujoco.mj_step1(m,d);length=d.ten_length[tids].copy();velocity=(length-last)/dt;last[:]=length;hill_step(length,velocity,actual_vm,p,internal,F,ia,ii,ib,dt)
        for b,a in act:d.ctrl[a]=F[b]
        d.ctrl[motor]=.2 if k*dt>=.5 else 0.;mujoco.mj_step2(m,d);mujoco.mj_forward(m,d)
        if k<=j:
            parityq=max(parityq,float(abs(d.qpos-sourceq[k]).max()));parityF=max(parityF,float(abs(F-sourceforces[k]).max()))
        if k>=i:
            if k==i or k==j or (k>j and (k-j)%candidate_steps==0):cycle_states.append(dict(cycle=0 if k==i else 1+(k-j)//candidate_steps,qpos=d.qpos.copy(),qvel=d.qvel.copy(),internal=internal.copy(),previous_length=last.copy()))
            measured=measure.sample(d);top=float(d.geom_xpos[belt,1]+m.geom_size[belt,1]);clearance=np.array([d.geom_xpos[g,1]-abs(d.geom_xmat[g].reshape(3,3)[1])@m.geom_size[g]-top for g in geoms]);records.append(dict(t=(k-i)*dt,qpos=d.qpos.copy(),qvel=d.qvel.copy(),force=F.copy(),voltage=actual_vm.copy(),foot_force=measured['foot_force_world'].copy(),clearance=clearance,com=measured['com_position'].copy(),cop=measured['foot_cop_world'].copy(),rail=measured['rail_force_world'].copy(),guide=measured['guide_force_world'].copy()))
        knee=np.degrees(d.qpos[knees])+[10.,55.]
        if not np.isfinite(d.qpos).all() or min(knee)<0 or d.qpos[1]<.65:failure=dict(step=k,relative_t_s=(k-i)*dt,reason='Nonfinite, protected knee hyperextension, or body fall',knee_flexion_deg=knee.tolist());break
        if k==j and (parityq>1e-7 or parityF>1e-6):failure=dict(step=k,reason='Original full-physics MN-voltage replay failed parity before periodic extension',max_qpos_error=parityq,max_force_error_N=parityF);break
    closure=[]
    for old,new in zip(cycle_states[:-1],cycle_states[1:]):
        tangent=np.zeros(m.nv);mujoco.mj_differentiatePos(m,tangent,1.,old['qpos'],new['qpos']);closure.append(dict(cycle=new['cycle'],root_yz_closure_m=float(np.linalg.norm(tangent[1:3])),root_orientation_closure_deg=float(np.degrees(np.linalg.norm(tangent[3:6]))),root_velocity_closure_m_s=(new['qvel'][:3]-old['qvel'][:3]).tolist(),maximum_joint_angle_closure_deg=float(abs(np.degrees(tangent[6:14])).max()),maximum_Hill_state_closure_N=float(abs(new['internal']-old['internal']).max()),world_root_translation_m=float(new['qpos'][0]-old['qpos'][0]),belt_relative_progress_m=float(new['qpos'][0]-old['qpos'][0]+.2*(N if new['cycle']==1 else candidate_steps)*dt)))
    if records:
        arrays={key:np.array([row[key] for row in records]) for key in records[0]};
        if save_trace:np.savez_compressed(folder/(artifact+'.npz'),**arrays,cycle_internal=np.array([row['internal'] for row in cycle_states]),cycle_previous_length=np.array([row['previous_length'] for row in cycle_states]),cycle_qpos=np.array([row['qpos'] for row in cycle_states]),cycle_qvel=np.array([row['qvel'] for row in cycle_states]))
        assistance=arrays['rail'][:,1]+arrays['guide'][:,1];normal=arrays['foot_force'][:,:,1];clear=np.minimum(arrays['clearance'][:,::2],arrays['clearance'][:,1::2]);unload_clear=((normal<.025*bw)&(clear>.01));min_knee=np.min(np.degrees(arrays['qpos'][:,knees])+[10.,55.],axis=0).tolist();summary=dict(minimum_knee_flexion_deg=min_knee,mean_upward_assistance_BW=float((np.maximum(arrays['rail'][:,1],0)+np.maximum(arrays['guide'][:,1],0)).mean()/bw),mean_downward_assistance_BW=float((np.maximum(-arrays['rail'][:,1],0)+np.maximum(-arrays['guide'][:,1],0)).mean()/bw),unloaded_and_cleared_fraction=unload_clear.mean(axis=0).tolist(),source_prefix_max_qpos_error=parityq,source_prefix_max_force_error_N=parityF)
        mask=arrays['t']>=N*dt;events=[];counts=[];peak_clear=[]
        for side in range(2):
            swing=mask&(clear[:,side]>.001)&(normal[:,side]<.025*bw);edges=np.diff(np.r_[False,swing,False].astype(int));starts=np.flatnonzero(edges==1);ends=np.flatnonzero(edges==-1);qualified=[]
            for first,lastindex in zip(starts,ends):
                if (lastindex-first)*dt>=.08 and clear[first:lastindex,side].max()>=.01 and np.mean(normal[first:lastindex,1-side]>.05*bw)>=.5:qualified.append(float(arrays['t'][first]));events.append((float(arrays['t'][first]),side))
            counts.append(len(qualified));peak_clear.append(float(clear[mask,side].max()) if mask.any() else 0.)
        events.sort();alternation=sum(a[1]!=b[1] for a,b in zip(events[:-1],events[1:]))/max(1,len(events)-1);summary.update(qualified_swings=counts,alternating_fraction=alternation,supported_stepping_gate=bool(min(counts)>=3 and alternation>=.8 and failure is None),maximum_whole_foot_clearance_m=peak_clear,peak_foot_normal_BW=float(np.maximum(normal[mask],0).max()/bw) if mask.any() else None,p95_foot_normal_BW=float(np.percentile(np.maximum(normal[mask],0),95)/bw) if mask.any() else None,belt_relative_progress_m=float(arrays['qpos'][-1,0]-arrays['qpos'][0,0]+.2*(arrays['t'][-1]-arrays['t'][0])))
        stops={}
        for joint in range(m.njnt):
            if m.jnt_type[joint]!=mujoco.mjtJoint.mjJNT_HINGE or not m.jnt_limited[joint]:continue
            angle=arrays['qpos'][:,m.jnt_qposadr[joint]];stops[mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,joint)]=float(np.degrees(max(0.,m.jnt_range[joint,0]-angle.min(),angle.max()-m.jnt_range[joint,1])))
        summary['nominal_stop_overshoot_deg']=stops
        candidate=arrays['t']>N*dt;summary['candidate_nominal_stop_overshoot_deg']={}
        for joint in range(m.njnt):
            if m.jnt_type[joint]!=mujoco.mjtJoint.mjJNT_HINGE or not m.jnt_limited[joint] or not candidate.any():continue
            angle=arrays['qpos'][candidate,m.jnt_qposadr[joint]];summary['candidate_nominal_stop_overshoot_deg'][mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,joint)]=float(np.degrees(max(0.,m.jnt_range[joint,0]-angle.min(),angle.max()-m.jnt_range[joint,1])))
        positive=np.maximum(normal,0.);weighted=np.nansum(arrays['cop'][:,:,0]*positive,axis=1)/np.maximum(positive.sum(axis=1),1e-9);loaded=candidate&(positive.sum(axis=1)>.05*bw);summary['median_loaded_CoP_ahead_COM_m']=float(np.median((weighted-arrays['com'][:,0])[loaded])) if loaded.any() else None
    else:summary={}
    report=dict(machine='easteregg2',source=str(base),source_sha256=rec['source_sha256'],period_s=candidate_steps*dt,source_basis_period_s=N*dt,phase_parameters=phase_parameters,requested_cycles=cycles,completed_physical_steps=k+1,warmed_wall_s=time.perf_counter()-started,failure=failure,summary=summary,cycle_closure=closure,accepted_complete_cycle=False,scope='Offline forward mechanical-control replay, NOT a new SNS controller or Hunt fit. Actual recorded source MN voltages excite unchanged pooled source sigmoid/force-length/velocity/passive/series-elastic Hill equations, actual contacts/toes, and original passive guides/rails. Whole source physical prefix reconstructed, then exact complete measured voltage period repeated; source dynamic replay parity required before extension. For explicit shared-wave candidates, left/right sigmoid recruitment fractions are pooled after half-cycle registration, then inverted through each unchanged source sigmoid; only shared hip/knee phase and period change. This is offline voltage teaching, not a neural architecture or learned conductance. No prescribed body/foot trajectory, Root wrench, force allocation, changed muscle strength, neural-feedback correction, or trajectory clamp. Finite original belt preserved; at most seven registered periods with50mm stroke safety guards. Unloaded+10mm-clear fraction is a diagnostic, not an alternating-event gate or independent walking claim. Complete-cycle feasibility/periodic closure and all source ROM/current gates still required before promotion.')
    
    if save_trace:(folder/(artifact+'_GATE.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report




def audit_integrated_cycle(root):
    """Signed full-cycle control check; periodic closure is a separate gate."""
    from native_kernel import hill_step
    root=Path(root);folder=root/'periodic_reference';r=json.loads((folder/'SHARED_MOTOR_SHOOTING_GATE.json').read_text());z=np.load(folder/'SHARED_MOTOR_SHOOTING.npz');base=Path(r['source']);bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);inv=mujoco.MjData(m);tids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+b['name']) for b in bank]);act=[(b,mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,row['name'])) for b,row in enumerate(bank) if row['type']=='LinearHillMuscle'];motor=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_ACTUATOR,'motor_Treadmill_Slide');driven=[m.jnt_dofadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n+'_'+s)] for s in 'LR' for n in ('hip','knee','ankle')];free=[i for i in range(m.nv) if i not in driven];dt=float(m.opt.timestep);N=int(round(r['period_s']/dt));start=len(z['t'])-N-1;q=z['qpos'];v=z['qvel'];F=z['force'];voltage=z['voltage'];internal=z['cycle_internal'][-2].copy();last=z['cycle_previous_length'][-2].copy();force=np.zeros(len(bank));ia=force.copy();ii=force.copy();ib=force.copy();errors=[];signed=[];started=time.perf_counter()
    for k in range(start,len(q)-1):
        mujoco.mj_resetData(m,d);d.qpos[:]=q[k];d.qvel[:]=v[k];mujoco.mj_step1(m,d);length=d.ten_length[tids].copy();velocity=(length-last)/dt;last[:]=length;hill_step(length,velocity,voltage[k+1],p,internal,force,ia,ii,ib,dt)
        for b,a in act:d.ctrl[a]=force[b]
        d.ctrl[motor]=.2;mujoco.mj_step2(m,d);a=d.qacc.copy();applied=d.qfrc_actuator.copy();delta=np.zeros(m.nv);mujoco.mj_differentiatePos(m,delta,1.,q[k+1],d.qpos);mujoco.mj_resetData(m,inv);inv.qpos[:]=q[k];inv.qvel[:]=v[k];inv.qacc[:]=a;mujoco.mj_inverse(m,inv);errors.append([float(abs(force-F[k+1]).max()),float(abs(delta).max()),float(abs(d.qvel-v[k+1]).max()),float(abs(inv.qfrc_inverse-applied).max()),float(abs(applied[free][:8]).max())]);signed.append(applied[driven].copy())
    errors=np.array(errors);joint_stops={};knees=[]
    for jid in range(m.njnt):
        if m.jnt_type[jid]!=mujoco.mjtJoint.mjJNT_HINGE or not m.jnt_limited[jid]:continue
        name=mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_JOINT,jid);angle=q[start:,m.jnt_qposadr[jid]];joint_stops[name]=float(np.degrees(max(0.,m.jnt_range[jid,0]-angle.min(),angle.max()-m.jnt_range[jid,1])))
        if name.startswith('knee_'):knees.append(float(np.degrees(angle.min())+(10 if name.endswith('L') else 55)))
    actual=np.array([b for b,a in act]);force_valid=bool((F[start:,actual]>=-1e-10).all() and (F[start:,actual]<=p[actual,5]+1e-10).all());replay=bool((errors.max(axis=0)<np.array([1e-6,1e-7,1e-7,1e-4,1e-8])).all());rom=bool(min(knees)>=0 and max(joint_stops.values())<=.6);report=dict(machine='easteregg2',source=r['source'],source_sha256=r['source_sha256'],condition='Explicitly mirrored offline MN recruitment waveform, unchanged source physics, full final cycle',physical_intervals=len(errors),period_s=r['period_s'],maximum_errors=dict(zip(('dynamic_Hill_tension_N','next_pose_tangent','next_velocity','inverse_applied_generalized_force','applied_Root_and_toe_wrench'),errors.max(axis=0).tolist())),signed_actual_joint_torque_range_Nm=np.column_stack((np.min(signed,axis=0),np.max(signed,axis=0))).tolist(),nominal_joint_stop_overshoot_deg=joint_stops,minimum_actual_knee_flexion_deg=knees,dynamic_replay_pass=replay,source_force_bounds_pass=force_valid,original_ROM_with_0p6_degree_soft_stop_tolerance_pass=rom,physical_cycle_feasible=bool(replay and force_valid and rom),periodic_state_closure_established=False,periodic_state_residual=r['cycle_closure'][-1],ready_for_neural_fit=False,ready_for_optimization=False,warmed_wall_s=time.perf_counter()-started,scope='Actual complete integrated461-step final cycle with source Hill internal and previous-length state. Source Hill tensions dynamically replayed, source muscle actuation only, zero applied Root/toe wrench, same floor/rails/guides. Original nominal joint ranges retained; explicit0.6degree engineering soft-stop tolerance, with strict anatomical knee>=0. No body clamp, allocator, changed geometry/strength, or neural fit. Full-cycle mechanical feasibility is distinct from an established periodic full-state orbit; the latter remains unaccepted. Automatic phase-search rejection is not bypassed by this diagnostic.')
    (folder/'INTEGRATED_CYCLE_PHYSICAL_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def audit_fixed_wave_drift(root):
    """Saved-data, full-state/load-phase Poincare diagnosis; no fit or new run."""
    from native_kernel import NativeController,polynomial,hill_step
    root=Path(root);folder=root/'periodic_reference';stem='SHARED_MOTOR_SHOOTING_LONG7';r=json.loads((folder/(stem+'_GATE.json')).read_text());protocol=json.loads((folder/'FIXED_WAVE_CLOSURE_PROTOCOL.json').read_text());base=Path(r['source']);rec=json.loads((base/'SPIKING_COMPARISON.json').read_text());bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];p=np.load(base/'POOLED_FORCE_PARAMETERS.npz')['pooled'];m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);dt=float(m.opt.timestep);N=int(round(r['period_s']/dt));tids=np.array([mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+b['name']) for b in bank]);actual=np.array([i for i,b in enumerate(bank) if b['type']=='LinearHillMuscle']);alpha=1-dt*p[:,2]/p[:,4]*(1+p[:,3]/p[:,2]);beta=dt*p[:,2]/p[:,4];zr=np.load(folder/(stem+'.npz'));q=zr['qpos'];v=zr['qvel'];vm=zr['voltage'];F=zr['force'];bw=float(m.body_subtreemass[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root')])*9.81
    def reconstruct(q,vm,force,seed_internal=None,seed_length=None,home=None):
        S=np.zeros((len(q),len(bank)));L=S.copy();V=S.copy();active=S.copy();internal=np.zeros(len(bank)) if seed_internal is None else seed_internal.copy();last=None if seed_length is None else seed_length.copy();f=internal.copy();ia=f.copy();ii=f.copy();ib=f.copy();err=0.
        if seed_internal is not None:S[0]=internal;L[0]=last
        for k in range(0 if seed_internal is None else 1,len(q)):
            d.qpos[:]=home if k==0 else q[k-1];mujoco.mj_forward(m,d);length=d.ten_length[tids].copy();velocity=np.zeros(len(bank)) if last is None else (length-last)/dt;last=length.copy();hill_step(length,velocity,vm[k],p,internal,f,ia,ii,ib,dt);S[k]=internal;L[k]=length;V[k]=velocity;fl=np.maximum(0.,1-((length-p[:,0])/p[:,1])**2);active[k]=np.maximum(0.,p[:,7]/(1+np.exp(np.clip(p[:,8]*(p[:,6]-vm[k]),-700,700)))+p[:,9])*fl;err=max(err,float(abs(f-force[k]).max()))
        return S,L,V,active,err
    S,L,V,active,error=reconstruct(q,vm,F,zr['cycle_internal'][0],zr['cycle_previous_length'][0]);assert error<1e-7
    def omega(qq,vv):
        R=np.empty(9);mujoco.mju_quat2Mat(R,qq[3:7]);return R.reshape(3,3)@vv[3:6]
    carriage=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'Passive_guide_X_slider');ca=m.jnt_qposadr[carriage];cv=m.jnt_dofadr[carriage];direction=float(m.jnt_axis[carriage,0]);belt=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,'Treadmill_Slide');ba=m.jnt_qposadr[belt];bv=m.jnt_dofadr[belt];belt_direction=float(m.jnt_axis[belt,0])
    def compare(q,v,S,L,a,b):
        delta=np.zeros(m.nv);mujoco.mj_differentiatePos(m,delta,1.,q[a],q[b]);dv=v[b]-v[a]
        return dict(root_yz_m=float(np.linalg.norm(delta[1:3])),root_orientation_deg=float(np.degrees(np.linalg.norm(delta[3:6]))),root_linear_velocity_m_s=float(abs(dv[:3]).max()),root_world_angular_velocity_deg_s=float(abs(np.degrees(omega(q[b],v[b])-omega(q[a],v[a]))).max()),limb_angle_deg=float(abs(np.degrees(delta[6:14])).max()),limb_velocity_deg_s=float(abs(np.degrees(dv[6:14])).max()),actual_muscle_internal_state_N=float(abs(S[b,actual]-S[a,actual]).max()),actual_previous_length_m=float(abs(L[b,actual]-L[a,actual]).max()),belt_velocity_m_s=float(abs(dv[bv])),carriage_relative_root_position_m=float(abs(direction*(q[b,ca]-q[a,ca])-delta[0])),carriage_relative_root_velocity_m_s=float(abs(direction*dv[cv]-dv[0])),actual_belt_relative_progress_m=float(q[b,0]-q[a,0]-belt_direction*(q[b,ba]-q[a,ba])))
    boundaries=[dict(cycle=k,**compare(q,v,S,L,(k-1)*N,k*N)) for k in range(1,len(zr['cycle_internal']))];thresholds={k:value for k,value in protocol['closure_gate'].items() if k!='consecutive_last_cycles'};primary=boundaries[-protocol['closure_gate']['consecutive_last_cycles']:];gate={key:bool(max(row[key] for row in primary)<=value) for key,value in thresholds.items()};phasebins=[]
    for phase in np.arange(8)/8:
        offset=int(round(phase*N));last_start=(len(q)-1-N-offset)//N;rows=[compare(q,v,S,L,k*N+offset,(k+1)*N+offset) for k in range(last_start-2,last_start+1)];indices=[k*N+offset for k in range(last_start-2,last_start+1)];phasebins.append(dict(phase=float(phase),indices=indices,maximum_residuals={key:max(row[key] for row in rows) for key in thresholds},mean_foot_normal_BW=np.mean(zr['foot_force'][indices,:,1]/bw,axis=0).tolist(),minimum_whole_foot_clearance_m=np.mean(np.minimum(zr['clearance'][indices,::2],zr['clearance'][indices,1::2]),axis=0).tolist()))
    start=(len(zr['cycle_internal'])-3)*N;aa=np.arange(start,start+N);bb=aa+N;ds=S[bb]-S[aa];dl=L[bb]-L[aa];dv=V[bb]-V[aa];da=active[bb]-active[aa];pred=alpha*ds[:-1]+beta*(p[:,3]*dl[1:]+p[:,4]*dv[1:]+da[1:]);recurrence=float(abs(pred-ds[1:]).max());parts=[alpha*ds[:-1],beta*p[:,3]*dl[1:],beta*p[:,4]*dv[1:],beta*da[1:]];dominance={name:float(np.mean(abs(part[:,actual]))) for name,part in zip(('previous_step_state','parallel_elastic_length','length_velocity','active_force_length_and_voltage'),parts)};peak=np.unravel_index(abs(ds[:,actual]).argmax(),ds[:,actual].shape);peak_i=int(aa[peak[0]]);peak_muscle=int(actual[peak[1]])
    original=np.load(base/'trajectory.npz');bq=original['q'];bforce=original['forces'];bvdata=np.zeros((len(bq),m.nv));home=np.array(rec['source_configuration']['initial_qpos'])
    for k in range(len(bq)):mujoco.mj_differentiatePos(m,bvdata[k],dt,home if k==0 else bq[k-1],bq[k])
    c=NativeController(rec['source_configuration']['source'],bank);bvm=np.full((len(bq),len(bank)),-.15);names=original['neuron_names'].tolist();neuro=original['neuro']
    for pre,target,coeff in c.outputs:bvm[:,target]=polynomial(coeff,neuro[:,names.index(c.names[pre])]*.001)
    BS,BL,_,_,berror=reconstruct(bq,bvm,bforce,home=home);assert berror<1e-7;case=next(row for row in json.loads((root/'INTERLOCK_GATES.json').read_text())['cases'] if Path(row['path'])==base);strike=[int(round(row['touchdown_s']/dt)) for row in case['landings'] if row['side']=='L' and row['touchdown_s']>=2];baseline=[dict(start_s=a*dt,end_s=b*dt,**compare(bq,bvdata,BS,BL,a,b)) for a,b in zip(strike[:-1],strike[1:])]
    terrain={};margin=float('inf');railmargin=float('inf');footgeom=[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+x) for x in ('foot_L','toe_L','foot_R','toe_R')];beltgeom=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');rails=[g for g in range(m.ngeom) if 'support_rail' in (mujoco.mj_id2name(m,mujoco.mjtObj.mjOBJ_GEOM,g) or '')]
    for qq in q:
        d.qpos[:]=qq;mujoco.mj_forward(m,d);beltx=d.geom_xpos[beltgeom,0];half=abs(d.geom_xmat[beltgeom].reshape(3,3)[0])@m.geom_size[beltgeom]
        for g in footgeom:
            extent=abs(d.geom_xmat[g].reshape(3,3)[0])@m.geom_size[g];margin=min(margin,d.geom_xpos[g,0]-extent-(beltx-half),(beltx+half)-(d.geom_xpos[g,0]+extent))
        for g in rails:
            extent=abs(d.geom_xmat[g].reshape(3,3)[0])@m.geom_size[g];railmargin=min(railmargin,d.qpos[0]-(d.geom_xpos[g,0]-extent),(d.geom_xpos[g,0]+extent)-d.qpos[0])
    for name,jj in [('belt',belt),('carriage',carriage)]:
        pos=q[:,m.jnt_qposadr[jj]];terrain[name]=dict(range_m=m.jnt_range[jj].tolist(),observed_m=[float(pos.min()),float(pos.max())],minimum_stroke_margin_m=float(min(pos.min()-m.jnt_range[jj,0],m.jnt_range[jj,1]-pos.max())))
    terrain.update(minimum_all_foot_box_longitudinal_belt_margin_m=float(margin),minimum_Root_center_longitudinal_rail_margin_m=float(railmargin));memory=[dict(route=bank[k]['name'],one_step_internal_decay=float(alpha[k]),time_constant_s=float(-dt/np.log(abs(alpha[k]))),one_period_initial_state_retention=float(abs(alpha[k])**N)) for k in actual];report=dict(machine='easteregg2',condition=stem,source_hash=r['source_sha256'],protocol=str(folder/'FIXED_WAVE_CLOSURE_PROTOCOL.json'),source_Hill_replay_error_N=error,live_baseline_Hill_replay_error_N=berror,period_boundaries=boundaries,preregistered_last3_criterion_pass=gate,preregistered_near_periodic_gate=bool(all(gate.values())),eight_prespecified_phase_bins=phasebins,live_baseline_actual_L_strike_residuals=baseline,Hill_memory=memory,last_window_exact_difference_recurrence_error_N=recurrence,last_window_mean_absolute_recurrence_terms_N=dominance,peak_Hill_state_difference=dict(route=bank[peak_muscle]['name'],phase=float((peak_i-start)/N),difference_N=float(ds[peak[0],peak_muscle]),paired_internal_states_N=[float(S[peak_i,peak_muscle]),float(S[peak_i+N,peak_muscle])],paired_actual_tension_N=[float(F[peak_i,peak_muscle]),float(F[peak_i+N,peak_muscle])],paired_length_velocity_m_s=[float(V[peak_i,peak_muscle]),float(V[peak_i+N,peak_muscle])],foot_normal_BW=(zr['foot_force'][peak_i,:,1]/bw).tolist(),paired_next_foot_normal_BW=(zr['foot_force'][peak_i+N,:,1]/bw).tolist(),whole_foot_clearance_m=np.minimum(zr['clearance'][peak_i,::2],zr['clearance'][peak_i,1::2]).tolist()),terrain=terrain,optimizer_started=False,ready_for_neural_fit=False,scope='Saved-data analysis only following ONE authorized unchanged longer run. Preregistered last3period-boundary closure gate uses every state component; eight uniformly declared phases are diagnostic and cannot replace a failed primary gate. Live SNS baseline uses actual force-qualified left strikes, whereas replay boundaries follow its fixed teaching clock; distributions are not matched interventions. Hill memory separates decayed initial state from ongoing kinematic/voltage forcing. No phase search or workaround of its rejection, no newly fitted neural controller, no support promotion.')
    report['actual_cumulative_belt_relative_progress_m']=float(q[-1,0]-q[0,0]-belt_direction*(q[-1,ba]-q[0,ba]))
    report['commanded_belt_reference_progress_m']=r['summary']['belt_relative_progress_m']
    report['progress_scope']='Actual displacement uses measured slider position/axis. The earlier mechanical-replay summary key belt_relative_progress_m used commanded0.2m/s integration; retain that value as an explicitly labelled approximation, not measured travel. Preferred live-SNS baseline progress0.516969809m was already computed from actual belt travel and is unchanged.'
    report['last_window_max_actual_tension_difference_N']=float(abs(F[bb][:,actual]-F[aa][:,actual]).max())
    report['last_window_max_command_voltage_difference_mV']=float(abs(vm[bb][:,actual]-vm[aa][:,actual]).max()*1000)
    report['live_baseline_strike_residual_medians']={key:float(np.median([row[key] for row in baseline])) for key in thresholds}
    report['fixed_last3_boundary_residual_maxima']={key:float(max(row[key] for row in primary)) for key in thresholds}
    report['review_action']='Do not promote a periodic reference or retry the rejected phase search. Existing live SNS supported stepping remains valid and preserved. Residual forcing is kinematic/contact-cycle variation, not slowly forgotten initial Hill state; further warmup or generic motor-strength sweeps are not justified by this audit.'
    (folder/'FIXED_WAVE_POINCARE_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def handle_request(request):
    """Bounded Rustuna adapter for OFFLINE periodic motor-wave shooting only."""
    import hashlib
    cfg=request['manifest'];root=Path(cfg['reference_root']);folder=root/'periodic_reference'
    source=Path(json.loads((folder/'MEASURED_CYCLE_REGISTRATION.json').read_text())['baseline'])
    rec=json.loads((source/'SPIKING_COMPARISON.json').read_text())
    assert rec['source_sha256']==cfg['expected_source_sha256']
    assert hashlib.sha256(Path(rec['source_configuration']['source']).read_bytes()).hexdigest()==cfg['expected_source_sha256']
    if request['op']=='preflight':
        control=json.loads((folder/'MEASURED_DYNAMICS_CONTROL_GATE.json').read_text());seed=json.loads((folder/'SHARED_MOTOR_SHOOTING_GATE.json').read_text());ready=control['inverse_control_pass'] and seed['summary']['supported_stepping_gate'] and seed['summary']['source_prefix_max_qpos_error']==0 and seed['summary']['source_prefix_max_force_error_N']==0
        return dict(ok=True,ready_for_optimization=False,reason='Automatic review rejected phase search while complete-cycle acceptance remains unresolved. Source replay/stepping pass alone does not authorize optimization.',physical_seed_ready=bool(ready),scope='Offline phase-only mechanical teacher search, no new SNS success',baseline_upward_BW=seed['summary']['mean_upward_assistance_BW'])
    assert request['op']=='evaluate'
    r=periodic_motor_shooting(root,phase_parameters=request['params'],save_trace=False);s=r['summary'];last=r['cycle_closure'][-1] if r['cycle_closure'] else None
    physical=s.get('supported_stepping_gate',False) and s.get('belt_relative_progress_m',-1)>0 and min(s.get('minimum_knee_flexion_deg',[-1]))>=0
    stops=s.get('candidate_nominal_stop_overshoot_deg',{});stop_ok=bool(stops) and max(stops.values())<=.6
    if last:
        closure=last['root_yz_closure_m']/.001+last['root_orientation_closure_deg']/.25+max(abs(x) for x in last['root_velocity_closure_m_s'])/.005+last['maximum_joint_angle_closure_deg']/.25+last['maximum_Hill_state_closure_N']/10
    else:closure=100.
    if physical and stop_ok:
        value=s['mean_upward_assistance_BW']+s['mean_downward_assistance_BW']+.005*min(closure,100)+.005*max(0,s['peak_foot_normal_BW']-5.65)
    else:value=10.+(0 if physical else 1)+.01*min(closure,100)
    return dict(ok=True,value=float(value),diagnostics=dict(scope='Offline physical teacher, not neural fit',physical_step_gate=bool(physical),nominal_stop_engineering_gate=bool(stop_ok),summary=s,last_cycle_closure=last,closure_penalty=closure,failure=r['failure'],warmed_wall_s=r['warmed_wall_s']))


def geometry(root,period=.912,duty=.65,height=.011,pitch_deg=-2.8895680781,heel_up_deg=0.,samples=48,speed=.1):
    root=Path(root);anchor=json.loads((root/'LANDING_PHASE_PAIR_TARGETS.json').read_text());base=Path(anchor['baseline']);m=mujoco.MjModel.from_xml_path(str(base/'model.xml'));d=mujoco.MjData(m);pair=np.load(root/'LANDING_PHASE_PAIR_TARGETS.npz')['qpos'];qbase=pair.mean(axis=0);qbase[3:7]=np.roll(Rotation.from_euler('z',-pitch_deg,degrees=True).as_quat(),1);qbase[2]=0.;rid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');belt=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_Treadmill_Belt');signs=np.array(list(itertools.product((-1,1),repeat=3)));dt=period/samples
    ids={s:[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,n+'_'+s) for n in ('hip','knee','ankle','toe')] for s in 'LR'};adrs={s:m.jnt_qposadr[ids[s]] for s in 'LR'};foot={s:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_foot_'+s) for s in 'LR'};toe={s:mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_toe_'+s) for s in 'LR'}
    def vertices(g):return d.geom_xpos[g]+(signs*m.geom_size[g])@d.geom_xmat[g].reshape(3,3).T
    d.qpos[:]=pair[0];mujoco.mj_forward(m,d);comx=float(d.subtree_com[rid,0]);heel=vertices(foot['R']);heel=heel[signs[:,0]<0];heelx=float(heel[np.argmin(heel[:,1]),0]);tv=vertices(toe['R']);toe_rear_offset=float(tv[np.argmin(tv[:,1]),0]-heelx);heel_relative=-.015
    roll_amplitude=abs(Rotation.from_quat(np.roll(pair[0,3:7],-1)).as_euler('xyz')[0]);qbase[3:7]=np.roll(Rotation.from_euler('xyz',[-roll_amplitude,0.,np.radians(-pitch_deg)]).as_quat(),1)
    d.qpos[:]=qbase;mujoco.mj_forward(m,d);gaps=[]
    for side in ('Left','Right'):
        g=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+side+'_hip_wing');rail=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_GEOM,'geom_'+side+'_support_rail');gaps.append(float(vertices(g)[:,1].min()-d.geom_xpos[rail,1]-m.geom_size[rail,1]))
    qbase[1]-=min(gaps);d.qpos[:]=qbase;mujoco.mj_forward(m,d);top=float(d.geom_xpos[belt,1]+m.geom_size[belt,1]);land_heel=comx+heel_relative+speed*period*duty/2
    bounds_lo=[qbase[0]-.08];bounds_hi=[qbase[0]+.08]
    for s in 'LR':
        bounds_lo.extend([*m.jnt_range[ids[s][:3],0],-.7]);bounds_hi.extend([*m.jnt_range[ids[s][:3],1],.9])
    x=np.r_[qbase[0],qbase[adrs['L']],qbase[adrs['R']]];x=np.clip(x,np.array(bounds_lo)+1e-8,np.array(bounds_hi)-1e-8);qs=[];rows=[];started=time.perf_counter()
    def desired(phase):
        if phase<duty:
            theta=np.radians(.5-(heel_up_deg+.5)*smooth((phase-.5)/(duty-.5))) if heel_up_deg else np.radians(.5)
            pivot=bool(heel_up_deg and phase>=.5);cx=land_heel-speed*period*phase+(toe_rear_offset if pivot else 0.);return dict(phase=phase,stance=True,theta=theta,toe_pivot=pivot,contact_x=cx)
        u=(phase-duty)/(1-duty);theta=-np.radians(heel_up_deg)+(np.radians(5+heel_up_deg))*smooth(u);a=-np.radians(heel_up_deg) if heel_up_deg else np.radians(.5);b=np.radians(.5);x0=land_heel+.12-speed*period*duty-.06*np.cos(a)-.01*np.sin(a);x1=land_heel+.06*np.cos(b)-.01*np.sin(b);v=-speed*period*(1-duty);cx=x0+v*u+(x1-x0-v)*smooth(u);return dict(phase=phase,stance=False,theta=theta,center_x=cx,height=height*smooth(u/.3)*smooth((1-u)/.3))
    for k in range(samples):
        phase=k/samples;targets={s:desired((phase+(.5 if s=='L' else 0))%1) for s in 'LR'}
        ds=duty-.5
        roll=roll_amplitude*(1-2*smooth(phase/ds)) if phase<ds else -roll_amplitude if phase<.5 else roll_amplitude*(-1+2*smooth((phase-.5)/ds)) if phase<duty else roll_amplitude
        phase_base=qbase.copy();phase_base[3:7]=np.roll(Rotation.from_euler('xyz',[roll,0.,np.radians(-pitch_deg)]).as_quat(),1)
        def residual(values,details=False):
            q=phase_base.copy();q[0]=values[0]
            for j,s in enumerate('LR'):q[adrs[s]]=values[1+4*j:5+4*j]
            d.qpos[:]=q;d.qvel[:]=0;mujoco.mj_forward(m,d);r=[float(d.subtree_com[rid,0]-comx)];stats={}
            for s in 'LR':
                a=targets[s];fv=vertices(foot[s]);tv=vertices(toe[s]);theta=np.arctan2(d.geom_xmat[foot[s]].reshape(3,3)[1,0],d.geom_xmat[foot[s]].reshape(3,3)[0,0]);errtheta=np.arctan2(np.sin(theta-a['theta']),np.cos(theta-a['theta']))
                if a['stance']:
                    if a['toe_pivot']:
                        j=int(np.argmin(tv[:,1]));point=tv[j];xgoal=a['contact_x']+(.04 if signs[j,0]>0 else 0.);toe_goal=np.radians(1.5571518638)+.5*max(0.,-a['theta']);r.extend([point[0]-xgoal,point[1]-top,.05*errtheta,.05*(q[adrs[s][3]]-toe_goal)])
                    else:
                        edge=fv[signs[:,0]<0];point=edge[np.argmin(edge[:,1])];r.extend([point[0]-a['contact_x'],point[1]-top,.05*errtheta,tv[:,1].min()-top-.00005])
                else:
                    tdof=m.jnt_dofadr[ids[s][3]];r.extend([d.geom_xpos[foot[s],0]-a['center_x'],min(fv[:,1].min(),tv[:,1].min())-top-a['height'],.05*(q[adrs[s][2]]-m.jnt_range[ids[s][2],1]),.01*(d.qfrc_bias[tdof]-d.qfrc_passive[tdof])])
                stats[s]=dict(**a,foot_min_gap_m=float(fv[:,1].min()-top),toe_min_gap_m=float(tv[:,1].min()-top),actual_foot_pitch_deg=float(np.degrees(theta)))
            return (np.array(r),q,stats) if details else np.array(r)
        lo=np.array(bounds_lo);hi=np.array(bounds_hi)
        for j,s in enumerate('LR'):
            if targets[s]['stance']:lo[4+4*j]=0.
            else:lo[4+4*j]=-.15;hi[4+4*j]=.15
        seeds=[np.clip(x,lo+1e-8,hi-1e-8)];alternate=x.copy()
        for j,s in enumerate('LR'):
            a=targets[s];hip,knee,ankle=(12.,30.,3.) if not a['stance'] else (5.,20.,-5.) if a['toe_pivot'] else (0.,4.,1.);alternate[1+4*j:5+4*j]=np.radians([hip+(10 if s=='L' else -20),knee-(10 if s=='L' else 55),ankle,1.5])
        seeds.append(np.clip(alternate,lo+1e-8,hi-1e-8));fits=[]
        for seed in seeds:
            f=least_squares(residual,seed,bounds=(lo,hi),diff_step=1e-5,ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=160);fits.append(f)
            if np.max(abs(residual(f.x)))<1e-5:break
        fit=min(fits,key=lambda f:np.sum(residual(f.x)**2));errors,q,stats=residual(fit.x,True);x=fit.x;qs.append(q);penetration=max([-v[n] for v in stats.values() for n in ('foot_min_gap_m','toe_min_gap_m')]+[0.]);rows.append(dict(phase=phase,max_scaled_IK_error_m=float(abs(errors).max()),penetration_m=penetration,feasible=bool(abs(errors).max()<1e-4 and penetration<1e-4),legs=stats))
    report=dict(machine='easteregg2',source_sha256=anchor['source_sha256'],baseline=str(base),parameters=dict(period_s=period,belt_speed_m_s=speed,stance_duty=duty,swing_clearance_m=height,body_anterior_pitch_deg=pitch_deg,late_stance_heel_up_deg=heel_up_deg,samples=samples),root_height_m=float(qbase[1]),heel_midstance_relative_COM_m=heel_relative,geometric_feasible_samples=sum(v['feasible'] for v in rows),samples=rows,warmed_wall_s=time.perf_counter()-started,accepted_complete_cycle=False,scope='Offline full-cycle geometric prototype only. Nominal original hip/knee/ankle limits unchanged; toe search bounded within its originally unlimited range. Main stance material point follows the actualâˆ’0.2m/s belt, heelâ†’forefoot rocking; swing position/height C2profiles. BodyCOM-x held only in reference generation, no model/controller clamp or simulation. Rootroll/yaw0, wings at original rail height. Toe gravity initialization is NOT a passive dynamic solution; acceleration/contact/Hill/unloading gates remain mandatory. No neural fitting or SNS success claimed.')
    report['root_roll_rule']='Opposite-wing plateaus during single support; smooth transitions during double support. Height from original wing/rail datum.';report['scope']=f'Offline geometric prototype only. Original nominal hip/knee/ankle ROM unchanged; original toe unlimited, search bounded. Actual belt speedâˆ’{speed}m/s is an explicit teacher condition, not neural improvement. Stance material points follow belt; default short step stays nearly flat, optional toe pivot explicitly selected by heel-up parameter. Original opposite-wing roll/contact datum reused, axial yaw0. Swing height has a C2 plateau; ankle uses its unchanged+3degree DF bound while foot pitch is free. Loaded toe clearance50micrometers matches the accepted static anchor contact tolerance, avoiding an artificial exact-zero constraint selecting the wrong passive branch. BodyCOM-x constrained only during reference generation. Toe gravity seed is NOT a dynamic passive solution. Acceleration, contact/CoP, actual Hill-state reachability, true unloading, and full-cycle periodicity still mandatory. No controller/physical-body clamp, neural fitting, or SNS success.'
    folder=root/'periodic_reference';folder.mkdir(exist_ok=True);ledger=folder/'GEOMETRY_ATTEMPTS.json';attempts=json.loads(ledger.read_text()) if ledger.exists() else []
    if (folder/'GEOMETRIC_REFERENCE.json').exists():
        old=json.loads((folder/'GEOMETRIC_REFERENCE.json').read_text());attempts.append({k:old[k] for k in ('parameters','root_height_m','geometric_feasible_samples','warmed_wall_s')});ledger.write_text(json.dumps(attempts,indent=2),encoding='utf-8')
    np.savez_compressed(folder/'GEOMETRIC_REFERENCE.npz',qpos=np.array(qs),phase=np.arange(samples)/samples);(folder/'GEOMETRIC_REFERENCE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report


if __name__=='__main__':
    import argparse
    operations={'register':register_measured_cycle,'validate':validate_measured_dynamics,'close':close_and_audit,'project':project_free_acceleration,'contact-warps':contact_preserving_warps,'motor-shooting':periodic_motor_shooting,'audit-cycle':audit_integrated_cycle,'drift-audit':audit_fixed_wave_drift}
    p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('--operation',choices=['geometry','shared-shooting','shared-shooting-long7',*operations],default='geometry');p.add_argument('--samples',type=int,default=48);a=p.parse_args()
    if a.operation=='geometry':r=geometry(a.root,samples=a.samples)
    elif a.operation=='shared-shooting-long7':
        assert (Path(a.root)/'periodic_reference'/'FIXED_WAVE_CLOSURE_PROTOCOL.json').is_file(),'Require saved preregistered closure protocol'
        r=periodic_motor_shooting(a.root,cycles=7,phase_parameters={'period_scale':1.,'hip_lead_phase':0.,'knee_delay_phase':0.},artifact='SHARED_MOTOR_SHOOTING_LONG7')
    elif a.operation=='shared-shooting':r=periodic_motor_shooting(a.root,phase_parameters={'period_scale':1.,'hip_lead_phase':0.,'knee_delay_phase':0.},artifact='SHARED_MOTOR_SHOOTING')
    else:r=operations[a.operation](a.root)
    print(json.dumps({k:v for k,v in r.items() if k!='samples'},indent=2))
