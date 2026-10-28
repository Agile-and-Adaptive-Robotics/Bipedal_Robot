"""Closed-loop SNS robot with explicit support and reference tension actuators.

SNS membrane outputs drive all 54 actuator controls. Postural torque allocation
sets descending per-MN currents, and Ia/II/Ib plus heel/toe sensors feed the
existing layered spinal circuit. Stepping references are generated from RG
activity, not a time-based sinusoid. This is a hybrid SNS controller, not an
original Deng replication or a calibrated pneumatic hardware simulation.
"""
import os,sys,time,json,argparse,hashlib
from pathlib import Path
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
os.environ.pop('AARL_NET',None)
import numpy as np
import mujoco
from scipy.optimize import lsq_linear

HERE=Path(__file__).resolve().parent
SPINAL=HERE.parent/'spinal'
sys.path.insert(0,str(SPINAL))
sys.path.insert(0,str(SPINAL/'w2l_cpg'))
import build_network as bn
import params
from dynamic_routes import DynamicRoutes
from robot_kinematics import knee_pose,anatomical_knee_tangent
from exact_fast_eval import ExactFastForward


def contact_loads(model,data):
    loads=np.zeros((2,2))
    for index in range(data.ncon):
        contact=data.contact[index]
        names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,int(g)) or '' for g in (contact.geom1,contact.geom2)]
        if 'ground' not in names:continue
        foot=next((n for n in names if n.startswith('foot_')),None)
        if foot is None:continue
        side=0 if foot.startswith('foot_r_') else 1
        wrench=np.zeros(6);mujoco.mj_contactForce(model,data,index,wrench)
        world=contact.frame.reshape(3,3).T@wrench[:3]
        body=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+('r' if side==0 else 'l'))
        local=data.xmat[body].reshape(3,3).T@(contact.pos-data.xpos[body])
        loads[side,int(local[0]>.035)]+=abs(world[1])
    return loads


def rotation_error(quaternion):
    q=np.array(quaternion)
    if q[0]<0:q=-q
    length=np.linalg.norm(q[1:])
    return np.zeros(3) if length<1e-12 else 2*np.arctan2(length,q[0])*q[1:]/length


def allocate(matrix,target,initial,iterations=80,activation_floor=0.,minimum_activation=None):
    # Exact bounded allocation avoids ill-conditioning between large hip
    # capacities and small knee moment arms. Omit welded-joint columns.
    enabled=np.linalg.norm(matrix,axis=0)>1e-7
    lower=np.full(matrix.shape[1],activation_floor) if minimum_activation is None else np.asarray(minimum_activation)
    result=lower.copy()
    active=matrix[:,enabled]
    augmented=np.vstack((active,.03*np.eye(active.shape[1])))
    rhs=np.r_[target,np.zeros(active.shape[1])]
    solution=lsq_linear(augmented,rhs,bounds=(lower[enabled],np.ones(active.shape[1])),method='bvls',tol=1e-7,max_iter=100)
    if solution.status<=0:raise RuntimeError('Muscle force allocation failed: '+solution.message)
    result[enabled]=solution.x
    return result

def motor_target_current(compiled,mn,desired,tau=.005):
    # Engineering MN current clamp: cancel the instantaneous background synaptic
    # current and impose a 30 ms membrane target. RG timing and afferent-central
    # interactions remain active; this is explicitly a hybrid control mode.
    voltage=compiled.V
    gmax=compiled.g_max_non[mn]
    span=compiled.e_hi[mn]-compiled.e_lo[mn]
    conductance=np.clip(gmax*(voltage[None,:]-compiled.e_lo[mn])/np.where(span==0,1.,span),0,gmax)
    synaptic=np.sum(conductance*compiled.del_e[mn],axis=1)-voltage[mn]*np.sum(conductance,axis=1)
    correction=(compiled.dt/tau)/compiled.time_factor_membrane[mn]*(params.E_HI*desired-voltage[mn])
    return compiled.g_m[mn]*(voltage[mn]-compiled.V_rest[mn])-compiled.i_b[mn]-synaptic+correction


def run(args):
    if not 0<=args.activation_floor<1:raise ValueError('Activation floor must be in [0,1).')
    source_paths=(args.model,HERE/'robotbody_routes.json',Path(__file__),HERE/'dynamic_routes.py',
                  SPINAL/'build_network.py',SPINAL/'params.py',SPINAL/'w2l_cpg/exact_fast_eval.py')
    source_hashes={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    snapshot_dir=HERE/'runs/source_snapshots';snapshot_dir.mkdir(parents=True,exist_ok=True)
    runner_bytes=Path(__file__).read_bytes();runner_hash=hashlib.sha256(runner_bytes).hexdigest()
    snapshot=snapshot_dir/(runner_hash+'.py')
    if not snapshot.exists():snapshot.write_bytes(runner_bytes)
    config=json.loads((HERE/'robotbody_routes.json').read_text())
    model=mujoco.MjModel.from_xml_path(str(args.model))
    data=mujoco.MjData(model);dynamic=DynamicRoutes(model,config)
    if args.initial_knee_flexion:
        for side in ('r','l'):
            knee_pose(model,data,side,-args.initial_knee_flexion)
            joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint')
            qa=model.jnt_qposadr[joint];angle=args.initial_knee_flexion/2
            data.qpos[qa:qa+4]=[np.cos(angle/2),0,0,np.sin(angle/2)]
        mujoco.mj_forward(model,data)
        minimum=np.inf
        for gid in range(model.ngeom):
            if not (mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_GEOM,gid) or '').startswith('foot_'):continue
            mid=model.geom_dataid[gid]
            if mid<0:continue
            vertices=model.mesh_vert[model.mesh_vertadr[mid]:model.mesh_vertadr[mid]+model.mesh_vertnum[mid]]
            minimum=min(minimum,float((vertices@data.geom_xmat[gid].reshape(3,3)[1]+data.geom_xpos[gid,1]).min()))
        data.qpos[1]-=minimum+.0002
    dynamic.update(data)
    acts=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
    # Keep production defaults untouched; these settings exist only in this
    # fresh robot-control process and are saved in the run summary.
    params.G['c1_gain']=.15
    params.G['v3_gain']=0.
    params.G['heel_rge']=.15
    params.G['toe_rge']=.05
    params.TAU['rg_nap_h']=args.rg_tau_h
    net=bn.build(acts,dt=model.opt.timestep,interleg=True)
    if args.backend=='exact_fast':net.compiled=ExactFastForward(net.compiled)
    currents=net.make_inputs();ports={name:i for i,name in enumerate(net.inputs)}
    mn=np.array([net.idx['MN_'+name] for name in acts])
    capacities=-model.actuator_gainprm[:,0]
    route_ids=model.actuator_trnid[:,0]
    pelvis=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'pelvis')
    torso=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'lumbar_ball')
    root_height=float(data.qpos[1]);root_x=float(data.qpos[0]);root_z=float(data.qpos[2])
    weight=float(model.body_mass.sum()*9.81)
    q0=data.qpos.copy();target=q0.copy()
    grid=np.linspace(-max(.65,args.knee_amplitude,args.knee_amplitude_right or 0.,args.knee_amplitude_left or 0.),.15,111);knee_lookup={}
    for side in ('r','l'):
        joints=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,x+'_'+side+'_joint') for x in ('back','shank','front')]
        addresses=model.jnt_qposadr[joints]
        values=[];seed=np.zeros(3)
        for angle in grid:
            scratch=mujoco.MjData(model)
            try:seed=knee_pose(model,scratch,side,float(angle),seed)
            except (ValueError,np.linalg.LinAlgError):seed=knee_pose(model,scratch,side,float(angle))
            values.append(seed.copy())
        knee_lookup[side]=(addresses,np.array(values))
    kp=np.full(model.nv,args.hip_kp);kd=np.full(model.nv,args.hip_kd)
    for side in ('r','l'):
        for prefix in ('back','shank','front'):
            j=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,prefix+'_'+side+'_joint')
            kp[model.jnt_dofadr[j]]=args.knee_kp;kd[model.jnt_dofadr[j]]=args.knee_kd
    address=model.jnt_dofadr[torso];kp[address:address+3]=65.;kd[address:address+3]=10.
    desired_activation=np.zeros(model.nu);integral=np.zeros(model.nu)
    length0=data.ten_length[route_ids].copy()
    error=np.zeros(model.nv);phase=np.zeros(2)
    rg_swing=np.zeros(2,dtype=bool);event_time=np.zeros(2)
    rhythm_scale=args.rg_tau_h/.35
    flex_duration=np.full(2,.45*rhythm_scale);extension_duration=np.full(2,.95*rhythm_scale)
    log={name:[] for name in ('time','qpos','qvel','activation','control','desired_activation','forces','rig_force','rig_torque','contact_loads','rg','knee_coordinate','route_lengths','active_counts')}
    start=time.perf_counter();failed=None
    previous_vertical_force=args.support*weight
    foot_ids=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+side) for side in ('r','l')]
    hip_dofs=[int(model.jnt_dofadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint')]) for side in ('r','l')]
    steps=round(args.duration/model.opt.timestep)
    for step in range(steps):
        t=step*model.opt.timestep
        if args.disable_muscles:data.act[:]=0
        dynamic.update(data)
        loads=contact_loads(model,data)
        v=net.compiled.V
        walking=args.mode=='step' and t>args.stand_time
        if walking:
            for si,side in enumerate(('r','l')):
                difference=float(v[net.idx['RG_F_'+side]]-v[net.idx['RG_E_'+side]])
                requested=float(np.clip(difference/5.,0.,1.))
                phase[si]+=(requested-phase[si])*model.opt.timestep/.1
                hip=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint')
                hip_fraction=phase[si];knee_fraction=phase[si]
                if args.pattern=='event':
                    swing=(difference>-.5 if rg_swing[si] else difference>.5)
                    if swing!=rg_swing[si]:
                        elapsed=t-event_time[si]
                        if .15<elapsed<max(2.5,2.5*rhythm_scale):
                            if rg_swing[si]:flex_duration[si]=.5*(flex_duration[si]+elapsed)
                            else:extension_duration[si]=.5*(extension_duration[si]+elapsed)
                        rg_swing[si]=swing;event_time[si]=t
                    elapsed=t-event_time[si]
                    if swing:
                        u=np.clip(elapsed/flex_duration[si],0.,1.)
                        hip_fraction=3*u*u-2*u*u*u
                        knee_fraction=np.sin(np.pi*u)
                    else:
                        u=np.clip(elapsed/extension_duration[si],0.,1.)
                        hip_fraction=1-(3*u*u-2*u*u*u)
                        knee_fraction=0.
                hip_angle=args.hip_flexion*hip_fraction+args.hip_extension*(1-hip_fraction)
                qa=model.jnt_qposadr[hip]
                target[qa:qa+4]=[np.cos(hip_angle/2),0,0,np.sin(hip_angle/2)]
                amplitude=(args.knee_amplitude_right if side=='r' else args.knee_amplitude_left)
                if amplitude is None:amplitude=args.knee_amplitude
                angle=-amplitude*knee_fraction
                addresses,table=knee_lookup[side]
                target[addresses]=[np.interp(angle,grid,table[:,i]) for i in range(3)]
        if step%5==0:
            mujoco.mj_differentiatePos(model,error,1.,data.qpos,target)
            rows=np.flatnonzero(data.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
            jacobian=data.efc_J.reshape(data.nefc,model.nv)[rows,6:]
            _,singular,vh=np.linalg.svd(jacobian,full_matrices=True)
            rank=int(np.sum(singular>1e-8));null=vh[rank:].T
            if args.anatomical_knee_tangent:null=anatomical_knee_tangent(model,data)
            desired_torque=data.qfrc_bias+kp*error-kd*data.qvel
            if args.abduction_gain>0:
                pelvis_right=data.xmat[pelvis].reshape(3,3)[:,2]
                separation=float(pelvis_right@(data.xpos[foot_ids[0]]-data.xpos[foot_ids[1]]))
                jacobians=[]
                for bid in foot_ids:
                    jac=np.zeros((3,model.nv));angular=np.zeros_like(jac)
                    mujoco.mj_jacBody(model,data,jac,angular,bid);jacobians.append(jac)
                gradient=pelvis_right@(jacobians[0]-jacobians[1]);gradient[:6]=0.
                rate=float(gradient@data.qvel)
                separating_force=np.clip(args.abduction_gain*(args.stance_width-separation)-args.abduction_damping*rate,
                    -args.abduction_force_limit,args.abduction_force_limit)
                # Virtual separation feedback enters muscle torque allocation only.
                # No external hip force, qpos correction, or rig motion is added.
                desired_torque+=gradient*separating_force
            if args.foot_heading_gain>0:
                R=data.xmat[pelvis].reshape(3,3).T
                for bid in foot_ids:
                    direction=R@data.xmat[bid].reshape(3,3)[:,0]
                    denominator=max(float(direction[0]**2+direction[2]**2),1e-8)
                    heading=float(np.arctan2(direction[2],direction[0]))
                    angular_gradient=np.array([direction[0]*direction[1]/denominator,-1.,direction[2]*direction[1]/denominator])
                    jac=np.zeros((3,model.nv));angular=np.zeros_like(jac)
                    mujoco.mj_jacBody(model,data,jac,angular,bid)
                    gradient=angular_gradient@R@angular;gradient[:6]=0.
                    torque=np.clip(-args.foot_heading_gain*heading-args.foot_heading_damping*(gradient@data.qvel),
                        -args.foot_heading_torque_limit,args.foot_heading_torque_limit)
                    desired_torque+=gradient*torque
            if args.stable_pd_horizon:
                full_mass=np.zeros((model.nv,model.nv));mujoco.mj_fullM(model,full_mass,data.qM)
                mass=null.T@full_mass[6:,6:]@null
                stiffness=null.T@(kp[6:,None]*null);damping=null.T@(kd[6:,None]*null)
                velocity=null.T@data.qvel[6:];position_error=null.T@error[6:]
                h=args.stable_pd_horizon
                acceleration=np.linalg.solve(mass+h*damping+h*h*stiffness,
                    stiffness@position_error-(damping+h*stiffness)@velocity)
                stable=null@(mass@acceleration)
                original_trunk=desired_torque[address:address+3].copy()
                desired_torque[6:]=data.qfrc_bias[6:]+stable
                desired_torque[address:address+3]=original_trunk
            if args.load_compensation:
                # Static equilibrium: tau + J_foot.T F_ground - bias = 0.
                # Transfer the unsupported body load through stance muscles;
                # the rig remains an external pelvis force, unchanged here.
                shares=np.maximum(1.-phase,.05) if walking else np.ones(2)
                if walking and args.pattern=='event' and args.unload_swing:
                    shares=(~rg_swing).astype(float)
                    if shares.sum()==0:shares=np.ones(2)
                shares/=shares.sum()
                unsupported=max(0.,weight-previous_vertical_force)
                for si,bid in enumerate(foot_ids):
                    jac=np.zeros((3,model.nv));angular=np.zeros_like(jac)
                    mujoco.mj_jacBodyCom(model,data,jac,angular,bid)
                    desired_torque-=args.load_compensation*min(1.,t/2.)*unsupported*shares[si]*jac[1]
            matrix=null.T@dynamic.moment_arms[route_ids,6:].T*capacities[None,:]
            reduced_target=null.T@desired_torque[6:]
            weights=np.ones(model.nv-6)
            for dof in hip_dofs:weights[dof-6:dof-4]=args.hip_frontal_weight
            weighting=null.T@(weights[:,None]*null)
            lower=np.full(model.nu,args.activation_floor)
            if args.abductor_floor>args.activation_floor:
                for dof,sign in zip(hip_dofs,(-1.,1.)):
                    abducting=sign*dynamic.moment_arms[route_ids,dof]>.002
                    lower[abducting]=args.abductor_floor
            desired_activation=allocate(weighting@matrix,weighting@reduced_target,desired_activation,
                activation_floor=args.activation_floor,minimum_activation=lower)
        if step==0:
            net.compiled.V[mn]=params.E_HI*desired_activation
            for _ in range(100):
                actual=np.clip(net.compiled.V[mn]/params.E_HI,0.,1.)
                integral=np.clip(integral+(desired_activation-actual)*model.opt.timestep,-1.,1.)
                warm_post=motor_target_current(net.compiled,mn,desired_activation,args.motor_tau)
                currents[:]=0
                for i,name in enumerate(acts):currents[ports['POST_'+name]]=float(warm_post[i])
                currents[ports['POSTURE']]=.15
                net.step(currents)
            data.act[:]=np.clip(net.compiled.V[mn]/params.E_HI,0.,1.)
            if args.disable_muscles:data.act[:]=0
        if args.mode=='step' and step==round(args.stand_time/model.opt.timestep):
            # Antiphase kickoff once; cycle timing thereafter comes from RGs.
            net.compiled.V[net.idx['RG_F_r']]=params.E_HI
            net.compiled.V[net.idx['RG_E_l']]=params.E_HI
        actual_activation=np.clip(net.compiled.V[mn]/params.E_HI,0.,1.)
        currents[:]=0
        integral=np.clip(integral+(desired_activation-actual_activation)*model.opt.timestep,-1.,1.)
        motor_request=np.clip(desired_activation+args.activation_tracking_gain*(desired_activation-data.act),args.activation_floor,1.)
        post=motor_target_current(net.compiled,mn,motor_request,args.motor_tau)
        for i,name in enumerate(acts):
            currents[ports['POST_'+name]]=float(post[i])
            tid=route_ids[i]
            currents[ports['Ia_'+name]]=max(0.,.3*dynamic.velocities[tid]/.3)
            currents[ports['II_'+name]]=max(0.,.3*(data.ten_length[tid]-length0[i])/max(.1*length0[i],.01))
            currents[ports['Ib_'+name]]=max(0.,.15*-data.actuator_force[i]/capacities[i])
        currents[ports['DRIVE']]=args.walk_drive if walking else 0.
        currents[ports['POSTURE']]=.15
        for si,side in enumerate(('r','l')):
            for part,prefix in enumerate(('HEEL_c_','TOE_c_')):
                name=prefix+side
                if name in ports:currents[ports[name]]=min(2.,loads[si,part]/200.)
        net.step(currents)
        data.ctrl[:]=np.clip(net.compiled.V[mn]/params.E_HI,0.,1.)
        if args.diagnostic_direct:data.ctrl[:]=desired_activation
        if args.disable_muscles:data.ctrl[:]=0
        # Support is an explicit external pelvis rig. It never applies joint
        # torques or cancels gravity on the muscles/limbs.
        lateral_target=root_z+args.lateral_shift*(phase[1]-phase[0]) if walking else root_z
        force=np.array([-args.horizontal_stiffness*(data.qpos[0]-root_x)-30*data.qvel[0],
            args.support*weight+args.vertical_stiffness*(root_height-data.qpos[1])-120*data.qvel[1],
            -args.lateral_stiffness*(data.qpos[2]-lateral_target)-args.lateral_damping*data.qvel[2]])
        support_cap=args.support_cap
        if args.support_cap_final is not None:
            fraction=np.clip((t-args.support_wean_start)/args.support_wean_duration,0.,1.)
            support_cap+=(args.support_cap_final-support_cap)*fraction
        force[1]=np.clip(force[1],0,support_cap*weight)
        previous_vertical_force=force[1]
        torque=-args.rig_angular_stiffness*rotation_error(data.qpos[3:7])-args.rig_angular_damping*data.qvel[3:6]
        data.xfrc_applied[:]=0;data.xfrc_applied[pelvis,:3]=force;data.xfrc_applied[pelvis,3:]=torque
        data.qfrc_applied[:]=0
        dynamic.apply_force_correction(data)
        if args.frontal_assist:
            for side in ('r','l'):
                joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint')
                qa=model.jnt_qposadr[joint];da=model.jnt_dofadr[joint]
                hip_error=rotation_error(data.qpos[qa:qa+4])
                data.qfrc_applied[da:da+2]+=-args.frontal_assist*hip_error[:2]-8*data.qvel[da:da+2]
        if step%5==0:
            rg=[net.compiled.V[net.idx[x+'_'+s]] for s in ('r','l') for x in ('RG_E','RG_F')]
            values=(t,data.qpos.copy(),data.qvel.copy(),data.act.copy(),data.ctrl.copy(),desired_activation.copy(),data.actuator_force.copy(),
                    force,torque,loads,rg,list(dynamic.coordinates(data).values()),data.ten_length.copy(),dynamic.active_counts.copy())
            for key,value in zip(log,values):log[key].append(value)
        mujoco.mj_step(model,data)
        if not np.isfinite(data.qpos).all() or not np.isfinite(net.compiled.V).all():failed='nonfinite state';break
        if any(w.number for w in data.warning):failed='MuJoCo numerical warning';break
        if data.qpos[1]<root_height-.3:failed='pelvis height failure';break
        if np.linalg.norm(rotation_error(data.qpos[3:7]))>np.deg2rad(60):failed='pelvis tilt failure';break
        if np.linalg.norm(rotation_error(data.qpos[model.jnt_qposadr[torso]:model.jnt_qposadr[torso]+4]))>np.deg2rad(60):failed='lumbar posture failure';break
    arrays={key:np.array(value) for key,value in log.items()}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(args.output.with_suffix('.npz'),**arrays)
    trunk_address=model.jnt_qposadr[torso]
    trunk_angles=[np.linalg.norm(rotation_error(q[trunk_address:trunk_address+4])) for q in arrays['qpos']]
    summary=dict(machine='easteregg2',model_path=str(args.model),mode=args.mode,backend=args.backend,
        duration_simulated_s=float(data.time),requested_duration_s=args.duration,wall_s=time.perf_counter()-start,
        failed=failed,finite=failed not in ('nonfinite state','MuJoCo numerical warning'),
        neuron_count=len(net.idx),neural_inputs=len(net.inputs),actuators=model.nu,
        active_force_actuators=int(np.sum(np.max(np.abs(arrays['forces']),axis=0)>1)),
        actuators_with_joint_leverage=int(np.sum(np.linalg.norm(matrix,axis=0)>1e-7)),
        maximum_control=float(arrays['control'].max()),pelvis_height_min_m=float(arrays['qpos'][:,1].min()),
        mean_neural_tracking_error=float(np.mean(np.abs(arrays['control']-arrays['desired_activation']))),
        maximum_lumbar_rotation_deg=float(np.rad2deg(max(trunk_angles))),
        pelvis_forward_displacement_m=float(arrays['qpos'][-1,0]-arrays['qpos'][0,0]),
        knee_range_deg=np.rad2deg(np.ptp(arrays['knee_coordinate'],axis=0)).tolist(),
        mean_support_fraction=float(arrays['rig_force'][:,1].mean()/weight),nominal_support_fraction=args.support,
        mean_rig_torque_Nm=np.mean(np.abs(arrays['rig_torque']),axis=0).tolist(),
        mean_ground_load_N=np.mean(arrays['contact_loads'].sum(axis=2),axis=0).tolist(),
        diagnostic_direct_control=args.diagnostic_direct,motor_target_tau_s=args.motor_tau,
        rig_angular_stiffness_Nm_rad=args.rig_angular_stiffness,rig_angular_damping_Nms_rad=args.rig_angular_damping,
        frontal_hip_guide_stiffness_Nm_rad=args.frontal_assist,
        muscles_disabled=args.disable_muscles,vertical_rig_stiffness_N_m=args.vertical_stiffness,
        run_options={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
        source_hashes=source_hashes,runner_snapshot=str(snapshot),
        effective_neural_gains=dict(params.G),
        effective_neural_time_constants=dict(params.TAU),
        neural_parameters=dict(c1_gain=.15,v3_gain=0,heel_rge=.15,toe_rge=.05),
        force_model='30 ms first-order reference tension, capped by existing 27-actuator map MIFs; not calibrated BPA hardware',
        controller='Layered SNS with Ia/II/Ib/contact inputs, descending postural torque allocation, compensated MN target currents, and RG-driven stepping references',
        status='Supported neural-control diagnostic; walking classification requires trajectory/contact review')
    args.output.with_suffix('.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model',type=Path,default=HERE/'ah001_neural.xml')
    parser.add_argument('--initial-knee-flexion',type=float,default=0.,help='Explicit initial knee flexion, rad, with half-angle hip flexion and soles placed at ground')
    parser.add_argument('--anatomical-knee-tangent',action='store_true',help='Preserve hip authority when projecting planar knee-loop constraints')
    parser.add_argument('--mode',choices=['stand','step'],default='stand')
    parser.add_argument('--duration',type=float,default=10.)
    parser.add_argument('--stand-time',type=float,default=3.)
    parser.add_argument('--support',type=float,default=.8)
    parser.add_argument('--horizontal-stiffness',type=float,default=20.)
    parser.add_argument('--lateral-shift',type=float,default=0.,help='Explicit RG-timed lateral pelvis support target amplitude in metres')
    parser.add_argument('--lateral-stiffness',type=float,default=250.)
    parser.add_argument('--lateral-damping',type=float,default=60.)
    parser.add_argument('--knee-amplitude',type=float,default=.3)
    parser.add_argument('--knee-amplitude-right',type=float)
    parser.add_argument('--knee-amplitude-left',type=float)
    parser.add_argument('--hip-flexion',type=float,default=.18)
    parser.add_argument('--hip-extension',type=float,default=-.04)
    parser.add_argument('--pattern',choices=['direct','event'],default='direct',help='Direct RG mapping or RG-event-anchored swing/stance trajectories')
    parser.add_argument('--unload-swing',action='store_true',help='Allocate ground-load feedforward only to RG stance legs')
    parser.add_argument('--activation-floor',type=float,default=0.,help='Minimum allocated activation, including opposing muscles; compare .05 with zero')
    parser.add_argument('--abductor-floor',type=float,default=0.,help='Minimum activation on routes currently producing hip abduction, measured from their Jacobians')
    parser.add_argument('--abduction-gain',type=float,default=0.,help='Muscle-mediated ankle-separation feedback, N/m')
    parser.add_argument('--abduction-damping',type=float,default=8.)
    parser.add_argument('--abduction-force-limit',type=float,default=15.)
    parser.add_argument('--stance-width',type=float,default=.18)
    parser.add_argument('--hip-frontal-weight',type=float,default=1.,help='Allocation priority for hip frontal and axial torque accuracy')
    parser.add_argument('--foot-heading-gain',type=float,default=0.,help='Muscle-mediated foot-heading feedback, Nm/rad')
    parser.add_argument('--foot-heading-damping',type=float,default=1.)
    parser.add_argument('--foot-heading-torque-limit',type=float,default=10.)
    parser.add_argument('--walk-drive',type=float,default=3.)
    parser.add_argument('--rg-tau-h',type=float,default=.35,help='Intrinsic neural burst-termination time constant; changes cadence without a time-based oscillator')
    parser.add_argument('--hip-kp',type=float,default=25.)
    parser.add_argument('--hip-kd',type=float,default=5.)
    parser.add_argument('--knee-kp',type=float,default=10.)
    parser.add_argument('--knee-kd',type=float,default=1.5)
    parser.add_argument('--load-compensation',type=float,default=0.,help='Stance-load feedforward gain through muscle allocation; no extra rig force')
    parser.add_argument('--backend',choices=['reference','exact_fast'],default='exact_fast')
    parser.add_argument('--motor-tau',type=float,default=.005)
    parser.add_argument('--activation-tracking-gain',type=float,default=0.,help='Closed-loop activation tracking through MN current; physical activation filter remains unchanged')
    parser.add_argument('--stable-pd-horizon',type=float,default=0.,help='Optional implicit posture prediction horizon on the knee-loop constraint tangent space')
    parser.add_argument('--rig-angular-stiffness',type=float,default=300.)
    parser.add_argument('--rig-angular-damping',type=float,default=60.)
    parser.add_argument('--frontal-assist',type=float,default=0.,help='Optional external hip frontal/yaw posture guide; explicitly reported')
    parser.add_argument('--disable-muscles',action='store_true',help='Ablation diagnostic; do not count as neural-control success')
    parser.add_argument('--vertical-stiffness',type=float,default=800.)
    parser.add_argument('--support-cap',type=float,default=1.5,help='Hard cap on total actual vertical rig force, as fraction of body weight')
    parser.add_argument('--support-cap-final',type=float)
    parser.add_argument('--support-wean-start',type=float,default=5.)
    parser.add_argument('--support-wean-duration',type=float,default=10.)
    parser.add_argument('--diagnostic-direct',action='store_true',help='Diagnostic only: bypass neural motor output; never count as neural control success')
    parser.add_argument('--output',type=Path,default=HERE/'runs/neural_stand')
    run(parser.parse_args())
