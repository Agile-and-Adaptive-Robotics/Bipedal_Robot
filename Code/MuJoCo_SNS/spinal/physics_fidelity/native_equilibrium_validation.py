"""Five-second source-muscle endpoint validation,without live state clamps.

This reusable check distinguishes short bounded holds from stationary native
equilibria. It keeps every native coordinate,contact,and follower constraint.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
from time import perf_counter
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent),str(HERE.parent/'functional_gait2392')]
import mujoco
import numpy as np
from mechanics import CoupledArms
from contact_standing_solver import _external,_support_regions
from thelen_dynamic import ConstitutiveBackend,DynamicThelenBank


def validate(model_folder,endpoint,spec_path,out,allow_unaccepted_static=False,physics_dt=None,initial_audit_only=False):
    code_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if (out/'STATIONARITY_RESULT.json').exists():
        raise ValueError('Completed validation exists; preserve it and use a new folder')
    source=json.loads((endpoint/'SOURCE_NATIVE_STANCE.json').read_text())
    if not source['accepted_static_source_recruitment'] and not allow_unaccepted_static:
        raise ValueError('Static gates failed; an explicit diagnostic override is required')
    cfg=json.loads((Path(source['source_pose'])/'MANIFEST.json').read_text())
    with np.load(endpoint/'source_recruitment.npz') as z:
        q=z['qpos'].copy();alpha=z['activation'].copy();p=z['parameters'].copy()
        names=z['muscles'].tolist();enabled=z.get('force_channel_enabled',np.ones(len(names),bool))
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'))
    if physics_dt is not None:
        if not np.isfinite(physics_dt) or physics_dt<=0.:
            raise ValueError('Positive finite diagnostic physics timestep required')
        m.opt.timestep=physics_dt
    if not np.allclose(m.opt.gravity/np.linalg.norm(m.opt.gravity),[0.,0.,-1.]):
        raise ValueError('Explicit captured z-up human model required')
    d=mujoco.MjData(m);ar=CoupledArms(m);dt=float(m.opt.timestep)
    bw=float(m.body_mass.sum()*9.81);support=float(cfg['support_fraction'])
    d.qpos[:]=q;d.act[:]=alpha;d.ctrl[:]=alpha
    d.xfrc_applied[m.body('pelvis').id,2]=support*bw
    mujoco.mj_forward(m,d)
    backend=ConstitutiveBackend(m,d,p,mode='compliant')
    if not enabled.all():
        class MaskedBank(DynamicThelenBank):
            def evaluate(self,length,activation,dynamic=True):
                super().evaluate(length,activation,dynamic)
                self.force[~enabled]=0.;self.tangent[~enabled]=0.
                self.error[~enabled]=0.;self.bounds[~enabled]=0
                return self.force
        backend.bank=MaskedBank(p,dt)
        backend.bank.initialize(d.actuator_length,alpha)
    m.actuator_biasprm[:,0]=-backend.bank.force
    mujoco.mj_forward(m,d)
    probe=mujoco.MjData(m)
    static_bank=DynamicThelenBank(p,dt)
    mass=np.empty((m.nv,m.nv));mujoco.mj_fullM(m,mass,d.qM)
    probe.qpos[:]=d.qpos;probe.qvel[:]=d.qvel;probe.qacc[:]=d.qacc
    probe.act[:]=d.act;probe.ctrl[:]=d.ctrl;probe.xfrc_applied[:]=d.xfrc_applied
    mujoco.mj_inverse(m,probe)
    forward_inverse=d.qfrc_actuator+_external(m,d)-probe.qfrc_inverse
    momentum=mass@d.qacc-d.qfrc_smooth-d.qfrc_constraint
    initial_forward=dict(source_endpoint=str(endpoint),source_code_sha256=code_hash,full_acceleration_max_native=float(np.max(abs(d.qacc))),root_acceleration_native=d.qacc[:6].tolist(),anatomical_joint_names=[m.joint(int(j)).name for j in ar.jids],anatomical_acceleration_deg_s2=np.rad2deg(d.qacc[ar.rows]).tolist(),full_forward_inverse_force_residual_max=float(np.max(abs(forward_inverse))),full_momentum_balance_max=float(np.max(abs(momentum))),source_actuator_generalized_force_error_max=float(np.max(abs(d.qfrc_actuator+d.actuator_moment.T@backend.bank.force))),protocol_mass_sum_kg=float(m.body_mass.sum()),supported_pelvis_subtree_mass_kg=float(m.body_subtreemass[m.body('pelvis').id]),prescribed_support_protocol_BW=support,scope='Read-only initial native forward/inverse and full105 momentum checks after installing unchanged source static forces. No live state change,force fit,or physics modification.')
    if initial_audit_only:
        out.mkdir(parents=True,exist_ok=True);target=out/('INITIAL_NATIVE_AUDIT_'+endpoint.name+'.json')
        if target.exists():raise ValueError('Preserve existing initial native audit')
        target.write_text(json.dumps(initial_forward,indent=2));print(json.dumps(initial_forward,indent=2));return initial_forward

    def stationary_residual(position,force=None):
        probe.qpos[:]=position;probe.qvel[:]=0.;probe.qacc[:]=0.
        probe.act[:]=alpha;probe.ctrl[:]=alpha
        probe.xfrc_applied[:]=d.xfrc_applied
        mujoco.mj_inverse(m,probe)
        if force is None:
            static_bank.initialize(probe.actuator_length,alpha)
            force=static_bank.force.copy();force[~enabled]=0.
        residual=-probe.actuator_moment.T@force-(probe.qfrc_inverse-_external(m,probe))
        return residual

    initial_raw=stationary_residual(q)
    initial_regions,_=_support_regions(m,d)
    source_air=np.array([min(r['gap_m'] for r in initial_regions[2*s:2*s+2]) for s in range(2)])>.005
    spec=json.loads(spec_path.read_text())
    addresses=np.array([row['qpos_address'] for row in spec['joints'].values()])
    low=np.array([row['native_bounds_rad'][0] for row in spec['joints'].values()])
    high=np.array([row['native_bounds_rad'][1] for row in spec['joints'].values()])
    n=round(5./dt);late_n=round(.5/dt)
    positions=np.empty((n,m.nq));velocities=np.empty((n,m.nv))
    forces=np.empty((n,m.nu));fiber_length=np.empty_like(forces);fiber_velocity=np.empty_like(forces)
    gaps=np.empty((n,2));normal=np.zeros((n,2));contact_force=np.zeros(6)
    guard='';maximum_activation_error=0.;started=perf_counter()
    for i in range(n):
        d.ctrl[:]=alpha
        try:
            backend.step()
        except RuntimeError as error:
            guard=str(error)
        positions[i]=d.qpos;velocities[i]=d.qvel;forces[i]=backend.bank.force
        length=np.hypot(backend.bank.state,p[:,3])
        fiber_length[i]=length/p[:,1]
        fiber_velocity[i]=(backend.bank.state-backend.old_state)*backend.bank.state/(length*dt*p[:,1])
        regions,_=_support_regions(m,d)
        gaps[i]=[min(r['gap_m'] for r in regions[2*s:2*s+2]) for s in range(2)]
        for ci in range(d.ncon):
            c=d.contact[ci];a,b=int(c.geom1),int(c.geom2)
            if m.geom_type[a]==mujoco.mjtGeom.mjGEOM_PLANE:g=b
            elif m.geom_type[b]==mujoco.mjtGeom.mjGEOM_PLANE:g=a
            else:continue
            name=m.body(int(m.geom_bodyid[g])).name
            if not name.endswith(('_r','_l')):continue
            mujoco.mj_contactForce(m,d,ci,contact_force)
            normal[i,0 if name.endswith('_r') else 1]+=contact_force[0]/bw
        maximum_activation_error=max(maximum_activation_error,float(np.max(abs(d.act-alpha))))
        if np.any(d.qpos[addresses]<low-1e-10) or np.any(d.qpos[addresses]>high+1e-10):
            guard=guard or 'Frozen source ROM'
        if guard:break
    elapsed=perf_counter()-started;count=i+1;late=slice(max(0,count-late_n),count)
    positions=positions[:count];velocities=velocities[:count];forces=forces[:count]
    fiber_length=fiber_length[:count];fiber_velocity=fiber_velocity[:count]
    gaps=gaps[:count];normal=normal[:count]
    final_dynamic_raw=stationary_residual(d.qpos,backend.bank.force)
    final_static_raw=stationary_residual(d.qpos)
    projected=mujoco.MjData(m);projected.qpos[:]=d.qpos;ar.project(projected)
    depaddr=np.array([m.jnt_qposadr[m.dof_jntid[k]] for k in ar.dependent])
    joint_speed=float(np.max(abs(velocities[late][:,ar.rows])))
    joint_drift=float(np.max(abs(positions[:,ar.addresses]-q[ar.addresses])))
    root_drift=float(np.max(abs(positions[:,:3]-q[:3])))
    root_rotation=float(np.max(abs(positions[:,3:6]-q[3:6])))
    fiber_speed=float(np.max(abs(fiber_velocity[late][:,enabled])))
    force_variation=float(np.max(np.ptp(forces[late],axis=0)/p[:,0]))
    contact_ok=bool(gaps.min()>=-.002)
    for side in range(2):
        contact_ok=contact_ok and bool((gaps[:,side].min()>=.01 and normal[:,side].max()<=.025) if source_air[side] else normal[late,side].min()>=.05)
    physical_ok=bool(not guard and count==n and root_drift<=.0005 and root_rotation<=np.deg2rad(.25) and joint_drift<=np.deg2rad(.25) and joint_speed<=np.deg2rad(.05) and np.max(abs(velocities[late,:3]))<=1e-4 and np.max(abs(velocities[late,3:6]))<=np.deg2rad(.05) and fiber_speed<=.001 and force_variation<=1e-4 and np.max(abs(final_dynamic_raw))<=.05 and np.max(abs(final_static_raw))<=.05 and maximum_activation_error<=1e-12 and contact_ok)
    load_reference=cfg.get('stationary_load_split_reference')
    load_error=None
    if load_reference is not None:
        load_error=float(np.max(abs(normal-np.asarray(load_reference['normal_target_BW']))))
        physical_ok=physical_ok and load_error<=.002
    metrics=dict(duration_s=float(d.time),guard=guard,wall_s=elapsed,initial_full_static_residual_max=float(np.max(abs(initial_raw))),final_full_dynamic_force_static_residual_max=float(np.max(abs(final_dynamic_raw))),final_full_static_force_residual_max=float(np.max(abs(final_static_raw))),final_follower_force_residual_max=float(np.max(abs(final_static_raw[ar.dependent]))),final_follower_position_offset_max_m=float(np.max(abs(d.qpos[depaddr]-projected.qpos[depaddr]))),root_translation_drift_max_m=root_drift,root_rotation_drift_max_deg=float(np.rad2deg(root_rotation)),joint_drift_max_deg=float(np.rad2deg(joint_drift)),last_half_joint_velocity_max_deg_s=float(np.rad2deg(joint_speed)),last_half_joint_velocity_deg_s=np.rad2deg(np.max(abs(velocities[late][:,ar.rows]),axis=0)).tolist(),last_half_follower_velocity_max_mps=float(np.max(abs(velocities[late][:,ar.dependent]))),last_half_fiber_velocity_max_l0_s=fiber_speed,last_half_force_variation_max_F0=force_variation,minimum_foot_mesh_gaps_m=gaps.min(axis=0).tolist(),last_half_normal_min_BW=normal[late].min(axis=0).tolist(),maximum_activation_change=maximum_activation_error,constitutive_backend=backend.diagnostics(),physical_stationarity_gates_pass=physical_ok)
    accepted=bool(source['accepted_static_source_recruitment'] and physical_ok)
    metrics['initial_native_forward']=initial_forward
    if load_reference is not None:
        metrics['normal_target_BW']=load_reference['normal_target_BW'];metrics['normal_target_error_max_BW']=load_error
    out.mkdir(parents=True,exist_ok=True)
    manifest=dict(format='aarl_full_native_stationarity_validation_v2',source_endpoint=str(endpoint),source_model=str(model_folder/'runtime_model.mjb'),support_fraction=support,source_code_sha256=code_hash,source_controls_sha256=hashlib.sha256((endpoint/'source_recruitment.npz').read_bytes()).hexdigest(),duration_s=5.,timestep_s=dt,allow_unaccepted_static_diagnostic=allow_unaccepted_static,scope='Constant saved source activation/excitation,all native coordinates and constraints retained; no live state clamps or control retuning.')
    if load_reference is not None:manifest['stationary_load_split_reference']=load_reference
    (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2))
    result=dict(accepted_stationary_endpoint=accepted,accepted_static_source_recruitment=source['accepted_static_source_recruitment'],source_pose=str(out),metrics=metrics,manifest=manifest)
    (out/'STATIONARITY_RESULT.json').write_text(json.dumps(result,indent=2))
    (out/'SOURCE_NATIVE_STANCE.json').write_text(json.dumps(result,indent=2))
    shutil.copyfile(endpoint/'source_recruitment.npz',out/'source_recruitment.npz')
    stride=max(1,round(.01/dt));indices=np.unique(np.r_[np.arange(0,count,stride),count-1])
    np.savez_compressed(out/'stationarity_sparse.npz',t=(indices+1)*dt,qpos=positions[indices],qvel=velocities[indices],foot_mesh_gaps_m=gaps[indices],normal_side_BW=normal[indices],late_t=np.arange(count-len(forces[late]),count)*dt+dt,late_force_N=forces[late],late_fiber_length_l0=fiber_length[late],late_fiber_velocity_l0_s=fiber_velocity[late],late_joint_velocity=velocities[late][:,ar.rows])
    print(json.dumps(dict(accepted_stationary_endpoint=accepted,metrics=metrics),indent=2))
    return result


def summarize_stationary_campaign(model_folder,folder,spec_path):
    """Curate repeatability and failed transport evidence without refitting physics."""
    def numeric_parity(first,second):
        with np.load(first) as a,np.load(second) as b:
            return {k:float(np.max(abs(a[k].astype(np.result_type(a[k].dtype,np.float64))-b[k]))) if a[k].size else 0.
                    for k in a.files if a[k].dtype.kind in 'biufc' and k in b.files and a[k].shape==b[k].shape}
    spec=json.loads(spec_path.read_text());addresses=np.array([v['qpos_address'] for v in spec['joints'].values()])
    model=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));ar=CoupledArms(model)
    result=dict(format='aarl_stationary_loadsplit_checkpoint_v1',support_fraction=.6,
                source_actuation_policy='Original source92,including six legacy-pruned force channels',
                source_engine_sha256=hashlib.sha256((HERE/'thelen_dynamic.py').read_bytes()).hexdigest(),
                scope='Stationary references and separate failed full-source command transport; no neural fitting or walking claim.',sides={})
    baseline=folder.parent/'native_equilibrium_v2/double_support_source92_jacobian_scaled_5s'
    base=np.load(baseline/'source_recruitment.npz');q0=base['qpos']
    muscles=base['muscles'].tolist();families=sorted({name[:-2] for name in muscles})
    phase_basis=[];source_jacobians=[]
    default=folder/'defaultoff_oldDS/stationarity_sparse.npz'
    if default.exists():result['old_DS_defaultoff_numeric_parity']=numeric_parity(baseline/'stationarity_sparse.npz',default)
    for side in ('right','left'):
        reports={key:json.loads((folder/(side+'_'+key)/'STATIONARITY_RESULT.json').read_text())
                 for key in ('stationary_5s','stationary_fresh','stationary_1ms')}
        parity=numeric_parity(folder/(side+'_stationary_5s')/'stationarity_sparse.npz',folder/(side+'_stationary_fresh')/'stationarity_sparse.npz')
        row=dict(accepted_stationary_reference=bool(all(r['accepted_stationary_endpoint'] for r in reports.values()) and max(parity.values())==0.),
                 reference=str(folder/(side+'_stationary_5s')),fresh_numeric_errors=parity,
                 validations={k:dict(accepted=r['accepted_stationary_endpoint'],metrics=r['metrics']) for k,r in reports.items()})
        motion=folder/(side+'_direct_C2_motion');transition=json.loads((motion/'TRANSITION_RESULT.json').read_text())
        z=np.load(motion/'training_reference.npz');target=np.load(folder/(side+'_stationary_5s')/'source_recruitment.npz')
        q=z['qpos'];v=z['qvel'];times=z['t'];late=times>=times[-1]-.5
        normal_target=np.array(reports['stationary_5s']['manifest']['stationary_load_split_reference']['normal_target_BW'])
        basis=np.array([[target['activation'][k]-base['activation'][k] if name[:-2]==family else 0.
                         for k,name in enumerate(muscles)] for family in families])
        phase_basis.append(basis)
        jp=mujoco.MjData(model);jp.qpos[:]=target['qpos'];jp.act[:]=target['activation'];mujoco.mj_forward(model,jp);ar.update(jp)
        tangent=DynamicThelenBank(target['parameters'],float(model.opt.timestep))
        alo=np.maximum(.01,target['activation']-.0001);ahi=np.minimum(.98,target['activation']+.0001)
        tangent.initialize(jp.actuator_length,alo);flo=tangent.force.copy();tangent.initialize(jp.actuator_length,ahi)
        gain=(tangent.force-flo)/(ahi-alo);source_jacobians.append(ar.E.T@(-jp.actuator_moment.T*gain[None,:]))
        probe=mujoco.MjData(model);probe.qpos[:]=q[-1];probe.act[:]=z['activation'][-1]
        probe.xfrc_applied[model.body('pelvis').id,2]=.6*float(model.body_mass.sum()*9.81)
        mujoco.mj_inverse(model,probe);required=probe.qfrc_inverse-_external(model,probe)
        static=DynamicThelenBank(target['parameters'],float(model.opt.timestep));static.initialize(probe.actuator_length,probe.act)
        target_error=q[-1,ar.addresses]-target['qpos'][ar.addresses]
        summary=dict(accepted_load_bearing_motion=False,scope='Unmodified native-transition diagnostic,not an accepted swing teacher. One common C2 phase multiplies the complete equilibrium activation difference,with no allocator or independent neural gains.',
                     common_ramp_s=1.812,ramp_derivation='Six previously accepted302-ms load-transfer segments; no duration sweep. Native activation/fiber dynamics retained; no extra engineering command filter.',
                     source_transition_result=str(motion/'TRANSITION_RESULT.json'),physical_metrics=transition['metrics'],
                     normal_target_BW=normal_target.tolist(),final_normal_BW=z['normal_side_BW'][-1].tolist(),
                     final_normal_error_BW=(z['normal_side_BW'][-1]-normal_target).tolist(),
                     maximum_root_translation_from_DS_m=float(np.max(abs(q[:,:3]-q0[:3]))),
                     maximum_root_rotation_from_DS_deg=float(np.rad2deg(np.max(abs(q[:,3:6]-q0[3:6])))),
                     final_anatomical_error_deg=np.rad2deg(target_error).tolist(),anatomical_joints=[model.joint(int(j)).name for j in ar.jids],
                     final_sagittal_angles_deg=np.rad2deg(q[-1,addresses]).tolist(),target_sagittal_angles_deg=np.rad2deg(target['qpos'][addresses]).tolist(),
                     last_half_joint_speed_max_deg_s=float(np.rad2deg(np.max(abs(v[late][:,ar.rows])))),
                     final_dynamic_force_static_residual_max=float(np.max(abs(-probe.actuator_moment.T@z['source_force_N'][-1]-required))),
                     final_static_force_residual_max=float(np.max(abs(-probe.actuator_moment.T@static.force-required))),
                     interpretation='Endpoint source capacity is feasible,but interpolation does not reach the target equilibrium. Full contact/activation/fiber-consistent motion is required before neural training.')
        ref=np.load(spec['source_target']);rate=np.max(abs(ref['qvel'][:,addresses]),axis=0);acc=np.max(abs(ref['qacc'][:,addresses]),axis=0)
        tests=dict(root_translation_5mm=np.max(abs(q[:,:3]-q0[:3]),axis=1)>.005,
                   root_rotation_2deg=np.max(abs(q[:,3:6]-q0[3:6]),axis=1)>np.deg2rad(2.),
                   source_velocity=np.any(abs(v[:,addresses])>rate[None,:]+1e-10,axis=1),
                   source_acceleration=np.any(abs(z['qacc'][:,addresses])>acc[None,:]+1e-10,axis=1))
        summary['first_postprocessed_guard_times_s']={key:float(times[np.flatnonzero(mask)[0]]) if mask.any() else None for key,mask in tests.items()}
        (motion/'STATIONARY_MOTION_AUDIT.json').write_text(json.dumps(summary,indent=2));row['transport']=summary
        result['sides'][side]=row
    result['accepted_bilateral_stationary_references']=all(row['accepted_stationary_reference'] for row in result['sides'].values())
    result['accepted_dynamic_teacher']=False
    phase_basis=np.array(phase_basis);basis_error=float(np.max(abs(phase_basis.sum(axis=1)-np.array([np.load(folder/(s+'_stationary_5s')/'source_recruitment.npz')['activation']-base['activation'] for s in ('right','left')]))))
    np.savez_compressed(folder/'STATIONARY_STANCE_BASIS.npz',phases=np.array(['right_unload','left_unload']),families=np.array(families),muscles=base['muscles'],baseline_activation=base['activation'],phase_family_activation_delta=phase_basis,source_static_torque_activation_jacobian=np.array(source_jacobians),source_parameters=base['parameters'])
    result['proposed_shared_stance_basis']=dict(file='STATIONARY_STANCE_BASIS.npz',families=len(families),phase_order=['right_unload','left_unload'],sum_reconstruction_error=basis_error,
        definition='Each source-muscle family contributes its two homologous endpoint activation differences. A future coefficient must be common to both sides and both phases. All-one coefficients reproduce the failed direct interpolation exactly; no coefficients have been fitted.',
        source_jacobian='Equality-aware17-coordinate derivative of steady source torque versus activation at each accepted endpoint; centered1e-4 activation probe,one-sided at source activation bounds.',accepted_dynamic_primitive=False)
    (folder/'CHECKPOINT.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(accepted_bilateral_stationary_references=result['accepted_bilateral_stationary_references'],accepted_dynamic_teacher=False,
                         fresh_max_errors={s:max(r['fresh_numeric_errors'].values()) for s,r in result['sides'].items()}),indent=2))
    return result


def audit_predicted_branches(model_folder,start_folder,target_folder,out):
    """Read-only shadow of the source force predictor,with unchanged replay hash."""
    from predictive_source_recruitment import AugmentedTransition,NativeContactKinematics
    from source_force_comparison import load_recruitment
    import zlib
    if (out/'PREDICTOR_BRANCH_AUDIT.json').exists():raise ValueError('Preserve completed shadow audit')
    a=load_recruitment(start_folder);target=load_recruitment(target_folder);p=a['parameters'];names=a['muscles'].tolist()
    m=mujoco.MjModel.from_binary_path(str(model_folder/'runtime_model.mjb'));m.opt.timestep=.001
    r=AugmentedTransition(m,p,.05);m,d,b=r.model,r.data,r.backend;dt=float(m.opt.timestep);bw=float(m.body_mass.sum()*9.81)
    d.qpos[:]=a['qpos'];d.act[:]=a['activation'];d.ctrl[:]=a['activation'];d.xfrc_applied[m.body('pelvis').id,2]=.6*bw
    mujoco.mj_forward(m,d);b.bank.initialize(d.actuator_length,d.act);m.actuator_biasprm[:,0]=-b.bank.force;mujoco.mj_forward(m,d);r.commit_snapshot()
    shadow=DynamicThelenBank(p,dt);contact=NativeContactKinematics(m,d);mtp=np.array([int(m.joint('mtp_angle_'+s).qposadr[0]) for s in ('r','l')]);selected=[names.index(n) for n in ('glut_min1_r','med_gas_l','ext_dig_r','flex_dig_r','ext_hal_r','flex_hal_r')]
    n=round(7.312/dt);digests=[];streams=[];max_force_error=0.;started=perf_counter();rows=[];pred_codes=[];actual_codes=[];constraint_window=[]
    def fv_code(state,old,activation):
        length=np.hypot(state,p[:,3]);velocity=(state-old)*state/(length*dt*p[:,1]);c=.25+.75*np.maximum(.01,activation);k=(2.+2./p[:,6])/(p[:,7]-1.);edge=.95*p[:,7];threshold=c*(edge-1.)/(k*(p[:,7]-edge))
        return np.select([velocity/p[:,8]<=-c,velocity<=0.,velocity/p[:,8]<=threshold],[0,1,2],default=3).astype(np.uint8)
    for instrument in (False,True):
        r.reset();checksum=hashlib.sha256();phase=[]
        for i in range(n):
            fraction=np.clip((d.time-.5)/1.812,0.,1.);blend=fraction**3*(10.+fraction*(-15.+6.*fraction));command=a['activation']+blend*(target['activation']-a['activation']);r.excitation+=r.beta*(command-r.excitation);d.ctrl[:]=r.excitation
            if instrument:
                shadow.state[:]=b.bank.state;shadow.evaluate(d.actuator_length+dt*d.actuator_velocity,d.act);predicted=shadow.proposed.copy();pred_force=shadow.force.copy();old=shadow.state.copy();alpha=d.act.copy();pc=fv_code(predicted,old,alpha)
            b.step();measure=contact.measure(True)
            checksum.update(np.r_[d.qpos,d.qvel,d.act,b.bank.state,b.bank.force,d.ctrl].tobytes())
            phase.append(np.r_[d.time,d.qpos[mtp],d.qvel[mtp],d.qacc[mtp],measure['normal_by_body_N'],measure['signature_u32'],zlib.crc32(np.asarray(d.efc_state[:d.nefc],np.int32).tobytes())])
            if instrument:
                max_force_error=max(max_force_error,float(max(abs(pred_force-b.predicted_force))));ac=fv_code(b.bank.state,b.old_state,b.activation)
                if d.time>=6.812-1e-10:
                    rows.append(np.r_[d.time,b.predicted_force[selected],b.bank.force[selected],b.length_tangent[selected],pc[selected],ac[selected],phase[-1][-2:]])
                    pred_codes.append(pc);actual_codes.append(ac)
                    constraint_window.append(np.c_[d.efc_type[:d.nefc],d.efc_id[:d.nefc],d.efc_state[:d.nefc]].copy())
        digests.append(checksum.hexdigest());streams.append(np.array(phase))
    parity=float(np.max(abs(streams[0]-streams[1])));verified=bool(digests[0]==digests[1] and parity==0. and max_force_error<=1e-8)
    maxrows=max(len(row) for row in constraint_window);constraint_codes=np.full((len(constraint_window),maxrows,3),-1,np.int32)
    events=[]
    for k,row in enumerate(constraint_window):
        constraint_codes[k,:len(row)]=row
        if k and row.shape==constraint_window[k-1].shape and np.array_equal(row[:,:2],constraint_window[k-1][:,:2]):
            for j in np.flatnonzero(row[:,2]!=constraint_window[k-1][:,2]):
                kind,identity,state=map(int,row[j]);joint=m.joint(int(m.dof_jntid[identity])).name if kind==int(mujoco.mjtConstraint.mjCNSTR_FRICTION_DOF) else None
                events.append(dict(time_s=float(rows[k][0]),type=kind,identity=identity,joint=joint,previous_state=int(constraint_window[k-1][j,2]),state=state))
    out.mkdir(parents=True,exist_ok=True);np.savez_compressed(out/'shadow_branch_trace.npz',phase=streams[1],selected=np.array([names[j] for j in selected]),selected_trace=np.array(rows),predicted_FV_branch=np.array(pred_codes),corrected_FV_branch=np.array(actual_codes),constraint_type_id_state=constraint_codes,muscles=a['muscles'])
    pc=np.array(pred_codes);ac=np.array(actual_codes);mismatches=(pc!=ac).sum(axis=0)
    report=dict(observer_verified=verified,physical_trajectory_sha256=digests[0],shadow_trajectory_sha256=digests[1],phase_parity_max=parity,predicted_force_replay_error_N=max_force_error,source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),model=str(model_folder),start=str(start_folder),target=str(target_folder),timestep_s=dt,support_fraction=.6,native_steps=2*n,wall_s=perf_counter()-started,
                selected_trace_columns=['time_s']+[prefix+names[j] for prefix in ('predicted_force_N_','corrected_force_N_','predicted_length_tangent_N_m_','predicted_FV_branch_','corrected_FV_branch_') for j in selected]+['contact_topology_crc32','constraint_state_crc32'],
                late_predictor_corrector_FV_branch_mismatches=[dict(muscle=names[j],frames=int(mismatches[j])) for j in np.argsort(mismatches)[::-1] if mismatches[j]>0],late_constraint_state_events=events,
                scope='Independent shadow evaluates the unchanged constitutive predictor without feeding forces,tangents,or state back into physics. Off/on full-state+force stream hashes must match. Branch association does not establish numerical or physiological causality.')
    (out/'PREDICTOR_BRANCH_AUDIT.json').write_text(json.dumps(report,indent=2));return report


def plot_native_branches(folder):
    """Plot observed terminal toe/contact/FV histories without model changes."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'legend.fontsize':10,'svg.fonttype':'none'})
    fig,axes=plt.subplots(3,2,figsize=(10,7),sharex='col');colors=('#0077BB','#EE7733')
    for col,side in enumerate(('right','left')):
        z=np.load(folder/(side+'_baseline_native_branches.npz'));x=z['phase'];window=x[:,0]>=6.812;x=x[window];t=x[:,0];names=z['muscles'].tolist();codes=z['muscle_branch'][window]
        for j,s in enumerate(('r','l')):
            style='-' if j==0 else '--';marker='o' if j==0 else 's'
            axes[0,col].plot(t,np.rad2deg(x[:,3+j]),style,color=colors[j],label='Right MTP' if j==0 else 'Left MTP',marker=marker,markevery=60,markersize=3)
            axes[1,col].plot(t,x[:,8+2*j],style,color=colors[j],label='Right toes' if j==0 else 'Left toes',marker=marker,markevery=60,markersize=3)
        for j,name in enumerate(('glut_min1_r','med_gas_l')):
            axes[2,col].step(t,codes[:,names.index(name)]&3,where='post',color=colors[j],linestyle='-' if j==0 else '--',label=name)
        axes[0,col].set_title('Right-to-left load target' if col==0 else 'Left-to-right load target')
        axes[1,col].set_ylim(bottom=0.);axes[2,col].set_yticks([1,2]);axes[2,col].set_yticklabels(['Short./iso.','Lengthening']);axes[2,col].set_ylim(.8,2.2);axes[2,col].set_xlabel('Time [s]')
        for row in range(3):
            axes[row,col].grid(alpha=.25);axes[row,col].spines[['top','right']].set_visible(False)
        spike=t[abs(x[:,5])>5.]
        if len(spike):
            for row in range(3):axes[row,col].axvspan(float(spike[0]),float(t[-1]),color='#BBBBBB',alpha=.3)
    axes[0,0].set_ylabel('MTP velocity [deg/s]');axes[1,0].set_ylabel('Toe normal force [N]');axes[2,0].set_ylabel('Corrected FV branch')
    for row in range(3):axes[row,0].legend(loc='best',frameon=False)
    fig.suptitle('Terminal response of the 1% load-transfer reference: 1 ms physics',fontsize=12)
    fig.text(.5,.012,'Shading marks right-MTP acceleration spikes. Stationary motion criteria fail; toe contacts remain loaded.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.035,1,.965));fig.savefig(folder/'native_branch_response.png',dpi=220);fig.savefig(folder/'native_branch_response.svg');plt.close(fig)
    (folder/'native_branch_response.alt.txt').write_text('Three rows compare right-to-left and left-to-right small load targets at one-millisecond physics. MTP angular velocities and actual toe-body normal forces are plotted above corrected force-velocity branch codes for right gluteus minimus1 and left medial gastrocnemius. Toe forces stay positive. The left-to-right target develops right-MTP reversals near the end,coincident with source muscle branch transitions. Grey shading denotes acceleration spikes,not a causal attribution. Both mechanical trajectories fail stationarity.')


def summarize_native_branches(folder):
    """Compact phase/constitutive diagnostics from the bounded native audit."""
    history=json.loads((folder/'SEARCH_HISTORY.json').read_text());manifest=json.loads((folder/'MANIFEST.json').read_text())
    finish=folder/'MANIFOLD_SHOOTING.json';result=dict(status='complete' if finish.exists() else 'partial',evaluations=len(history),
        source=str(folder),task='Full shared-family derivative audit; no physical or neural changes.',
        row_coverage=dict(minimum=min(manifest['retained_squared_row_coverage']),maximum=max(manifest['retained_squared_row_coverage'])),baseline={},probes=[])
    for side in ('right','left'):
        z=np.load(folder/(side+'_baseline_native_branches.npz'));columns=z['phase_columns'].tolist();phase=z['phase'];late=phase[:,0]>=manifest['warm_s']+manifest['ramp_s']+manifest['hold_s']-.5
        x=phase[late];codes=z['muscle_branch'][late];names=z['muscles'].tolist();summary={}
        for channel in ('toes_r_N','toes_l_N','calcn_r_N','calcn_l_N','solver_iter','nefc'):
            v=x[:,columns.index(channel)];summary[channel]=dict(minimum=float(v.min()),maximum=float(v.max()),std=float(v.std()))
        summary['contact_topology_signatures']=int(len(np.unique(x[:,columns.index('contact_crc32')])));summary['wrapping_count_signatures']=int(len(np.unique(x[:,columns.index('wrap_crc32')])))
        summary['MTP']={}
        dt=float(np.median(np.diff(x[:,0])))
        for s in ('r','l'):
            velocity=x[:,columns.index('vmtp_'+s)];acc=x[:,columns.index('amtp_'+s)];time=x[:,0]-x[0,0]
            detrended=velocity-np.c_[np.ones(len(time)),time]@np.linalg.lstsq(np.c_[np.ones(len(time)),time],velocity,rcond=None)[0]
            power=abs(np.fft.rfft(detrended*np.hanning(len(time))))**2;freq=np.fft.rfftfreq(len(time),dt);peak=int(np.argmax(power[1:])+1)
            summary['MTP'][s]=dict(velocity_range_rad_s=[float(velocity.min()),float(velocity.max())],acceleration_range_rad_s2=[float(acc.min()),float(acc.max())],velocity_sign_changes=int(np.sum(velocity[1:]*velocity[:-1]<0.)),detrended_velocity_peak_bin_Hz=float(freq[peak]))
        summary['corrected_muscle_branch_transitions']={}
        for name,mask in (('force_velocity',3),('passive_fiber',4),('tendon',24),('activation_rising',32)):
            changed=(codes[1:]&mask)!=(codes[:-1]&mask);counts=changed.sum(axis=0)
            summary['corrected_muscle_branch_transitions'][name]=[dict(muscle=names[j],count=int(counts[j]),first_time_s=float(x[1:][np.flatnonzero(changed[:,j])[0],0])) for j in np.argsort(counts)[::-1] if counts[j]>0][:15]
        result['baseline'][side]=summary
    for row in history:
        b=row['branches'];phases=[]
        for v in b:
            comparison=dict(v.get('native_branch_comparison') or {})
            if comparison.get('late_changed_frames',{}).get('muscle_branch')==0:comparison['most_changed_muscles']=[]
            phases.append(dict(guard=v['guard'],native_branch_comparison=comparison,late_joint_speed_deg_s=v['late_joint_speed_deg_s'],final_static_residual_max=v['final_static_residual_max']))
        result['probes'].append(dict(evaluation=row['evaluation'],nonzero_mode_coordinates=np.flatnonzero(np.asarray(row['z'])!=0.).tolist(),cost=row['cost'],phases=phases))
    derivative=folder/'shooting_derivatives.npz'
    if derivative.exists():
        z=np.load(derivative);blocks=[('mean_pose',23),('mean_velocity',23),('mean_fiber',92),('mean_dynamic105',105),('mean_static105',105),('mean_normal',2),('std_velocity',23),('std_dynamic105',105),('std_static105',105),('guard',1)]
        result['physical_residual_J_singular_values_without_parameter_regularizer']=np.linalg.svd(z['J'][:-manifest['shared_modes']],compute_uv=False).tolist()
        width=sum(size for _,size in blocks);edges=np.r_[0,np.cumsum([size for _,size in blocks])];rows=[]
        scale=np.load(folder/'shooting_scale_checks.npz') if (folder/'shooting_scale_checks.npz').exists() else None
        for side,offset in zip(('right','left'),(0,width)):
            for (name,_),a,b in zip(blocks,edges[:-1],edges[1:]):
                measured=z['directional_measured'][offset+a:offset+b];predicted=z['directional_predicted'][offset+a:offset+b];signal=float(np.linalg.norm(measured))
                row=dict(side=side,block=name,measured_norm=signal,absolute_prediction_error=float(np.linalg.norm(predicted-measured)),relative_error=float(np.linalg.norm(predicted-measured)/signal) if signal>1e-8 else None)
                if scale is not None and signal>1e-8:
                    row.update(half_scale_relative_error=float(np.linalg.norm(2*scale['half'][offset+a:offset+b]-measured)/signal),double_scale_relative_error=float(np.linalg.norm(.5*scale['double'][offset+a:offset+b]-measured)/signal),opposite_scale_relative_error=float(np.linalg.norm(-scale['negative'][offset+a:offset+b]-measured)/signal))
                rows.append(row)
        result['directional_blocks']=rows
    shadow_path=folder/'predictor_shadow_v1/PREDICTOR_BRANCH_AUDIT.json'
    if shadow_path.exists():
        report=json.loads(shadow_path.read_text());shadow=np.load(shadow_path.parent/'shadow_branch_trace.npz');baseline=np.load(folder/'left_baseline_native_branches.npz')['phase'];phase=shadow['phase']
        report['fresh_baseline_phase_error_max']=float(np.max(abs(phase[:,:11]-baseline[:,:11])));report['fresh_baseline_contact_signature_error_max']=float(np.max(abs(phase[:,11]-baseline[:,15])))
        report['whole_trajectory_constraint_state_signatures']=int(len(np.unique(phase[:,12])));trace=shadow['selected_trace'];columns=report['selected_trace_columns'];j=columns.index('predicted_FV_branch_glut_min1_r');k=columns.index('corrected_FV_branch_glut_min1_r');mismatch=trace[:,j]!=trace[:,k]
        tangent=trace[:,columns.index('predicted_length_tangent_N_m_glut_min1_r')]
        report['glut_min1_r_predicted_tangent_range_N_m']=[float(tangent.min()),float(tangent.max())]
        report['glut_min1_r_branch_disagreement_times_s']=trace[mismatch,0].tolist();report['glut_min1_r_adapter_velocity_coefficient_magnitude_range_Ns_m']=(np.array(report['glut_min1_r_predicted_tangent_range_N_m'])*report['timestep_s']).tolist()
        result['predictor_shadow']=report
    result['scope']='Contact signature tracks active contact geometry/count,not every friction active-set face. Wrapping signature tracks cardinality,not wrap-side topology. Muscle codes describe the corrected constitutive solution,not the separate predicted-force branch. Frequency is a detrended finite-window bin,not an identified natural mode.'
    (folder/'NATIVE_BRANCH_AUDIT.json').write_text(json.dumps(result,indent=2));return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--model',required=True,type=Path)
    parser.add_argument('--endpoint',type=Path)
    parser.add_argument('--spec',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--allow-unaccepted-static',action='store_true')
    parser.add_argument('--physics-dt',type=float)
    parser.add_argument('--initial-audit-only',action='store_true')
    parser.add_argument('--campaign-summary',type=Path)
    parser.add_argument('--native-branch-audit',type=Path)
    parser.add_argument('--predictor-branch-target',type=Path)
    parser.add_argument('--predictor-branch-start',type=Path)
    parser.add_argument('--plot-native-branches',type=Path)
    args=parser.parse_args()
    if args.plot_native_branches:
        plot_native_branches(args.plot_native_branches)
    elif args.predictor_branch_target:
        if args.predictor_branch_start is None:parser.error('Predictor branch audit requires --predictor-branch-start')
        print(json.dumps(audit_predicted_branches(args.model,args.predictor_branch_start,args.predictor_branch_target,args.out),indent=2))
    elif args.native_branch_audit:
        print(json.dumps(summarize_native_branches(args.native_branch_audit),indent=2))
    elif args.campaign_summary:
        summarize_stationary_campaign(args.model,args.campaign_summary,args.spec)
    elif args.endpoint:
        validate(args.model,args.endpoint,args.spec,args.out,args.allow_unaccepted_static,args.physics_dt,args.initial_audit_only)
    else:parser.error('Provide --endpoint or --campaign-summary')
