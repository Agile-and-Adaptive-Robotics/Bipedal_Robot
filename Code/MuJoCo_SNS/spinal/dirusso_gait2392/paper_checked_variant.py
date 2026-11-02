"""Opt-in sensor/event comparison; preserves legacy core, W, and muscle physics.

This is a protected reconstruction, not an exact Di Russo reproduction. The
electrical probe runs without MuJoCo integration. A pilot must be requested
explicitly and always applies the existing Hunt knee bumper to BOTH cases.
"""
from pathlib import Path
import os, sys, json, hashlib, time, argparse
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):
    os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import numpy as np
import controller as C
from controller_numba import CompiledDiRusso

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load_weights(path=None):
    path=Path(path).resolve() if path is not None else HERE/'dirusso_W_stage1.npy'
    W=np.load(path)
    if W.shape!=(18,5) or not np.isfinite(W).all():raise ValueError('Expected finite18x5 weights')
    if path!=HERE/'dirusso_W_stage1.npy':
        if abs(W).max()>1.+1e-12 or not np.array_equal(W[:9],W[9:]):raise ValueError('New fit must be bounded[-1,1] and exactly mirrored')
    return path,W

def phase_initialization(policy):
    if policy=='source_default':return np.array([0.,.5]),dict(policy=policy)
    from bounded_primitive_fit import DATA,GRF,measured_onsets
    reference_t=float(np.load(DATA,allow_pickle=True)['t'][0]);onset=measured_onsets(GRF)['r']
    period=onset[1]-onset[0];right=((reference_t-onset[0])/period)%1.;left=(right+.5)%1.
    return np.array([left,right]),dict(policy=policy,reference_time_s=reference_t,right_heel_onset_s=onset[0],right_cycle_period_s=period,backward_periodic_phase_extrapolation=reference_t<onset[0],independent_left_phase_fitted=False,dataset_sha256=sha(DATA),GRF_sha256=sha(GRF))

def weights_probe(path,out,phase_init='source_default'):
    out.mkdir(parents=True,exist_ok=True);path,W=load_weights(path)
    a=CompiledDiRusso(W);b=CompiledDiRusso(np.zeros_like(W));assert np.array_equal(a.W,W)
    initial_phase,phase_audit=phase_initialization(phase_init);a.phi[:]=initial_phase;b.phi[:]=initial_phase
    sens=np.zeros((2,4,9));event=np.zeros(2,bool);maximum=0.
    first_phase=None
    for i in range(500):
        maximum=max(maximum,float(abs(a.step_array(sens,event,0.)-b.step_array(sens,event,0.)).max()))
        if i==0:first_phase=a.phi.copy()
    expected=(initial_phase+.002*5.11/(2*np.pi))%1.
    assert np.allclose(first_phase,expected,rtol=0.,atol=1e-12)
    assert maximum>0.
    result=dict(weights_path=str(path),weights_sha256=sha(path),compiled_W_array_equal=True,maximum_actual_MN_difference_vs_zero_weights=maximum,actual_steps=500,body_run=False,source_sha256=sha(__file__),initial_phase_cycles_LR=initial_phase.tolist(),actual_phase_after_first_step=first_phase.tolist(),expected_phase_after_first_step=expected.tolist(),phase_initialization=phase_audit)
    (out/'WEIGHTS_CONSUMPTION.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result
def receptors(length,velocity,force,variant):
    return C.DiRussoController.receptors(*( (length,velocity,force) if variant=='input_event_v1' else (velocity,length,force)))
def reset_signal(force,previous,last_strike,t,variant):
    edge=previous<=20.<force
    latest=t if edge else last_strike
    reset=edge if variant=='input_event_v1' else force>20. and (latest is None or t-latest>.45)
    return edge, bool(reset), latest

def electrical_probe(out):
    out.mkdir(parents=True,exist_ok=True)
    W=np.random.default_rng(42).uniform(-.2,.4,(18,5))
    gains={k:(v if np.isscalar(v) else v*.2) for k,v in C.default_gains().items()}
    a=C.DiRussoController(W,.002,omega=2*np.pi,gamma=1.,gains=gains)
    b=CompiledDiRusso(W,omega=2*np.pi,gamma=1.,gains=gains)
    sens=np.zeros((2,4,9)); strikes=np.zeros(2,dtype=np.bool_); error=0.
    for i in range(2000):
        t=i*.002; sens[:,0]=.15*(1.+np.sin(t*6.28));sens[:,1]=.5;sens[:,2]=.1;sens[:,3]=.4
        strikes[:]=False
        if i==300: strikes[0]=True
        sa={s:tuple(sens[j,:3])+(sens[j,3,0],) for j,s in enumerate('LR')}
        aa=a.step(sa,dict(zip('LR',strikes)),.01*np.sin(t)).copy()
        bb=b.step_array(sens,strikes,.01*np.sin(t)).copy()
        error=max(error,float(np.max(abs(aa-bb))))
        assert np.allclose(a.phi,b.phi,rtol=0,atol=1e-13)
    assert error<1e-12
    rec={}
    for variant in ('legacy_inputs','input_event_v1'):
        actual=np.array(receptors(np.full(9,.8),np.full(9,-.5),np.full(9,.3),variant))
        prev=20.1;last=None;events=[];resets=[]
        for i in range(1000):
            t=i*.002; f=100. if .1<=t<.9 or 1.1<=t<1.9 else 0.
            edge,reset,last=reset_signal(f,prev,last,t,variant);prev=f
            if edge:events.append(t)
            if reset:resets.append(t)
        rec[variant]=dict(receptors_first_channel=actual[:,0].tolist(),true_events=events,reset_count=len(resets),first_resets=resets[:4])
    assert rec['input_event_v1']['reset_count']==2
    assert np.array_equal(np.array(rec['input_event_v1']['receptors_first_channel']),[0.,.8,.3])
    assert rec['legacy_inputs']['reset_count']>2
    f6=np.array([100.,-4.,-2.,0.,0.,0.]);assert max(f6[0],0.)==100. and max(f6[2],0.)==0.
    # Both implementations are already warm; this is ONLY a neural-kernel benchmark.
    start=time.perf_counter()
    for _ in range(2000):a.step(sa,dict(L=False,R=False),0.)
    python_s=time.perf_counter()-start;start=time.perf_counter()
    for _ in range(2000):b.step_array(sens,strikes,0.)
    numba_s=time.perf_counter()-start
    result=dict(format='dirusso_input_event_electrical_v1',machine=os.environ.get('COMPUTERNAME'),core_max_difference=error,synthetic=rec,contact_force_fixture=dict(contact_frame_force_N=f6.tolist(),normal_N=100.,legacy_tangent2_N=0.),neural_only_warm_benchmark=dict(steps=2000,python_s=python_s,numba_s=numba_s,speedup=python_s/numba_s),source_sha256={p.name:sha(p) for p in (HERE/'controller.py',HERE/'controller_numba.py',HERE/'dirusso_fast.py',Path(__file__))},body_run=False)
    (out/'ELECTRICAL_PROOF.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2));return result

def initialize_dependencies(model,data,project=False):
    """Exact source joint-equality polynomials; no anatomical driver edits."""
    import mujoco
    links=[]
    for eq in range(model.neq):
        if model.eq_type[eq]==mujoco.mjtEq.mjEQ_JOINT and model.eq_obj2id[eq]>=0:
            links.append((int(model.eq_obj1id[eq]),int(model.eq_obj2id[eq]),model.eq_data[eq,:5].copy()))
    dependent={a for a,b,c in links};ordered=[];pending=list(links)
    while pending:
        unresolved={a for a,b,c in pending};ready=[r for r in pending if r[1] not in unresolved]
        if not ready:raise ValueError('Cyclic native joint-equality dependency')
        for r in ready:ordered.append(r);pending.remove(r)
    drivers=[j for j in range(model.njnt) if j not in dependent]
    qa=model.jnt_qposadr[drivers];va=model.jnt_dofadr[drivers]
    before_q=data.qpos[qa].copy();before_v=data.qvel[va].copy()
    def residuals():
        position=[];velocity=[]
        for dep,ind,c in ordered:
            a,b=int(model.jnt_qposadr[dep]),int(model.jnt_qposadr[ind]);av,bv=int(model.jnt_dofadr[dep]),int(model.jnt_dofadr[ind]);x=data.qpos[b]-model.qpos0[b]
            deriv=c[1]+x*(2*c[2]+x*(3*c[3]+4*x*c[4]))
            position.append(float(data.qpos[a]-model.qpos0[a]-np.polynomial.polynomial.polyval(x,c)))
            velocity.append(float(data.qvel[av]-deriv*data.qvel[bv]))
        return max(map(abs,position),default=0.),max(map(abs,velocity),default=0.)
    initial=residuals()
    if project:
        for dep,ind,c in ordered:
            a,b=int(model.jnt_qposadr[dep]),int(model.jnt_qposadr[ind]);av,bv=int(model.jnt_dofadr[dep]),int(model.jnt_dofadr[ind]);x=data.qpos[b]-model.qpos0[b]
            data.qpos[a]=model.qpos0[a]+np.polynomial.polynomial.polyval(x,c)
            data.qvel[av]=(c[1]+x*(2*c[2]+x*(3*c[3]+4*x*c[4])))*data.qvel[bv]
    final=residuals()
    if project:mujoco.mj_forward(model,data)
    eqmask=np.asarray(data.efc_type)==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY)
    audit=dict(projected=project,dependent_joint_count=len(ordered),initial_polynomial_position_max_m_or_rad=initial[0],initial_polynomial_velocity_max_mps_or_radps=initial[1],final_polynomial_position_max_m_or_rad=final[0],final_polynomial_velocity_max_mps_or_radps=final[1],native_equality_position_max=float(np.max(abs(data.efc_pos[eqmask]))) if eqmask.any() else 0.,driver_qpos_max_change=float(np.max(abs(data.qpos[qa]-before_q))),driver_qvel_max_change=float(np.max(abs(data.qvel[va]-before_v))))
    assert audit['driver_qpos_max_change']==audit['driver_qvel_max_change']==0.
    if project:assert final[0]<1e-12 and final[1]<1e-12
    return audit

def foot_plane_gaps(model,data):
    plane=model.geom('ground-plane').id;normal=data.geom_xmat[plane].reshape(3,3)[:,2];point=data.geom_xpos[plane];rows=[]
    for side in ('l','r'):
        for region in ('calcn','toes'):
            gid=model.geom(region+'_'+side+'_geom_1').id;mesh=int(model.geom_dataid[gid]);start=int(model.mesh_vertadr[mesh]);number=int(model.mesh_vertnum[mesh]);local=model.mesh_vert[start:start+number]
            world=local@data.geom_xmat[gid].reshape(3,3).T+data.geom_xpos[gid]
            rows.append(dict(side=side,region=region,gap_m=float(((world-point)@normal).min())))
    return rows,normal

def align_left_toe(model,data):
    import mujoco
    before,normal=foot_plane_gaps(model,data);gap=next(r['gap_m'] for r in before if r['side']=='l' and r['region']=='toes');j=model.joint('pelvis_ty');qa=int(j.qposadr[0]);axis=data.xaxis[j.id].copy();qbefore=data.qpos.copy();delta=(-.001-gap)/float(axis@normal)
    data.qpos[qa]+=delta;mujoco.mj_forward(model,data);after,_=foot_plane_gaps(model,data)
    target=next(r['gap_m'] for r in after if r['side']=='l' and r['region']=='toes');assert abs(target+.001)<1e-10
    changed=np.flatnonzero(data.qpos!=qbefore);assert np.array_equal(changed,[qa])
    toe=model.geom('toes_l_geom_1').id;plane=model.geom('ground-plane').id
    contacts=[dict(geom1=model.geom(c.geom1).name,geom2=model.geom(c.geom2).name,distance_m=float(c.dist)) for c in data.contact if {int(c.geom1),int(c.geom2)}=={toe,plane}]
    assert contacts,'Aligned toe must have an actual plane contact before stepping'
    return dict(policy='left_toe_m1mm',root_joint='pelvis_ty',root_coordinate_shift_m=float(delta),root_world_displacement_m=(axis*delta).tolist(),before_mesh_gaps=before,after_mesh_gaps=after,initial_toe_contacts=contacts,only_root_position_changed=True)

def anatomical_guards(model,policy):
    if policy=='none':return [],dict(policy=policy)
    import mujoco,xml.etree.ElementTree as ET
    path=HERE.parents[3]/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
    source={c.get('name'):np.fromstring(c.findtext('range',''),sep=' ') for c in ET.parse(path).getroot().iter('Coordinate')}
    rows=[];missing=[]
    for name in source:
        jid=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,name)
        if jid<0 or model.jnt_type[jid]!=mujoco.mjtJoint.mjJNT_HINGE or name.startswith('knee_angle_'):continue
        bounds=source[name]
        if len(bounds)!=2 or not np.isfinite(bounds).all() or bounds[1]<=bounds[0]:missing.append(name);continue
        rows.append(dict(name=name,address=int(model.jnt_qposadr[jid]),bounds_rad=bounds.tolist(),bounds_deg=np.degrees(bounds).tolist(),engine_limited=bool(model.jnt_limited[jid])))
    return rows,dict(policy=policy,source=str(path),source_sha256=sha(path),numeric_tolerance_deg=.05,coordinates=rows,unresolved=missing,note='These are the actual input-model bounds, including broad hips±120/ankles±90deg; they are rejection limits, not a normative healthy-walking envelope. No physics stop or state clamp is added.')

def check_anatomical_guard(data,rows):
    for row in rows:
        q=float(data.qpos[row['address']]);lo,hi=row['bounds_rad'];tol=np.radians(.05)
        if q<lo-tol or q>hi+tol:return dict(joint=row['name'],angle_deg=float(np.degrees(q)),bounds_deg=row['bounds_deg'],time_s=float(data.time),tolerance_deg=.05)
    return None

def environment(initializer='legacy_reference',ground_init='source_height',support=.7):
    import dirusso_fast as F
    e=F.get_env();m=e.model
    for side in ('r','l'):
        j=m.joint('knee_angle_'+side)
        m.jnt_limited[j.id]=True;m.jnt_range[j.id]=np.radians([-120.,0.])
        m.jnt_margin[j.id]=np.radians(2.);m.jnt_solref[j.id]=[.004,1.]
        m.jnt_solimp[j.id]=[.999,.9999,.001,.5,2.]
    e.set_reference_ic()
    source_velocity=e.data.qvel.copy()
    if initializer=='equality_static':e.data.qvel[:]=0.
    e.initialization_audit=initialize_dependencies(m,e.data,initializer in ('equality_consistent','equality_static'))
    e.initialization_audit['explicit_zero_velocity_condition']=initializer=='equality_static'
    e.initialization_audit['source_reference_knee_velocity_deg_s']=[float(np.degrees(source_velocity[m.joint('knee_angle_'+s).dofadr[0]])) for s in ('l','r')]
    e.initialization_audit['selected_knee_velocity_deg_s']=[float(np.degrees(e.data.qvel[m.joint('knee_angle_'+s).dofadr[0]])) for s in ('l','r')]
    e.ground_audit=dict(policy=ground_init)
    if ground_init=='left_toe_m1mm':
        # The same prescribed support must already be present when evaluating
        # initial contact forces, not introduced only after the first sensor read.
        e.data.xfrc_applied[:]=0.;e.data.xfrc_applied[e.pelvis,2]=support*e.weight
        e.ground_audit=align_left_toe(m,e.data)
    return e

def pilot(out,variant,duration=4.,support=.7,preflight=False,initializer='legacy_reference',weights_path=None,phase_init='source_default',ground_init='source_height',rom_guard='none'):
    import mujoco
    out.mkdir(parents=True,exist_ok=True)
    if (out/'RESULT.json').exists():raise FileExistsError('Use a new output directory to preserve the previous experiment.')
    started=time.perf_counter();e=environment(initializer,ground_init,support);m,d=e.model,e.data
    guard_rows,guard_audit=anatomical_guards(m,rom_guard)
    knees=np.array([int(m.joint('knee_angle_'+s).qposadr[0]) for s in ('l','r')])
    names=[joint+'_'+s for s in ('l','r') for joint in ('hip_flexion','knee_angle','ankle_angle')]
    addresses=np.array([int(m.joint(n).qposadr[0]) for n in names])
    Wpath,W=load_weights(weights_path);initial_phase,phase_audit=phase_initialization(phase_init)
    cfg=dict(format='dirusso_protected_input_event_v1',variant=variant,machine=os.environ.get('COMPUTERNAME'),duration_s=duration,support_fraction=support,omega_rad_s=5.11,gamma_rad_s=1.,dt_s=.002,knee_limit_deg=[-120.,0.],knee_tolerance_rad=1e-4,initial_knees_deg=np.degrees(d.qpos[knees]).tolist(),W_path=str(Wpath),W_sha256=sha(Wpath),W_min=float(W.min()),W_max=float(W.max()),W_LR_max_difference=float(abs(W[:9]-W[9:]).max()),source_sha256={p.name:sha(p) for p in (HERE/'controller.py',HERE/'controller_numba.py',HERE/'dirusso_fast.py',HERE/'plant.py',Path(__file__))},scope='Protected sagittal Gait2392 reconstruction; unchanged 34 grouped muscle actuators, wrapped primitives, W_stage1, gains, balance, and source force normalization. No training or paper-objective claim.',knee_patch_source=str(HERE.parent/'hunt_training/physical_eval.py')+':53-60',neural_output_units='dimensionless activation, not mV/nA',joint_names=names,control_side_order=['L','R'],physical_frame=dict(gravity=m.opt.gravity.tolist(),up=(-m.opt.gravity/np.linalg.norm(m.opt.gravity)).tolist()))
    cfg['compiled_knees']=[dict(name='knee_angle_'+s,limited=bool(m.jnt_limited[m.joint('knee_angle_'+s).id]),range_deg=np.degrees(m.jnt_range[m.joint('knee_angle_'+s).id]).tolist()) for s in ('l','r')]
    cfg['initializer']=initializer;cfg['initialization_audit']=e.initialization_audit
    cfg['phase_initialization']=phase_audit;cfg['initial_phase_cycles_LR']=initial_phase.tolist()
    cfg['ground_initialization']=e.ground_audit;cfg['anatomical_ROM_guard']=guard_audit
    cfg['weights_origin']='unchanged legacy W_stage1' if weights_path is None else 'explicit opt-in bounded mirrored weights, not a published-paper fit'
    if weights_path is not None:cfg['scope']='Protected sagittal Gait2392 reconstruction with explicitly loaded bounded mirrored weights. Original neural core, gains, support policy, and muscle/contact physics are retained. No physical success or paper-objective replication is implied.'
    cfg['trace_timing']='qpos/qvel/activation are integration state t; cached geometry/contact/force fields queried by the legacy loop are from the preceding force evaluation (approximately t-dt after the first sample). Preserved caller semantics, not a zero-delay sensor claim.'
    (out/'CONFIG.json').write_text(json.dumps(cfg,indent=2))
    # Save actual compiled model, including the identical knee protection.
    # FastEnv subsequently loads its full reference model. mj_saveLastXML would
    # therefore refer to a different last-loaded XML; the compiled binary is exact.
    mujoco.mj_saveModel(m,str(out/'runtime_model.mjb'),None)
    guard='initial knee hyperextension' if np.any(d.qpos[knees]>1e-4) else None
    guard_details=check_anatomical_guard(d,guard_rows)
    if guard_details:guard='initial anatomical ROM violation'
    if preflight or guard:
        result=dict(preflight_only=preflight,guard=guard,ready_for_pilot=guard is None,setup_wall_s=time.perf_counter()-started,**cfg)
        (out/'RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result
    ctrl=CompiledDiRusso(W,omega=5.11,gamma=1.)
    ctrl.phi[:]=initial_phase
    sens=np.zeros((2,4,9));previous=np.full(2,20.1);last=[None,None];strike=np.zeros(2,dtype=bool);f6=np.zeros(6)
    n=int(duration/.002);qlog=np.empty((n,m.nq));vlog=np.empty((n,m.nv));mlog=np.empty((n,2,9));plog=np.empty((n,2));prims=np.empty((n,2,5));snlog=np.empty((n,2,4,9));force=np.zeros((n,2,2));reset=np.empty((n,2),bool);ground=np.zeros((n,2,3));count=0
    controls=np.empty((n,m.nu));tensions=np.empty((n,m.nu));constraint=np.empty((n,m.nv));normal_alignment=1.
    activation=np.empty((n,m.na));lengths=np.empty((n,m.nu));velocities=np.empty((n,m.nu));group_torque=np.empty((n,2,9,6));ff_drive=np.empty((n,2,9));total_drive=np.empty((n,2,9));balance_output=np.empty((n,2,9));com_log=np.empty((n,3));regional=np.zeros((n,2,2,3))
    dof=np.array([int(m.joint(n).dofadr[0]) for n in names]);plane_ids=set(np.flatnonzero(m.geom_type==mujoco.mjtGeom.mjGEOM_PLANE).tolist());contact_pairs={}
    up=-m.opt.gravity/np.linalg.norm(m.opt.gravity);foot_bodies=[{m.body('calcn_'+s).id,m.body('toes_'+s).id} for s in ('l','r')]
    body_ids=m.geom_bodyid;setup_s=time.perf_counter()-started;start=time.perf_counter()
    # Compile the neural kernel before the measured integrated loop.
    warm=CompiledDiRusso(W);warm.step_array(sens,np.zeros(2,bool),0.)
    warmup_s=time.perf_counter()-start
    loop_start=time.perf_counter()
    for i in range(n):
        t=i*.002;signals=np.zeros((2,2))
        for k in range(d.ncon):
            co=d.contact[k];mujoco.mj_contactForce(m,d,k,f6)
            pair=tuple(sorted((int(co.geom1),int(co.geom2))))
            rec=contact_pairs.setdefault(pair,[0,0.,0.]);rec[0]+=1;rec[1]+=float(f6[0]);rec[2]=max(rec[2],float(f6[0]))
            cue=max(f6[0 if variant=='input_event_v1' else 2],0.)
            for si in range(2):
                if co.geom1==e.heel_gid[si] or co.geom2==e.heel_gid[si]:signals[si,0]+=cue
                elif co.geom1==e.toe_gid[si] or co.geom2==e.toe_gid[si]:signals[si,1]+=cue
                # Independent physical force, with correct contact-frame rotation.
                for own,other,sign in ((co.geom1,co.geom2,-1.),(co.geom2,co.geom1,1.)):
                    if int(body_ids[own]) in foot_bodies[si] and int(other) in plane_ids:
                        world_force=sign*np.asarray(co.frame).reshape(3,3).T@f6[:3]
                        ground[i,si]+=world_force
                        region=0 if int(body_ids[own])==m.body('calcn_'+('l','r')[si]).id else 1
                        regional[i,si,region]+=world_force
                        normal_alignment=min(normal_alignment,abs(float(np.dot(np.asarray(co.frame).reshape(3,3)[0],up))))
        for si in range(2):
            length=np.empty(9);velocity=np.empty(9);load=np.empty(9)
            for gsi,mi,idx in e.grp_slices:
                if gsi!=si:continue
                length[mi]=np.clip((np.mean(d.actuator_length[idx])-e.lo[si,mi])/e.span[si,mi],0.,2.)
                velocity[mi]=np.clip(np.mean(d.actuator_velocity[idx])/e.span[si,mi],-5.,5.)
                load[mi]=np.clip(np.mean(abs(d.actuator_force[idx]))/e.fmax[si,mi],0.,2.)
            sens[si,:3]=receptors(length,velocity,load,variant)
            sens[si,3,0]=np.clip(signals[si].sum()/700.,0.,1.)
            edge,strike[si],last[si]=reset_signal(signals[si,0],previous[si],last[si],t,variant);previous[si]=signals[si,0]
        old_mnstate=ctrl.mnstate.copy()
        mn=ctrl.step_array(sens,strike,float(d.qpos[e.tilt_adr]))
        d.ctrl[e.ctrl_flat]=mn.reshape(-1)[e.ctrl_mn_idx]
        d.xfrc_applied[:]=0.;d.xfrc_applied[e.pelvis,2]=support*e.weight
        qlog[i]=d.qpos;vlog[i]=d.qvel;mlog[i]=mn;plog[i]=ctrl.phi;snlog[i]=sens;force[i]=signals;reset[i]=strike
        controls[i]=d.ctrl;tensions[i]=d.actuator_force;constraint[i]=d.qfrc_constraint
        activation[i]=d.act;lengths[i]=d.actuator_length;velocities[i]=d.actuator_velocity;com_log[i]=d.subtree_com[0]
        total_drive[i]=old_mnstate+(ctrl.mnstate-old_mnstate)*C.TAU_A/.002
        balance_output[i]=mn-np.clip(ctrl.mnstate,0.,1.)
        for gsi,mi,idx in e.grp_slices:group_torque[i,gsi,mi]=d.actuator_force[idx]@d.actuator_moment[np.ix_(idx,dof)]
        prims[i]=np.stack([C.primitives(p) for p in ctrl.phi]);ff_drive[i]=(W@prims[i].T).reshape(2,9,2)[np.arange(2),:,np.arange(2)]*ctrl.ff;count=i+1
        mujoco.mj_step(m,d)
        if not np.isfinite(d.qpos).all():guard='nonfinite';break
        if np.any(d.qpos[knees]>1e-4):guard='knee hyperextension';break
        guard_details=check_anatomical_guard(d,guard_rows)
        if guard_details:guard='anatomical ROM violation';break
    elapsed=time.perf_counter()-loop_start
    np.savez_compressed(out/'TRACE.npz',t=np.arange(count)*.002,qpos=qlog[:count],qvel=vlog[:count],mn=mlog[:count],phase_cycles=plog[:count],primitives=prims[:count],receptor_inputs=snlog[:count],selected_contact_cues_N=force[:count],reset_event=reset[:count],foot_ground_force_world_N=ground[:count],support_force_world_N=np.tile([0.,0.,support*e.weight],(count,1)),joint_angles_rad=qlog[:count,addresses],joint_names=np.array(names),muscle_names=np.array(C.MUSCLES),ctrl=controls[:count],actuator_force_N=tensions[:count],qfrc_constraint=constraint[:count],activation=activation[:count],actuator_length_m=lengths[:count],actuator_velocity_mps=velocities[:count],group_joint_torque_Nm=group_torque[:count],feedforward_drive=ff_drive[:count],total_integrator_drive=total_drive[:count],balance_output=balance_output[:count],com_world_m=com_log[:count],calcn_toe_ground_force_world_N=regional[:count])
    result=dict(guard=guard,completed=guard is None and count==n,recorded_steps=count,simulated_s=float(d.time),integrated_wall_s=elapsed,setup_wall_s=setup_s,neural_warmup_wall_s=warmup_s,cold_total_wall_s=time.perf_counter()-started,real_time_factor=float(d.time)/elapsed,neural_warmup_excluded=True,knee_min_deg=np.degrees(np.minimum(qlog[:count,knees].min(axis=0),d.qpos[knees])).tolist(),knee_max_deg=np.degrees(np.maximum(qlog[:count,knees].max(axis=0),d.qpos[knees])).tolist(),reset_counts=reset[:count].sum(axis=0).tolist(),mean_foot_upward_N=(ground[:count]@up).mean(axis=0).tolist(),ground_contact_normal_min_abs_up_dot=normal_alignment,prescribed_upward_support_N=support*e.weight,body_weight_N=e.weight,config=str(out/'CONFIG.json'),trace=str(out/'TRACE.npz'),training_run=False)
    result['guard_details']=guard_details;result['post_step_terminal_joint_angles_deg']=np.degrees(d.qpos[addresses]).tolist()
    result['contact_pairs']=[dict(geom1=m.geom(pair[0]).name,geom2=m.geom(pair[1]).name,samples=val[0],mean_normal_N=val[1]/val[0],peak_normal_N=val[2]) for pair,val in contact_pairs.items()]
    (out/'RESULT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['electrical','preflight','pilot','weights_probe'],default='electrical');p.add_argument('--variant',choices=['legacy_inputs','input_event_v1','both'],default='input_event_v1');p.add_argument('--initializer',choices=['legacy_reference','equality_consistent','equality_static'],default='legacy_reference');p.add_argument('--weights',type=Path);p.add_argument('--phase-init',choices=['source_default','reference_antiphase'],default='source_default');p.add_argument('--ground-init',choices=['source_height','left_toe_m1mm'],default='source_height');p.add_argument('--rom-guard',choices=['none','source_opensim'],default='none');p.add_argument('--out',type=Path,required=True);p.add_argument('--duration',type=float,default=4.);p.add_argument('--support',type=float,default=.7);a=p.parse_args()
    if a.mode=='electrical':electrical_probe(a.out)
    elif a.mode=='weights_probe':weights_probe(a.weights,a.out,a.phase_init)
    elif a.variant=='both':
        for variant in ('legacy_inputs','input_event_v1'):pilot(a.out/variant,variant,a.duration,a.support,a.mode=='preflight',a.initializer,a.weights,a.phase_init,a.ground_init,a.rom_guard)
    else:pilot(a.out,a.variant,a.duration,a.support,a.mode=='preflight',a.initializer,a.weights,a.phase_init,a.ground_init,a.rom_guard)
