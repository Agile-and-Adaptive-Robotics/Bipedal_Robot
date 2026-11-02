"""Isolated matched Gait standing probes; existing production files are untouched."""
import os,sys,json,time,hashlib,argparse,copy
from pathlib import Path
HERE=Path(__file__).resolve().parent;SPINAL=HERE.parents[1]/'spinal';REPO=HERE.parents[3]
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
sys.path[:0]=[str(SPINAL),str(SPINAL/'w2l_cpg'),str(SPINAL/'hunt_training'),str(SPINAL/'functional_gait2392')]
import numpy as np,mujoco
from scipy.spatial import ConvexHull
import runner as R,build_network as B,params as P
from sparse_pairwise_eval import SparsePairwiseForward

def mirror(n):
 for a,b in (('_r','_l'),('_l','_r'),('_R','_L'),('_L','_R')):
  if n.endswith(a):return n[:-2]+b
 return n

def canonical_neuron_names(net):
 names=[None]*len(net.compiled.V)
 for name,index in net.idx.items():
  if names[index] is None:names[index]=name
 if any(n is None for n in names):raise ValueError('Missing canonical neuron name')
 return names

def tie_weights(net):
 names=canonical_neuron_names(net);lookup={n:i for i,n in enumerate(names)};c=net.compiled;changed=[];seen=set()
 for post,pre in zip(*np.nonzero(c.g_max_non)):
  a,b=names[pre],names[post];ma,mb=mirror(a),mirror(b)
  if (ma,mb)==(a,b) or ma not in lookup or mb not in lookup:continue
  mp,mq=lookup[mb],lookup[ma];pair=tuple(sorted(((int(post),int(pre)),(int(mp),int(mq)))))
  if pair in seen:continue
  seen.add(pair);g,h=c.g_max_non[post,pre],c.g_max_non[mp,mq]
  if abs(g-h)<=1e-14:continue
  if g and h and c.del_e[post,pre]!=c.del_e[mp,mq]:raise ValueError('Cannot tie mirrored synapses of opposite type')
  if not g:
   for key in ('e_lo','e_hi','del_e'):getattr(c,key)[post,pre]=getattr(c,key)[mp,mq]
  if not h:
   for key in ('e_lo','e_hi','del_e'):getattr(c,key)[mp,mq]=getattr(c,key)[post,pre]
  value=float((g+h)/2);c.g_max_non[post,pre]=value;c.g_max_non[mp,mq]=value
  changed.append(dict(source=a,target=b,source_g=float(g),mirror_g=float(h),pooled_g=value))
 return changed

def run(family,out,duration=8.,tag='baseline',gain=1.,backend='sparse',muscle_afferent_scale=1.,standing_bias_scale=1.,ib_group_gain=1.,bpa_configuration=None,sensory_replay=None,mechanical_push_bw=0.,compiled_transform=None,fast_motor_io=True,pinned_actuator_settling=False,settling_reference_configuration=None,response_onset_physics_s=0.,rig_rotational_scale=1.,rig_lateral_scale=1.,rig_vertical_scale=1.,release_pressure_reset=False,release_history_reset=False,pathway_group_scales=None,capture_neural_inputs=False):
 out.mkdir(parents=True,exist_ok=True);prefix=out/(family+'_'+tag)
 if prefix.with_suffix('.result.json').exists():raise FileExistsError('Preserve completed trial: '+str(prefix))
 os.environ['AARL_NET']=family;os.environ['AARL_NPZ']=str(prefix.with_suffix('.npz'))
 checkpoint=SPINAL/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/aligned_support_ladder_v1/CHECKPOINT.json';cp=json.loads(checkpoint.read_text());saved=np.load(Path(cp['configuration']).with_suffix('.npz'));q0=saved['qfull'][0].copy()
 forcepath=SPINAL/'physics_fidelity/shared_force_law_functional_v3.json';force=json.loads(forcepath.read_text());capture={};rows=[];neurons=[];pressure_rows=[];neural_inputs=[];extended_rows=[];motor_rows=[];stop=None
 replay=None
 if sensory_replay:
  with np.load(sensory_replay) as archive:replay={'inputs':archive['inputs'],'names':archive['names']}
 bpa_config=json.loads(Path(bpa_configuration).read_text()) if bpa_configuration else None
 if bpa_config is not None and not bpa_config.get('valid_complete_assignment'):raise ValueError('Incomplete BPA assignment; no silent biological/BPA mixture')
 reference_config=json.loads(Path(settling_reference_configuration).read_text()) if settling_reference_configuration else None
 if reference_config is not None and (not pinned_actuator_settling or bpa_config is None):raise ValueError('Reference pressure preparation requires opt-in pinned BPA settling')
 if (release_pressure_reset or release_history_reset) and (not pinned_actuator_settling or bpa_config is None):raise ValueError('Actuator initialization controls require pinned BPA preparation')
 if release_history_reset and not release_pressure_reset:raise ValueError('This bounded cold-start control resets pressure and history together')
 if min(rig_rotational_scale,rig_lateral_scale,rig_vertical_scale)<=0:raise ValueError('Rig scales must remain positive in this bounded experiment')
 warmup_steps=round(R.WARMUP/P.DT)
 oldpose=R.apply_start_pose;oldh=R.apply_harness;oldbuild=B.build;oldstep=mujoco.mj_step;oldsolver=R.solve_standing_activations
 def patch(m):
  for row in force['muscles'].values():
   for name in row['actuator_names']:
    i=m.actuator(name).id
    for key in ('gainprm','biasprm','dynprm','lengthrange'):getattr(m,'actuator_'+key)[i]=row[key]
  for side in ('r','l'):
   j=m.joint('knee_angle_'+side);m.jnt_limited[j.id]=True;m.jnt_range[j.id]=np.radians([-90.,0.]);m.jnt_margin[j.id]=np.radians(2.);m.jnt_solref[j.id]=[.004,1.];m.jnt_solimp[j.id]=[.999,.9999,.001,.5,2.]
 def pose(m,d):
  oldpose(m,d)
  if m.nq!=len(q0):raise ValueError('Corrected pose/model nq mismatch')
  d.qpos[:]=q0;patch(m);mujoco.mj_forward(m,d)
 R.apply_start_pose=pose
 def harness(*args,**kwargs):
  m=oldh(*args,**kwargs);patch(m);d=mujoco.MjData(m);mujoco.mj_forward(m,d)
  for axis in ('tx','ty','tz'):
   j=m.joint('pelvis_'+axis)
   if abs(d.xaxis[j.id,2])>.8:m.jnt_stiffness[j.id]*=.01*rig_vertical_scale;m.dof_damping[j.dofadr[0]]*=.1*4.2*np.sqrt(rig_vertical_scale)
   elif abs(d.xaxis[j.id,1])>.8:m.jnt_stiffness[j.id]*=rig_lateral_scale;m.dof_damping[j.dofadr[0]]*=np.sqrt(rig_lateral_scale)
  if rig_rotational_scale!=1.:
   for name in ('pelvis_tilt','pelvis_list','pelvis_rotation','lumbar_extension','lumbar_bending','lumbar_rotation','hip_rotation_r','hip_rotation_l'):
    j=m.joint(name);m.jnt_stiffness[j.id]*=rig_rotational_scale;m.dof_damping[j.dofadr[0]]*=np.sqrt(rig_rotational_scale)
  capture['model']=m;return m
 R.apply_harness=harness
 def solver(m,d,Fmax,preload=True):
  from physics_fidelity.standing_solver import solve_standing_activations
  activation,diag=solve_standing_activations(m,d,Fmax,preload=False,bilateral=True,return_diagnostics=True);capture['standing_solver']=diag;capture['standing_activation_targets']=activation.copy();return activation
 R.solve_standing_activations=solver
 def pressure_adapter(m,d):
  if 'pressure_adapter' not in capture:
   from native_bpa_force import NativeBPAForce
   cfg=copy.deepcopy(reference_config if reference_config is not None else bpa_config)
   if reference_config is not None:
    cfg['pressure_history_capacity_s']=max(cfg.get('pressure_history_capacity_s',0.),max(r['delay_s'] for r in cfg['actuators'].values()),max(r['delay_s'] for r in bpa_config['actuators'].values()))
   initial=np.zeros(m.nu)
   for aid in range(m.nu):
    if m.actuator_actadr[aid]>=0:initial[aid]=d.act[m.actuator_actadr[aid]]
   adapter=NativeBPAForce(m,cfg,float(m.opt.timestep),cfg.get('disabled_actuators',[]));capture['pressure_adapter']=adapter
   mujoco.mj_forward(m,d);adapter.precharge(d,initial);capture['initial_pressure_kpa']=adapter.bank.pressure.copy()
  return capture['pressure_adapter']
 def release_snapshot(n,label='release_state'):
  m,d=capture['bound_model'],capture['bound_data'];state=dict(qpos=d.qpos.copy(),qvel=d.qvel.copy(),act=d.act.copy(),ctrl=d.ctrl.copy(),qacc=d.qacc.copy(),qacc_warmstart=d.qacc_warmstart.copy(),physics_time=np.array(d.time))
  state.update(runner_neural_steps=np.array(capture['neural_steps']),current_neural_input=capture['current_input'].copy())
  for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last'):
   if hasattr(n.compiled,key):state['neural_'+key]=np.asarray(getattr(n.compiled,key)).copy()
  if compiled_transform is not None:
   transformed=n.compiled.export_state()
   if not isinstance(transformed,dict) or not transformed:raise ValueError('Transformed backend must export a nonempty numeric/string state dictionary')
   for key,value in transformed.items():
    array=np.asarray(value)
    if array.dtype.hasobject:raise ValueError('Object-valued transformed state cannot be audited: '+key)
    state['transformed_'+key]=array.copy()
  if bpa_config is not None:
   adapter=pressure_adapter(m,d);bank=adapter.bank
   state.update(pressure=bank.pressure.copy(),pressure_command=bank.command.copy(),pressure_history=bank.history.copy(),pressure_step_count=np.array(bank.step_count),direction_memory=adapter.direction_state.copy(),tension=bank.force.copy(),series_spring_energy_J=bank.force**2/(2*bank.parameters[:,4]))
  np.savez_compressed(prefix.with_suffix('.'+label+'.npz'),**state)
  if label=='release_state':capture['release_state_saved']=True
 def build(acts,dt,interleg):
  n=oldbuild(acts,dt=dt,interleg=interleg);capture['projection']=tie_weights(n)
  # One shared multiplier, never independent left/right fitting or MN cancellation.
  names=canonical_neuron_names(n)
  for post,pre in zip(*np.nonzero(n.compiled.g_max_non)):
   if names[post].startswith('MN_') and names[pre].startswith('PF_'):n.compiled.g_max_non[post,pre]*=gain
   if names[post].startswith('MN_') and names[pre].startswith('IBEXC'):n.compiled.g_max_non[post,pre]*=ib_group_gain
  if pathway_group_scales:
   from muscle_map import classify
   counts={key:0 for key in pathway_group_scales}
   for post,pre in zip(*np.nonzero(n.compiled.g_max_non)):
    if not names[post].startswith('MN_'):continue
    group=classify(names[post][3:]).groups[0];source=names[pre].split('_')[0]
    for key,value in pathway_group_scales.items():
     cls,target=key.split(':')
     if source==cls and target in ('ALL',group):
      if not 0<=value<=2:raise ValueError('Bounded pathway scales must be in[0,2]')
      n.compiled.g_max_non[post,pre]*=value;counts[key]+=1
   if any(value==0 for value in counts.values()):raise ValueError('Unused pathway/group scale: '+str(counts))
   capture['pathway_group_counts']=counts
  reference=copy.deepcopy(n.compiled);fast=SparsePairwiseForward(copy.deepcopy(reference));u=n.make_inputs()
  if 'DRIVE' in n.inputs:u[n.inputs.index('DRIVE')]=2.5
  err=0.
  for _ in range(64):reference.forward(u);fast.forward(u);err=max(err,float(np.max(abs(reference.V-fast.V))))
  if err>1e-10:raise ValueError('Projected-network backend parity failed')
  capture['parity_max_mV']=err;capture['network']=n;capture['neural_steps']=0;n.compiled=SparsePairwiseForward(n.compiled) if backend=='sparse' else n.compiled
  if compiled_transform is not None:
   if backend!='sparse':raise ValueError('Compiled transform requires the verified sparse baseline')
   n.compiled=compiled_transform(n.compiled,names)
   if pinned_actuator_settling and not all(callable(getattr(n.compiled,k,None)) for k in ('export_state','import_state')):raise ValueError('Pinned transformed networks require explicit complete export_state/import_state methods')
  n.fast_io=bool(fast_motor_io)
  if pinned_actuator_settling:
   previous_bind=getattr(n,'bind_mujoco',None)
   def bind(m,d):
    if previous_bind is not None:previous_bind(m,d)
    capture['bound_model']=m;capture['bound_data']=d;capture['model']=m
   n.bind_mujoco=bind
  native=n.step;current=np.empty_like(n.make_inputs());afferent_ids=np.array([i for i,p in enumerate(n.inputs) if p.startswith(('Ia_','II_','Ib_'))],dtype=int);bias_ids=np.array([i for i,p in enumerate(n.inputs) if p.startswith('POST_')],dtype=int)
  warm_controls=np.empty(len(acts));mn_indices=np.array([n.idx[n.mn_names[a]] for a in acts])
  feedback_ids=np.array([i for i,p in enumerate(n.inputs) if not(p.startswith(('DRIVE','POST_')) or p=='POSTURE')],dtype=int)
  if replay is not None and replay['names'].tolist()!=list(n.inputs):raise ValueError('Sensory replay input identity/order mismatch')
  def advance(x):
   index=capture['neural_steps'];capture['neural_steps']+=1;np.copyto(current,x);current[afferent_ids]*=muscle_afferent_scale;current[bias_ids]*=standing_bias_scale;capture['current_input']=current
   if pinned_actuator_settling and index==warmup_steps:
    if release_pressure_reset:
     release_snapshot(n,'prepared_state');m,d=capture['bound_model'],capture['bound_data'];adapter=pressure_adapter(m,d);bank=adapter.bank;bank.pressure[:]=0.;bank.force[:]=0.;adapter.nominal_force[:]=0.
     if release_history_reset:bank.history[:]=0.;bank.command[:]=0.
     mujoco.mj_forward(m,d)
    release_snapshot(n)
   if replay is not None and (not pinned_actuator_settling or index>=warmup_steps):
    if index>=len(replay['inputs']):raise ValueError('Sensory replay exhausted; no silent padding')
    current[feedback_ids]=replay['inputs'][index,feedback_ids]
   if bpa_config is not None or capture_neural_inputs:neural_inputs.append(current.copy())
   v=native(current)
   if pinned_actuator_settling and index<warmup_steps:
    if 'bound_data' not in capture:raise ValueError('Opt-in settling requires runner bind_mujoco hook')
    m,d=capture['bound_model'],capture['bound_data'];np.take(v,mn_indices,out=warm_controls);np.divide(warm_controls,R.E_HI,out=warm_controls);np.clip(warm_controls,0.,1.,out=warm_controls);d.ctrl[:]=warm_controls
    if bpa_config is not None:
     adapter=pressure_adapter(m,d);mujoco.mj_forward(m,d);adapter.advance(d,warm_controls)
    else:
     if m.na!=m.nu or np.any(m.actuator_dyntype!=mujoco.mjtDyn.mjDYN_MUSCLE):raise ValueError('Activation settling supports native scalar muscle activation states only')
     mujoco.mj_forward(m,d);d.act[:]+=P.DT*d.act_dot;np.clip(d.act,0.,1.,out=d.act)
   return v
  n.step=advance;return n
 B.build=build
 def step(m,d,*args,**kwargs):
  capture['runtime_data']=d
  if bpa_config is not None and capture.get('neural_steps'):
   adapter=pressure_adapter(m,d)
   if reference_config is not None and not capture.get('response_switched') and d.time>=response_onset_physics_s:
    adapter.configure_response(bpa_config);capture['response_switched']=True;capture['response_switch_physics_time_s']=float(d.time)
   # Pressure advances exactly once, outside the force callback. Refresh path
   # geometry explicitly; callbacks may be queried repeatedly without a tick.
   mujoco.mj_forward(m,d);adapter.advance(d,d.ctrl)
  pulse=0.;pulse_body=None
  if mechanical_push_bw and 3.<=d.time<3.15:
   pulse=float(mechanical_push_bw)*float(np.sum(m.body_mass)*9.81);pulse_body=m.body('torso').id;d.xfrc_applied[pulse_body,0]+=pulse
  applied=d.xfrc_applied.copy()
  try:oldstep(m,d,*args,**kwargs)
  finally:
   if pulse_body is not None:d.xfrc_applied[pulse_body,0]-=pulse
  if not capture.get('neural_steps'):return
  if capture_neural_inputs:motor_rows.append((float(d.time),d.act.copy(),d.ctrl.copy(),d.actuator_force.copy()))
  knees=np.array([d.qpos[m.joint('knee_angle_'+s).qposadr[0]] for s in ('r','l')]);sample=np.r_[d.time,knees]
  if not np.isfinite(d.qpos).all() or np.any(knees>1e-8):raise RuntimeError('Nonfinite/hyperextension safety stop')
  feet=np.zeros((2,3));buffer=np.zeros(6)
  for c in range(d.ncon):
   co=d.contact[c];ids=[m.body(int(m.geom_bodyid[g])).name for g in (co.geom1,co.geom2)]
   for side,suffix in enumerate(('_r','_l')):
    signs=[i for i,n in enumerate(ids) if n in ('calcn'+suffix,'toes'+suffix,'talus'+suffix)]
    if len(signs)!=1:continue
    mujoco.mj_contactForce(m,d,c,buffer);feet[side]+=(1 if signs[0]==1 else -1)*(co.frame.reshape(3,3).T@buffer[:3])
  rig=np.zeros(3)
  for axis in ('tx','ty','tz'):
   j=m.joint('pelvis_'+axis);rig+=d.xaxis[j.id]*(d.qfrc_passive[j.dofadr[0]]+d.qfrc_applied[j.dofadr[0]])
  if 'rig_joint_names' not in capture:
   capture['rig_joint_names']=['pelvis_tx','pelvis_ty','pelvis_tz','pelvis_tilt','pelvis_list','pelvis_rotation','lumbar_extension','lumbar_bending','lumbar_rotation','hip_rotation_r','hip_rotation_l']
   capture['rig_joint_parameters']={name:dict(stiffness=float(m.jnt_stiffness[m.joint(name).id]),damping=float(m.dof_damping[m.joint(name).dofadr[0]]),springref=float(m.qpos_spring[m.joint(name).qposadr[0]])) for name in capture['rig_joint_names']}
  generalized=np.array([d.qfrc_passive[m.joint(name).dofadr[0]]+d.qfrc_applied[m.joint(name).dofadr[0]] for name in capture['rig_joint_names']])
  axes=np.column_stack([d.xaxis[m.joint(name).id] for name in ('pelvis_tilt','pelvis_list','pelvis_rotation')]);root_torque=np.linalg.solve(axes.T,generalized[3:6])
  reference=d.xipos[m.body('pelvis').id];external_force=applied[:,:3].sum(0);external_torque=(applied[:,3:]+np.cross(d.xipos-reference,applied[:,:3])).sum(0)
  extended_rows.append(np.r_[d.time,rig,root_torque,generalized,external_force,external_torque])
  if 'foot_meshes' not in capture:
   geometry=[]
   for side in ('r','l'):
    shapes=[]
    for g in range(m.ngeom):
     if m.body(int(m.geom_bodyid[g])).name not in ('calcn_'+side,'toes_'+side) or m.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
     mesh=m.geom_dataid[g];a=m.mesh_vertadr[mesh];v=m.mesh_vert[a:a+m.mesh_vertnum[mesh]].copy()
     if len(v)>4:v=v[ConvexHull(v).vertices]
     shapes.append((g,v))
    if not shapes:raise ValueError('Actual foot meshes missing')
    geometry.append(shapes)
   capture['foot_meshes']=geometry
  clearance=[];plane=m.geom('ground-plane').id;normal=d.geom_xmat[plane].reshape(3,3)[:,2].copy()
  if normal[2]<0:normal=-normal
  for shapes in capture['foot_meshes']:clearance.append(min(float(np.min(v@(d.geom_xmat[g].reshape(3,3).T@normal))+(d.geom_xpos[g]-d.geom_xpos[plane])@normal) for g,v in shapes))
  rows.append(np.r_[sample,feet.ravel(),rig,clearance]);neurons.append(capture['network'].compiled.V.copy())
  if 'pressure_adapter' in capture:
   adapter=capture['pressure_adapter'];bank=adapter.bank
   pressure_rows.append(np.stack((bank.command,bank.pressure,bank.force,bank.flags,d.actuator_force[adapter.ids],adapter.nominal_force,adapter.effective_gain,adapter.direction_state)))
 mujoco.mj_step=step
 started=time.perf_counter();metrics={}
 try:
  if P.G.get('pm_gain',0.)!=0:raise ValueError('Disable external phase-machine motor overrides for neural-family comparisons')
  metrics=R.main(['--stand-eval',str(duration),'--reflex-stand'])
 except Exception as e:stop=str(e)
 finally:
  if 'pressure_adapter' in capture:capture['pressure_adapter'].close()
  mujoco.mj_step=oldstep;R.apply_start_pose=oldpose;R.apply_harness=oldh;B.build=oldbuild;R.solve_standing_activations=oldsolver
 samples=np.asarray(rows);m=capture.get('model');bw=float(np.sum(m.body_mass)*9.81) if m is not None else None
 if stop and 'runtime_data' in capture:
  failed=capture['runtime_data'];np.savez_compressed(prefix.with_suffix('.failure_state.npz'),qpos=failed.qpos,qvel=failed.qvel,act=failed.act,physics_time=failed.time,actuator_length=failed.actuator_length,actuator_force_last_forward=failed.actuator_force,neuro=capture['network'].compiled.V)
 summary=dict(machine='easteregg2',family=family,condition='Biological Hill muscles, corrected-pose initialization, shared static fit and pooled mirrored neural weights',duration_requested_s=duration,completed_samples=len(rows),failure=stop,backend=backend,wall_s=time.perf_counter()-started,shared_pf_mn_gain=gain,muscle_afferent_scale=muscle_afferent_scale,standing_bias_scale=standing_bias_scale,ib_group_gain=ib_group_gain,
  neurons=len(capture['network'].compiled.V) if 'network'in capture else None,projected_weight_pairs=len(capture.get('projection',[])),backend_parity_max_error_mV=capture.get('parity_max_mV'),metrics=metrics,
  pathway_group_scales=pathway_group_scales,pathway_group_scaled_edge_counts=capture.get('pathway_group_counts',{}),capture_neural_inputs=bool(capture_neural_inputs),compiled_transform=getattr(compiled_transform,'__qualname__',None),backend_parity_scope='Pre-transform graded reference vs sparse backend; no transformed/spiking equivalence claim' if compiled_transform else 'Native graded reference vs accelerated sparse backend',
  fast_motor_io=bool(fast_motor_io),
  pinned_actuator_settling=bool(pinned_actuator_settling),settling_reference_configuration=str(settling_reference_configuration) if settling_reference_configuration else None,response_onset_physics_s=response_onset_physics_s,response_switch_physics_time_s=capture.get('response_switch_physics_time_s'),release_state_saved=capture.get('release_state_saved',False),
  release_state_file=str(prefix.with_suffix('.release_state.npz')) if capture.get('release_state_saved') else None,
  release_pressure_reset=bool(release_pressure_reset),release_history_reset=bool(release_history_reset),
  rig_scales=dict(rotational=rig_rotational_scale,lateral=rig_lateral_scale,vertical=rig_vertical_scale),
  body_weight_N=bw,mean_ground_load_bw=float(np.mean(samples[:,[5,8]].sum(axis=1))/bw) if len(rows) else None,
  mean_upward_rig_support_bw=float(np.mean(np.maximum(samples[:,11],0))/bw) if len(rows) else None,
  mean_downward_rig_force_bw=float(np.mean(np.maximum(-samples[:,11],0))/bw) if len(rows) else None,mean_net_rig_force_bw=float(np.mean(samples[:,11])/bw) if len(rows) else None,
  actual_visual_foot_mesh_clearance_range_m=np.stack((samples[:,12:14].min(0),samples[:,12:14].max(0))).tolist() if len(rows) else None,
  standing_activation_targets=capture.get('standing_activation_targets',np.array([])).tolist(),standing_solver_diagnostics=capture.get('standing_solver'),
  max_knee_angle_deg=np.degrees(samples[:,1:3].max(0)).tolist() if len(rows) else None,
  mechanical_push=dict(body='torso',world_axis='+x' if mechanical_push_bw>=0 else '-x',force_bw=abs(mechanical_push_bw),start_sim_time_s=3.,duration_s=.15),
  apparatus='Existing pelvis springs; vertical stiffness0.01, damping multiplier4.2. Not harness-free.',frame='Legacy converted Gait z-up; distinct from y-up CAD.',
  time_conventions='Native NPZ time includes0.6s pinned initialization; apparatus/pressure data.time begins with free physics integration, excluding that initialization.',
  checkpoint=str(checkpoint),checkpoint_sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest(),force_law=str(forcepath),force_law_sha256=hashlib.sha256(forcepath.read_bytes()).hexdigest(),source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),SPINAL/'runner.py',SPINAL/'params.py',SPINAL/'build_network.py',SPINAL/'build_network_w2lvar.py',SPINAL/'build_network_syn6.py')},
  scope='Untrained family standing diagnostic, not the Hunt winner, gait/curriculum pass, or BPA case. Actual visual mesh clearance remains required before any ground-step promotion.')
 if bpa_config is not None:
  summary['condition']='Virtual two20mm-per-route BPA pressure-state force replacement, common biological pinned initializer, same neural controller; not mounted CAD hardware identification'
  summary['scope']='Declared virtual BPA standing comparison; static surface characterized, dimensions/pressure dynamics hypothetical. No calibrated CAD BPA or gait/curriculum claim.'
  summary['bpa_configuration']=str(Path(bpa_configuration).resolve());summary['bpa_configuration_sha256']=hashlib.sha256(Path(bpa_configuration).read_bytes()).hexdigest()
  summary['sensory_replay']=str(sensory_replay) if sensory_replay else None
  summary['sensory_replay_scope']='Every state-dependent neural input replayed from nominal; DRIVE/POSTURE/POST_ bias commands remain live. Local static pressure inversion still uses actual length in both conditions.' if sensory_replay else 'Live neural state feedback'
  summary['pressure_model_source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'native_bpa_force.py',HERE/'cad_interface.py',HERE/'actuator_interface.py')}
  if pressure_rows:
   data=np.asarray(pressure_rows);summary['pressure_peak_kpa']=float(np.max(data[:,1,:]));summary['inverse_saturation_fraction']=float(np.mean(data[:,3,:]!=0));summary['native_pressure_force_max_error_N']=float(np.max(abs(data[:,4,:]+data[:,2,:])))
   summary['initial_pressure_kpa']=capture['initial_pressure_kpa'].tolist()
   np.savez_compressed(prefix.with_suffix('.pressure.npz'),time=samples[:,0],samples=data,names=np.array(capture['pressure_adapter'].names),channels=['reference_pressure_kpa','actual_pressure_kpa','tension_N','inverse_flag','native_signed_force_N','nominal_tension_N','effective_force_gain','direction_memory'])
  np.savez_compressed(prefix.with_suffix('.inputs.npz'),inputs=np.asarray(neural_inputs),names=np.array(capture['network'].inputs))
 if capture_neural_inputs and bpa_config is None:np.savez_compressed(prefix.with_suffix('.inputs.npz'),inputs=np.asarray(neural_inputs),names=np.array(capture['network'].inputs))
 if motor_rows:np.savez_compressed(prefix.with_suffix('.motor_state.npz'),time=np.array([a[0] for a in motor_rows]),activation=np.array([a[1] for a in motor_rows]),excitation=np.array([a[2] for a in motor_rows]),engine_force_cache_N=np.array([a[3] for a in motor_rows]),actuator_names=np.array([m.actuator(i).name for i in range(m.nu)]),activation_addresses=m.actuator_actadr)
 prefix.with_suffix('.result.json').write_text(json.dumps(summary,indent=2,default=float));np.savez_compressed(prefix.with_suffix('.apparatus.npz'),samples=samples,neuro=np.asarray(neurons),columns=['time','knee_r','knee_l','right_GRFx','right_GRFy','right_GRFz','left_GRFx','left_GRFy','left_GRFz','rigFx','rigFy','rigFz','right_visual_mesh_clearance','left_visual_mesh_clearance'])
 np.savez_compressed(prefix.with_suffix('.apparatus_extended.npz'),samples=np.asarray(extended_rows),columns=['time','rig_world_Fx','rig_world_Fy','rig_world_Fz','pelvis_rig_world_Tx','pelvis_rig_world_Ty','pelvis_rig_world_Tz']+['generalized_passive_'+n for n in capture.get('rig_joint_names',[])]+['applied_world_Fx','applied_world_Fy','applied_world_Fz','applied_world_Tx_about_pelvis','applied_world_Ty_about_pelvis','applied_world_Tz_about_pelvis'],joint_parameters_json=json.dumps(capture.get('rig_joint_parameters',{})))
 (out/(family+'_'+tag+'_projection.json')).write_text(json.dumps(capture.get('projection',[]),indent=2));print(json.dumps(summary,indent=2,default=float));return summary

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--family',required=True,choices=['s3k','w2lvar','syn6']);p.add_argument('--out',required=True);p.add_argument('--duration',type=float,default=8.);p.add_argument('--tag',default='baseline');p.add_argument('--gain',type=float,default=1.);p.add_argument('--backend',choices=['sparse','reference'],default='sparse');p.add_argument('--muscle-afferent-scale',type=float,default=1.);p.add_argument('--standing-bias-scale',type=float,default=1.);p.add_argument('--ib-group-gain',type=float,default=1.);p.add_argument('--parameters-json');p.add_argument('--bpa-configuration')
 p.add_argument('--pinned-actuator-settling',action='store_true');p.add_argument('--settling-reference-configuration');p.add_argument('--response-onset-physics-s',type=float,default=0.);p.add_argument('--sensory-replay');p.add_argument('--release-pressure-reset',action='store_true');p.add_argument('--release-history-reset',action='store_true');p.add_argument('--rig-rotational-scale',type=float,default=1.);p.add_argument('--rig-lateral-scale',type=float,default=1.);p.add_argument('--rig-vertical-scale',type=float,default=1.);p.add_argument('--scalar-motor-io',action='store_true');a=p.parse_args()
 if a.parameters_json:
  if any(k in sys.argv for k in ('--gain','--ib-group-gain','--standing-bias-scale')):p.error('Use JSON parameters or scalar CLI overrides, not both')
  data=json.loads(Path(a.parameters_json).read_text());parameters=data.get('trial',data);parameters=parameters.get('params',parameters.get('parameters',parameters))
  if set(parameters)!={'pf_mn_gain','ib_group_gain','standing_bias_scale'}:p.error('JSON must specify exactly all three consumed shared standing parameters')
  a.gain=float(parameters['pf_mn_gain']);a.ib_group_gain=float(parameters['ib_group_gain']);a.standing_bias_scale=float(parameters['standing_bias_scale'])
 run(a.family,Path(a.out),a.duration,a.tag,a.gain,a.backend,a.muscle_afferent_scale,a.standing_bias_scale,a.ib_group_gain,a.bpa_configuration,sensory_replay=a.sensory_replay,fast_motor_io=not a.scalar_motor_io,pinned_actuator_settling=a.pinned_actuator_settling,settling_reference_configuration=a.settling_reference_configuration,response_onset_physics_s=a.response_onset_physics_s,rig_rotational_scale=a.rig_rotational_scale,rig_lateral_scale=a.rig_lateral_scale,rig_vertical_scale=a.rig_vertical_scale,release_pressure_reset=a.release_pressure_reset,release_history_reset=a.release_history_reset)
