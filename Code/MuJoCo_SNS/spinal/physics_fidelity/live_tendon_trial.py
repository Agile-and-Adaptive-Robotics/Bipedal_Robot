"""Opt-in subprocess-only source muscle backend for existing live SNS protocols.

No production file is changed. Native, rigid, and compliant conditions share
the selected saved controller and its original posture/apparatus/initializer.
"""
import os,sys,json,time,runpy,argparse,hashlib
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent;S=HERE.parent;REPO=S.parents[2]
sys.path[:0]=[str(HERE),str(S),str(S/'functional_gait2392')]
import numpy as np,mujoco
from thelen_dynamic import parameters_from_osim,ConstitutiveBackend

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def assess_syn6(out,physics):
 bpa=REPO/'Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination';sys.path[:0]=[str(bpa),str(S/'functional_gait2392')];from gait_standing_objective import score
 prefix=out/('syn6_source_'+physics);r=json.loads(prefix.with_suffix('.result.json').read_text());_,a=score(r,prefix);(out/'OFFLINE_ASSESSMENT.json').write_text(json.dumps(a,indent=2));return a
def summarize(root):
 records={};numeric=('t','act','ctrl','q','qfull','com','neuro','contact','footz')
 def compare(a,b,keys):
  x,y=np.load(a,allow_pickle=True),np.load(b,allow_pickle=True);return {k:dict(bit_exact=bool(np.array_equal(x[k],y[k])),max_abs_difference=float(np.max(abs(x[k]-y[k])))) for k in keys if k in x and k in y}
 for mode in ('native','rigid','compliant'):
  folder=root/f'live_syn6_{mode}_v1';r=json.loads((folder/f'syn6_source_{mode}.result.json').read_text());a=json.loads((folder/'OFFLINE_ASSESSMENT.json').read_text());b=json.loads((folder/'BACKEND_RESULT.json').read_text());records[mode]=dict(accepted=a['accepted'],tilt_deg=a['tilt_deg'],sway_m=a['sway_m'],mean_ground_load_bw=r['mean_ground_load_bw'],upward_rig_bw=r['mean_upward_rig_support_bw'],downward_rig_bw=r['mean_downward_rig_force_bw'],knee_max_deg=r['max_knee_angle_deg'],protocol_wall_s=r['wall_s'],backend_diagnostics=b['backend_diagnostics'],source=str(folder))
 cp=json.loads(Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\STANDING_CHECKPOINT.json').read_text());source=Path(cp['fresh_result'].replace('.result.json','.npz'))
 parity=compare(source,root/'live_syn6_native_v1/syn6_source_native.npz',numeric)
 fresh=compare(root/'live_syn6_compliant_v1/backend_trace.npz',root/'live_syn6_compliant_fresh_v2/backend_trace.npz',('t','qpos','force','activation','ctrl'))
 answer=dict(family='Constrained syn6 standing candidate, not Hunt',standing=records,native_saved_candidate_numeric_parity=parity,compliant_fresh_numeric_parity=fresh,native_saved_candidate_bit_exact=all(v['bit_exact'] for v in parity.values()),compliant_fresh_bit_exact=all(v['bit_exact'] for v in fresh.values()),source_validation=str(root/'dynamic_source_v2/VALIDATION.json'),scope='Existing assisted standing protocol, same neural parameters and initializer. No gait or independent standing claim. Inner runner force-profile metadata describes initialization; runtime backend is explicitly identified here.')
 hunt={}
 for mode in ('native','rigid','compliant'):
  folder=root/f'live_hunt_{mode}_v1'
  if (folder/'controller.result.json').exists():
   r=json.loads((folder/'controller.result.json').read_text());hunt[mode]=dict(assessment=r['assessment'],foot_clearance=r.get('foot_clearance'),ground_contacts=r.get('ground_contacts'),source=str(folder))
 answer['hunt_ground_transfer']=hunt
 coarse=np.load(root/'source_dynamic_full_compliant_v1/compliant_damp1_beta0.npz');medium=np.load(root/'source_dynamic_compliant_1ms/compliant_damp1_beta0.npz');fine=np.load(root/'source_dynamic_compliant_0p5ms/compliant_damp1_beta0.npz')
 import runner as R
 from mechanics import JOINTS
 model=mujoco.MjModel.from_xml_path(str(R.MODEL));ix=[model.joint(n).qposadr[0] for n in JOINTS];convergence={}
 for label,z,stride in [('2ms',coarse,4),('1ms',medium,2)]:
  n=len(fine['t'])//stride;dq=z['qpos'][:n,ix]-fine['qpos'][stride-1::stride,ix];df=z['force'][:n]-fine['force'][stride-1::stride];convergence[label]=dict(joint_RMS_deg=float(np.degrees(np.sqrt(np.mean(dq*dq)))),force_RMS_N=float(np.sqrt(np.mean(df*df))))
 answer['source_dynamic_full_model_refinement_vs_0p5ms']=convergence
 (root/'SOURCE_DYNAMIC_CHECKPOINT.json').write_text(json.dumps(answer,indent=2));print(json.dumps({k:answer[k] for k in ('native_saved_candidate_bit_exact','compliant_fresh_bit_exact','source_dynamic_full_model_refinement_vs_0p5ms')},indent=2))
def main():
 parser=argparse.ArgumentParser();parser.add_argument('--protocol',choices=['hunt-ground','hunt-stand','hunt-air1','hunt-air2','syn6-stand'],required=True);parser.add_argument('--physics',choices=['native','rigid','compliant'],required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--beta',type=float,default=0.);parser.add_argument('--duration',type=float);parser.add_argument('--config',type=Path);parser.add_argument('--source-task-basis',action='store_true');parser.add_argument('--pulse-N',type=float,default=0.);parser.add_argument('--pulse-onset',type=float,default=4.);parser.add_argument('--pulse-duration',type=float,default=.15);parser.add_argument('--pulse-axis',default='+x');parser.add_argument('--assess-only',action='store_true');parser.add_argument('--summarize',action='store_true');parser.add_argument('--export-model-only',action='store_true');args=parser.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
 if args.summarize:return summarize(out)
 if args.assess_only:print(json.dumps(assess_syn6(out,args.physics),indent=2));return
 if (out/'BACKEND_RESULT.json').exists():raise FileExistsError('Preserve completed trial')
 expected=8. if args.protocol in ('hunt-stand','syn6-stand') else 16.
 duration=args.duration if args.duration is not None else expected
 if args.protocol.startswith('hunt-') and duration!=expected:raise ValueError('Preserve the original Hunt stage duration and scoring window')
 if args.protocol=='hunt-stand' and args.config is None:raise ValueError('Standing requires an explicit same-controller stage3 config')
 from single_pulse import SinglePulse
 pulse=SinglePulse(args.pulse_N,args.pulse_onset,args.pulse_duration,args.pulse_axis)
 source=S/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/contact_feasible_start_v1/ground_guide.json'
 if args.protocol in ('hunt-air1','hunt-air2'):source=S/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/zero_hyperextension_air_v1/shared.json'
 if args.config is not None:source=args.config.resolve()
 cp=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\STANDING_CHECKPOINT.json')
 physical_source=REPO/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
 if args.source_task_basis:
  if args.physics=='native':raise ValueError('Source task-basis diagnostic requires a source force backend')
  sys.path.insert(0,str(S/'hunt_training'))
  from balance_diagnostics import install_source_task_basis
  install_source_task_basis(physical_source,args.physics)
 manifest=dict(protocol=args.protocol,physics=args.physics,beta=args.beta,duration=duration,controller_source=str(source if args.protocol.startswith('hunt-') else cp),controller_source_sha256=sha(source if args.protocol.startswith('hunt-') else cp),physical_source_sha256=sha(physical_source),module_sha256=sha(HERE/'thelen_dynamic.py'),driver_sha256=sha(__file__),runner_sha256=sha(S/'runner.py'),scope='Live SNS, explicit saved/experimental Hunt config. Native initializer and native sensory-normalization constants retained identically. Source dynamic laws have .01 minimum activation and shared bilateral physical parameters.')
 manifest['isolated_pulse']=dict(force_N=args.pulse_N,onset_s=args.pulse_onset,duration_s=args.pulse_duration,native_axis=args.pulse_axis,clock='actual MuJoCo data.time before step; excludes frozen neural warmup',source_sha256=sha(HERE/'single_pulse.py'))
 manifest['source_task_basis']=dict(enabled=args.source_task_basis,source_sha256=sha(HERE/'balance_diagnostics.py'))
 full_dynamics=bool(json.loads(source.read_text()).get('trace_full_dynamics',False)) if args.protocol.startswith('hunt-') else False;manifest['trace_full_dynamics']=full_dynamics
 (out/'MANIFEST.json').write_text(json.dumps(manifest,indent=2));original=mujoco.mj_step;models={};log={};guard=None;protocol_result=None;ib_session=None;aff_session=None
 def step(m,d,*a,**kw):
  if a or kw:raise RuntimeError('Stateful backend requires one ordinary step; no hidden multi-step/adaptive rollback')
  if args.export_model_only:
   from tendon_replay_campaign import R
   fields=('qpos0','qpos_spring','body_mass','body_inertia','body_ipos','body_iquat','body_parentid','dof_armature','dof_damping','dof_frictionloss','jnt_type','jnt_pos','jnt_axis','jnt_range','jnt_margin','jnt_solref','jnt_solimp','jnt_limited','jnt_stiffness','geom_type','geom_size','geom_pos','geom_quat','geom_contype','geom_conaffinity','geom_friction','geom_solref','geom_solimp','geom_margin','pair_geom1','pair_geom2','pair_friction','pair_solref','pair_solimp','pair_margin','eq_type','eq_obj1id','eq_obj2id','eq_data','eq_active','actuator_gainprm','actuator_biasprm','actuator_dynprm','actuator_lengthrange','mesh_vert')
   mujoco.mj_saveLastXML(str(out/'runtime_model.xml'),m);mujoco.mj_saveModel(m,str(out/'runtime_model.mjb'),None);np.savez_compressed(out/'runtime_model_arrays.npz',**{name:getattr(m,name).copy() for name in fields},initial_qpos=d.qpos.copy(),initial_qvel=d.qvel.copy(),initial_activation=d.act.copy())
   (out/'MODEL_EXPORT.json').write_text(json.dumps(dict(scope='Exact live-model export before the first advancing physical step; initializer and pose-held neural initialization retained. This is an audit,not a failed walking simulation.',timestep=m.opt.timestep,integrator=int(m.opt.integrator),gravity=m.opt.gravity.tolist(),solver=int(m.opt.solver),cone=int(m.opt.cone),iterations=m.opt.iterations,tolerance=m.opt.tolerance,physical_time_s=d.time,source_model=str(R.MODEL)),indent=2));raise RuntimeError('Export-only audit complete before physical advancement')
  key=(id(m),id(d))
  if key not in models:
   if m.nu!=92:raise RuntimeError('Expected actual92-muscle Gait2392 model')
   names=[m.actuator(i).name for i in range(m.nu)];p=parameters_from_osim(physical_source,names,args.beta);backend=None
   if args.physics!='native':backend=ConstitutiveBackend(m,d,p,args.physics)
   if aff_session is not None:aff_session.bind_backend(m,d,backend)
   models[key]=(backend,names,p);capacity=round(duration/m.opt.timestep)+2
   log[key]=dict(t=np.empty(capacity),qpos=np.empty((capacity,m.nq)),force=np.empty((capacity,m.nu)),activation=np.empty((capacity,m.nu)),ctrl=np.empty((capacity,m.nu)),fiber_state=np.zeros((capacity,m.nu)),fiber_force_components=np.zeros((capacity,3,m.nu)),applied_pulse=np.zeros((capacity,3)),external_pelvis_wrench=np.zeros((capacity,6)),n=0,knee_max=d.qpos[[m.joint('knee_angle_'+s).qposadr[0] for s in ('r','l')]].copy())
   if full_dynamics:
    for field in ('qvel','qacc','qfrc_bias','qfrc_passive','qfrc_actuator','qfrc_constraint','qfrc_smooth','qfrc_contact','qfrc_limit','qfrc_equality'):log[key][field]=np.zeros((capacity,m.nv))
  backend,names,p=models[key]
  failure=None
  z=log[key];i=z['n'];pb=m.body('pelvis').id
  try:
   with pulse.applied(d,pb,d.time) as delivered:
    z['applied_pulse'][i]=delivered;z['external_pelvis_wrench'][i]=d.xfrc_applied[pb]
    if backend is None:original(m,d)
    else:
     if abs(m.opt.timestep-backend.dt)>1e-12:raise RuntimeError('Timestep changed without fiber-state rollback')
     backend.step()
  except Exception as exc:failure=exc
  z=log[key];i=z['n'];z['t'][i]=d.time;z['qpos'][i]=d.qpos;z['force'][i]=-d.actuator_force;z['activation'][i]=d.act;z['ctrl'][i]=d.ctrl;z['n']=i+1
  if backend is not None:z['fiber_state'][i]=backend.bank.state;z['fiber_force_components'][i]=backend.components
  if full_dynamics:
   for field in ('qvel','qacc','qfrc_bias','qfrc_passive','qfrc_actuator','qfrc_constraint','qfrc_smooth'):z[field][i]=getattr(d,field)
   for field,types in [('qfrc_contact',(5,6,7)),('qfrc_limit',(3,4)),('qfrc_equality',(0,))]:mujoco.mj_mulJacTVec(m,d,z[field][i],np.where(np.isin(d.efc_type,types),d.efc_force,0.))
  z['knee_max']=np.maximum(z['knee_max'],d.qpos[[m.joint('knee_angle_'+s).qposadr[0] for s in ('r','l')]])
  if failure is not None:raise failure
  if np.any(z['knee_max']>0):raise RuntimeError('Strict zero-hyperextension gate')
 mujoco.mj_step=step;tic=time.perf_counter()
 try:
  if args.protocol.startswith('hunt-'):
   cfg=json.loads(source.read_text())
   if cfg.get('isolated_posture_scales'):
    sys.path.insert(0,str(S/'hunt_training'))
    from balance_diagnostics import install_posture_scales
    install_posture_scales(cfg['isolated_posture_scales'])
   if cfg.get('isolated_ib_load'):
    sys.path.insert(0,str(S/'hunt_training'))
    from stance_ib import install
    stage=int(args.protocol[-1]) if args.protocol in ('hunt-air1','hunt-air2') else int(cfg['stage'])
    ib_session=install(cfg['isolated_ib_load'],physical_source,stage,cfg.get('reflex_scale',1.))
   if cfg.get('isolated_afferents'):
    if cfg.get('backend','sparse')!='sparse':raise ValueError('Current-accounting audit requires sparse forward wrapper')
    sys.path.insert(0,str(S/'hunt_training'))
    from source_afferents import install
    stage=int(args.protocol[-1]) if args.protocol in ('hunt-air1','hunt-air2') else int(cfg['stage'])
    aff_session=install(cfg['isolated_afferents'],physical_source,stage)
   if args.protocol in ('hunt-air1','hunt-air2'):cfg['stage']=int(args.protocol[-1])
   if args.protocol=='hunt-stand' and cfg['stage']!=3:raise ValueError('Standing config must explicitly select stage3')
   if args.physics!='native':cfg['diagnose_swing']=False # native zero-activation probe is invalid for stateful affine forces
   target=out/'controller.json';target.write_text(json.dumps(cfg,indent=2));sys.argv=[str(S/'hunt_training/physical_eval.py'),str(target)];runpy.run_path(sys.argv[0],run_name='__main__')
   protocol_result=json.loads((out/'controller.result.json').read_text())
  else:
   bpa=REPO/'Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination';sys.path.insert(0,str(bpa));import gait_standing_family as G
   checkpoint=json.loads(cp.read_text());p=checkpoint['parameters'];protocol_result=G.run('syn6',out,duration,tag='source_'+args.physics,gain=p['pf_mn_gain'],standing_bias_scale=p['standing_bias_scale'],ib_group_gain=p['ib_group_gain']);protocol_result['assessment']=assess_syn6(out,args.physics)
 except Exception as exc:guard=type(exc).__name__+': '+str(exc)
 finally:mujoco.mj_step=original
 ib_metadata=ib_session.save(out) if ib_session is not None else None
 aff_metadata=aff_session.save(out) if aff_session is not None else None
 diagnostics=[]
 for k,(backend,names,p) in models.items():
  z=log[k];n=z['n'];np.savez_compressed(out/'backend_trace.npz',**{name:value[:n] for name,value in z.items() if name not in ('n','knee_max')},muscle_names=np.array(names));r=dict(samples=n,maximum_knee_extension_deg=np.degrees(z['knee_max']).tolist(),strict_no_hyperextension=bool(np.all(z['knee_max']<=0)),physical_parameters=p.tolist())
  if backend is not None:
   r.update(backend.diagnostics());bound=np.flatnonzero(backend.bank.bounds);r['final_bound_actuators']=[names[i] for i in bound];r['final_bound_projected_fiber_over_l0']=[float(backend.bank.proposed[i]/p[i,1]) for i in bound]
  diagnostics.append(r)
 result=dict(manifest=manifest,guard=guard,total_wall_s=time.perf_counter()-tic,backend_diagnostics=diagnostics,protocol_assessment=None if protocol_result is None else protocol_result.get('assessment'),complete=bool(guard is None and protocol_result is not None and protocol_result.get('failure') is None and protocol_result.get('assessment',{}).get('complete',True)),scope='Experimental live SNS force-backend comparison; no retuning and no new walking inference beyond the actual protocol gates.')
 result['source_ib_pathway']=ib_metadata
 result['source_afferent_pathway']=aff_metadata
 (out/'BACKEND_RESULT.json').write_text(json.dumps(result,indent=2,default=float));print('BACKEND_RESULT',json.dumps(result,default=float),flush=True)
 if guard:raise RuntimeError(guard)
if __name__=='__main__':main()
