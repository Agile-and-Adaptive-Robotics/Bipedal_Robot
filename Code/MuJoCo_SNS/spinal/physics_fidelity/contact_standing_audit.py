"""Actual forward holding tests for root/contact-aware static initialization."""
import os,sys,json,hashlib
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parents[1];ROOT=H.parents[2];sys.path[:0]=[str(H),str(H/'functional_gait2392')]
import numpy as np,mujoco,runner as R
from mechanics import CoupledArms
from standing_solver import solve_standing_activations as oldsolve
from contact_standing_solver import solve_contact_standing_activations as newsolve
from physics_fidelity.ground_alignment import align_via_repaired_probe
HERE=Path(__file__).resolve().parent;OUT=HERE/'runs/contact_standing_v3';OUT.mkdir(parents=True,exist_ok=True)
long_hold='--aligned-long' in sys.argv
if long_hold:OUT=HERE/'runs/contact_standing_v4';OUT.mkdir(parents=True,exist_ok=True)
duration=8. if long_hold else 2.
source=H/'hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/ground_winner_validation_v1/standing.json'
cfg=json.loads(source.read_text());z=np.load(source.with_suffix('.npz'),allow_pickle=True);q0=z['qfull'][0].copy()

def build(force_support,pose_mode):
 raw=mujoco.MjModel.from_xml_path(str(R.MODEL));initial=mujoco.MjData(raw);initial.qpos[:]=q0;mujoco.mj_forward(raw,initial)
 alignment=None
 if pose_mode!='original':alignment=align_via_repaired_probe(raw,initial,R.apply_harness,flatten_ankles=pose_mode.startswith('flat_aligned'),target_gap_m=-.002 if pose_mode.endswith('2mm') else -.001)
 aligned_q=initial.qpos.copy()
 m=R.apply_harness(raw,initial,kxy=1500.)
 profile=json.loads(Path(cfg['force_law_parameters_file']).read_text())
 for family,row in profile['muscles'].items():
  for name in row['actuator_names']:
   i=m.actuator(name).id
   for key in ('gainprm','biasprm','dynprm','lengthrange'):getattr(m,'actuator_'+key)[i]=row[key]
 d=mujoco.MjData(m);d.qpos[:]=aligned_q
 for name in ('pelvis_list','pelvis_rotation','hip_adduction_r','hip_adduction_l','hip_rotation_r','hip_rotation_l','subtalar_angle_r','subtalar_angle_l','lumbar_bending','lumbar_rotation'):
  j=m.joint(name);m.jnt_stiffness[j.id]+=cfg.get('restriction',0.)*800.;m.dof_damping[j.dofadr[0]]+=cfg.get('restriction',0.)*40.
 for side in ('r','l'):
  j=m.joint('knee_angle_'+side);m.jnt_limited[j.id]=True;m.jnt_range[j.id]=np.radians([cfg['knee_training_limits_deg'][0],0.]);m.jnt_margin[j.id]=np.radians(2.);m.jnt_solref[j.id]=[.004,1.];m.jnt_solimp[j.id]=[.999,.9999,.001,.5,2.]
 ar=CoupledArms(m);ar.project(d);mujoco.mj_forward(m,d)
 if force_support:
  for axis in ('tx','ty','tz'):
   j=m.joint('pelvis_'+axis)
   if abs(d.xaxis[j.id,2])>.8:m.jnt_stiffness[j.id]=0.;m.dof_damping[j.dofadr[0]]=0.
  d.xfrc_applied[m.body('pelvis').id,2]=.7*m.body_mass.sum()*9.81
 mujoco.mj_forward(m,d);return m,d,alignment

def force(m,d):
 f=np.zeros(6);total=0.
 for i in range(d.ncon):
  c=d.contact[i]
  if m.geom_type[c.geom1]!=mujoco.mjtGeom.mjGEOM_PLANE and m.geom_type[c.geom2]!=mujoco.mjtGeom.mjGEOM_PLANE:continue
  mujoco.mj_contactForce(m,d,i,f);total+=max(float(f[0]),0.)
 return total

def hold(m,initial,x):
 d=mujoco.MjData(m);d.qpos[:]=initial.qpos;d.xfrc_applied[:]=initial.xfrc_applied;d.qfrc_applied[:]=initial.qfrc_applied;d.act[:]=x;d.ctrl[:]=x;mujoco.mj_forward(m,d)
 ar=CoupledArms(m);initial_acc={m.joint(n).name:float(d.qacc[m.joint(n).dofadr[0]]) for n in ('pelvis_tx','pelvis_ty','pelvis_tz','pelvis_tilt','hip_flexion_r','hip_flexion_l','knee_angle_r','knee_angle_l')}
 bw=float(m.body_mass.sum()*9.81);initial_force=force(m,d)/bw;trace=[];aborted=None
 com_reference=d.subtree_com[0,:2].copy();max_sway=0.;max_pelvis_tilt=abs(float(np.degrees(d.qpos[m.joint('pelvis_tilt').qposadr[0]])))
 runner_reference=np.array([np.mean([d.xpos[m.body('talus_'+s).id,0] for s in ('r','l')]),com_reference[1]])
 runner_sway=0.
 knee_max=np.array([d.qpos[m.joint('knee_angle_'+s).qposadr[0]] for s in ('r','l')])
 for k in range(round(duration/m.opt.timestep)):
  mujoco.mj_step(m,d)
  knee_max=np.maximum(knee_max,[d.qpos[m.joint('knee_angle_'+s).qposadr[0]] for s in ('r','l')])
  max_sway=max(max_sway,float(np.linalg.norm(d.subtree_com[0,:2]-com_reference)))
  runner_sway=max(runner_sway,float(np.linalg.norm(d.subtree_com[0,:2]-runner_reference)))
  max_pelvis_tilt=max(max_pelvis_tilt,abs(float(np.degrees(d.qpos[m.joint('pelvis_tilt').qposadr[0]]))))
  if not np.all(np.isfinite(d.qpos)) or not np.all(np.isfinite(d.qvel)) or d.time<(k+1)*m.opt.timestep-1e-6:aborted='Nonfinite/reset';break
  up=d.xmat[m.body('torso').id].reshape(3,3)[:,1];tilt=float(np.degrees(np.arccos(np.clip(up[2],-1,1))))
  trace.append([d.time,tilt,force(m,d)/bw,*d.qpos[ar.addresses[:4]],d.qpos[m.joint('hip_flexion_l').qposadr[0]]])
  if np.max(abs(d.qvel))>300. or tilt>80.:aborted='Velocity/tilt guard';break
 arr=np.asarray(trace)
 return dict(duration_s=float(d.time),complete=bool(d.time>=duration-1e-6 and aborted is None),aborted=aborted,initial_actual_foot_force_bw=initial_force,
  max_abs_pelvis_tilt_deg=max_pelvis_tilt,max_horizontal_COM_displacement_from_initial_m=max_sway,
  runner_style_COM_reference_xy_m=runner_reference.tolist(),initial_COM_offset_from_runner_reference_m=float(np.linalg.norm(com_reference-runner_reference)),runner_style_sway_max_m=runner_sway,
  maximum_knee_extension_deg=np.degrees(knee_max).tolist(),strict_no_hyperextension=bool(np.all(knee_max<=0.)),
  initial_coordinate_accelerations=initial_acc,max_torso_tilt_deg=float(arr[:,1].max()) if len(arr) else None,
  late_mean_actual_foot_force_bw=float(arr[arr[:,0]>=1.,2].mean()) if np.any(arr[:,0]>=1.) else None,
  final_native_hip_deg=[float(np.degrees(d.qpos[m.joint('hip_flexion_'+s).qposadr[0]])) for s in ('r','l')]),arr

rows={}
for pose_mode in (('flat_aligned','flat_aligned_2mm') if long_hold else ('original','height_aligned','flat_aligned')):
 for support in (False,True):
  m,d,alignment=build(support,pose_mode);old,od=oldsolve(m,d,preload=False,bilateral=True,return_diagnostics=True);new,nd=newsolve(m,d,preload=False,bilateral=True,return_diagnostics=True)
  for label,x,diagnostic in [('internal_only',old,od),('root_contact',new,nd)]:
   name=pose_mode+('_constant70_' if support else '_stage3_spring_')+label;h,trace=hold(m,d,x)
   rows[name]=dict(alignment=alignment,solver=diagnostic,holding=h,activation_by_muscle={m.actuator(i).name:float(v) for i,v in enumerate(x)},npair=m.npair,disableflags=int(m.opt.disableflags),timestep=float(m.opt.timestep))
   np.savez_compressed(OUT/(name+'.npz'),trace=trace);print(name,h,flush=True)
report=dict(source_config=str(source),pose_source='Recorded first standing qpos; zero initial velocity; equality projection enforced.',scope='No SNS. Constant activation/control open-loop hold on reconstructed ground-on solve-style rig. Full support settings specified, no optimized normal forces injected. Not a curriculum evaluation.',cases=rows)
(OUT/'RESULTS.json').write_text(json.dumps(report,indent=2))
