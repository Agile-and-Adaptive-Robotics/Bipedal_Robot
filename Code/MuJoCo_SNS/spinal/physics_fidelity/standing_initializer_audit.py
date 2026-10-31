"""Isolated standing initializer sign/row/projection audit; no runner edits.

Corrected prototype is bounded nonnegative least squares in anatomical joint
space with equality projection, actual fixed-pose active capacity, native
muscle passive forces, and optional bilateral activation ties. Constraint
forces are frozen for the linear solve; recomputed dynamics are reported too.
"""
import os,sys,json
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parents[1];ROOT=H.parents[2];sys.path[:0]=[str(H),str(H/'functional_gait2392')]
import numpy as np,mujoco,runner as R
from scipy.optimize import lsq_linear
from scipy.optimize import nnls
from mechanics import CoupledArms,JOINTS
from build_network import _group_weight
from muscle_map import classify
from params import POSTURE_OVERRIDE,W_POSTURE
OUT=Path(__file__).resolve().parent/'runs/audit_20261030';OUT.mkdir(parents=True,exist_ok=True)

def foot_floor_gap(m,d):
 pg=m.geom('ground-plane').id;normal=d.geom_xmat[pg].reshape(3,3)[:,2].copy()
 if normal[2]<0:normal=-normal
 gaps=[]
 for g in range(m.ngeom):
  if m.body(int(m.geom_bodyid[g])).name not in ('calcn_r','calcn_l','toes_r','toes_l'):continue
  if m.geom_type[g]!=mujoco.mjtGeom.mjGEOM_MESH:continue
  mesh=m.geom_dataid[g];start=m.mesh_vertadr[mesh];num=m.mesh_vertnum[mesh]
  vert=m.mesh_vert[start:start+num];gaps.append(float(np.min(vert@(d.geom_xmat[g].reshape(3,3).T@normal))+(d.geom_xpos[g]-d.geom_xpos[pg])@normal))
 return min(gaps)

def build(pose):
 raw=mujoco.MjModel.from_xml_path(str(R.MODEL));d=mujoco.MjData(raw);mujoco.mj_resetDataKeyframe(raw,d,0);ar=CoupledArms(raw)
 d.qpos[ar.addresses]=0.
 for name in ('pelvis_tilt','pelvis_list','pelvis_rotation'):d.qpos[raw.joint(name).qposadr[0]]=0.
 if pose=='symmetric_bent':
  for side in ('r','l'):
   for name,value in [('hip_flexion',-7.),('knee_angle',-10.),('ankle_angle',0.),('hip_adduction',np.degrees(.04))]:d.qpos[raw.joint(name+'_'+side).qposadr[0]]=np.radians(value)
 ar.project(d);mujoco.mj_forward(raw,d)
 slides=[raw.joint('pelvis_'+s).id for s in ('tx','ty','tz')];j=max(slides,key=lambda j:abs(d.xaxis[j,2]))
 d.qpos[raw.jnt_qposadr[j]]+=(-.001-foot_floor_gap(raw,d))/d.xaxis[j,2]
 ar.project(d);mujoco.mj_forward(raw,d);pose_values=R.capture_pose(raw,d)
 m=R.apply_harness(raw,d,kxy=1500.);d=mujoco.MjData(m);R.seed_pose(m,d,pose_values);mujoco.mj_forward(m,d)
 return m,d

report={}
for pose in ('neutral','symmetric_bent'):
 m,d=build(pose);names=[m.actuator(i).name for i in range(m.nu)];ar=CoupledArms(m)
 d.qvel[:]=0.;d.act[:]=0.;mujoco.mj_forward(m,d);ar.update(d);E=ar.E.copy();moment=d.actuator_moment.copy()
 passive=d.actuator_force.copy();passive_general=moment.T@passive
 demand=(d.qfrc_bias-d.qfrc_constraint-d.qfrc_passive).copy();target=E.T@(demand-passive_general)
 rawA=moment.T;selected=np.flatnonzero(np.max(abs(rawA),axis=1)>1e-9);preload_dest=[]
 for side in ('r','l'):
  for joint,value in [('hip_flexion_'+side,-40.),('knee_angle_'+side,-60.),('ankle_angle_'+side,40.)]:
   index=int(m.joint(joint).qposadr[0]);dest=int(selected[index]) if index<len(selected) else None
   preload_dest.append(dict(intended_joint=joint,legacy_added_Nm=value,cropped_index=index,actual_full_dof=dest,actual_joint=m.joint(int(m.dof_jntid[dest])).name if dest is not None else None))
 joint='lumbar_extension';index=int(m.joint(joint).qposadr[0]);dest=int(selected[index]) if index<len(selected) else None
 preload_dest.append(dict(intended_joint=joint,legacy_added_Nm=-20.,cropped_index=index,actual_full_dof=dest,actual_joint=m.joint(int(m.dof_jntid[dest])).name if dest is not None else None))
 d.act[:]=1.;mujoco.mj_forward(m,d);capacity=-(d.actuator_force-passive);A=E.T@(-moment.T*capacity[None,:])
 prior=np.array([min(.5,.6*_group_weight(classify(a),W_POSTURE,POSTURE_OVERRIDE)) for a in names])
 variants={};activations={}
 for preload in (False,True):
  d.act[:]=0.;mujoco.mj_forward(m,d)
  x=R.solve_standing_activations(m,d,np.maximum(m.actuator_gainprm[:,2],5.),preload=preload)
  activations['legacy_preload' if preload else 'legacy_no_preload']=x
 # Exact legacy algebra with only the pulling-force sign corrected.
 Fmax=np.maximum(m.actuator_gainprm[:,2],5.)
 signed=-rawA[selected]*Fmax[None,:]
 xx,_=nnls(np.vstack([.01*signed,.25*np.eye(m.nu)]),np.r_[.01*demand[selected],.25*prior])
 activations['legacy_sign_only_no_preload']=np.clip(xx,0.,1.)
 for tied in (False,True):
  groups=[];seen=set()
  for i,name in enumerate(names):
   if name in seen:continue
   other=name[:-2]+('_l' if name.endswith('_r') else '_r') if name.endswith(('_r','_l')) else name
   indices=[i]+([names.index(other)] if other!=name and other in names else [])
   for j in indices:seen.add(names[j])
   groups.append(indices)
  S=np.zeros((m.nu,len(groups))) if tied else np.eye(m.nu)
  if tied:
   for k,ids in enumerate(groups):S[ids,k]=1.
  # Same per-muscle ridge as legacy; bounds solved directly rather than clipped.
  fit=lsq_linear(np.vstack([.01*A@S,.25*S]),np.r_[.01*target,.25*prior],bounds=(0.,1.),method='bvls')
  activations['corrected_tied' if tied else 'corrected_free']=S@fit.x
 from standing_solver import solve_standing_activations as corrected_callable
 prototype,prototype_diagnostics=corrected_callable(m,d,preload=False,bilateral=True,return_diagnostics=True)
 assert np.max(abs(prototype-activations['corrected_tied']))<1e-8
 saved_wrench=d.xfrc_applied.copy();d.xfrc_applied[m.body('pelvis').id,2]=.7*float(m.body_mass.sum())*9.81
 support_x,support_diagnostics=corrected_callable(m,d,preload=False,bilateral=True,return_diagnostics=True)
 assert support_diagnostics['incoming_external_wrench_present']
 d.xfrc_applied[:]=saved_wrench;mujoco.mj_forward(m,d)
 for name,x in activations.items():
  predicted=A@x;residual=predicted-target
  d.act[:]=x;d.qvel[:]=0.;mujoco.mj_forward(m,d)
  net=d.qfrc_actuator+d.qfrc_passive+d.qfrc_constraint-d.qfrc_bias
  actual=E.T@net
  asym=[abs(x[i]-x[names.index(n[:-2]+'_l')]) for i,n in enumerate(names) if n.endswith('_r') and n[:-2]+'_l'in names]
  variants[name]=dict(frozen_constraint_residual_L2_Nm=float(np.linalg.norm(residual)),frozen_constraint_residual_max_Nm=float(np.max(abs(residual))),
   recomputed_projected_net_force_L2_Nm=float(np.linalg.norm(actual)),max_qacc=float(np.max(abs(d.qacc))),max_selected_anatomical_qacc=float(np.max(abs(d.qacc[ar.rows]))),max_bilateral_activation_difference=float(max(asym)),
   activation_range=[float(x.min()),float(x.max())],activation_by_muscle=dict(zip(names,x.tolist())),
   residual_by_joint=dict(zip(JOINTS,residual.tolist())),recomputed_net_force_by_joint=dict(zip(JOINTS,actual.tolist())))
 # How misleading is the legacy positive sign under its own fitted activations?
 old=activations['legacy_no_preload'];legacy_design=rawA*np.maximum(m.actuator_gainprm[:,2],5.)[None,:]
 design_residual=legacy_design[selected]@old-demand[selected]
 report[pose]=dict(ground_contacts=int(d.ncon),minimum_foot_mesh_gap_m=foot_floor_gap(m,d),fitted_rows=len(selected),nv=m.nv,nu=m.nu,
  callable_diagnostics=prototype_diagnostics,existing_external_support_test=support_diagnostics,
  preload_destinations=preload_dest,legacy_own_design_residual_L2=float(np.linalg.norm(design_residual)),
  anatomical_target_Nm=dict(zip(JOINTS,target.tolist())),native_passive_muscle_torque_Nm=dict(zip(JOINTS,(E.T@passive_general).tolist())),variants=variants)
 print(pose,{n:(v['frozen_constraint_residual_L2_Nm'],v['recomputed_projected_net_force_L2_Nm']) for n,v in variants.items()},flush=True)
result=dict(scope='Isolated fixed-pose initializer audit, ground-on rig with 1mm initial mesh penetration. Corrected solve freezes zero-activation constraint forces, then reports independently recomputed dynamics. No integration, walking, or causal curriculum claim. Bilateral activation ties do not imply identical forces under asymmetric surrogate laws.',poses=report)
(OUT/'STANDING_INITIALIZER_AUDIT.json').write_text(json.dumps(result,indent=2))
