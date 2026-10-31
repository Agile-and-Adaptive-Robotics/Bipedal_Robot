"""Equality-aware armature reflection, without changing follower physics."""
import os,sys,json
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parents[1];ROOT=H.parents[2]
sys.path[:0]=[str(H),str(H/'functional_gait2392')]
import numpy as np,mujoco,bsolve_ik as B
from mechanics import CoupledArms,JOINTS
OUT=Path(__file__).resolve().parent/'runs/audit_20261030';OUT.mkdir(parents=True,exist_ok=True)
rows={}
for mode in (False,True):
 m,d=B.build_air_model(mode);ar=CoupledArms(m)
 follower=np.unique(ar.dependent);mass=np.zeros((m.nv,m.nv));velocity=np.zeros(m.nv)
 cases={}
 for label,pose in [('neutral',{}),('swing_right',{'hip_flexion_r':30.,'knee_angle_r':-60.,'ankle_angle_r':10.})]:
  d.qpos[ar.addresses]=0.
  for name,value in pose.items():d.qpos[m.joint(name).qposadr[0]]=np.radians(value)
  ar.project(d);mujoco.mj_forward(m,d);ar.update(d)
  E=ar.E.copy();mask=np.zeros(m.nv);mask[follower]=m.dof_armature[follower]
  reflected=E.T@(mask[:,None]*E)
  mujoco.mj_fullM(m,mass,d.qM);full=E.T@mass@E
  saved=d.qpos.copy();finite=np.zeros_like(E);eps=1e-6
  for k,adr in enumerate(ar.addresses):
   d.qpos[:]=saved;d.qpos[adr]-=eps;ar.project(d);minus=d.qpos.copy()
   d.qpos[:]=saved;d.qpos[adr]+=eps;ar.project(d);plus=d.qpos.copy()
   mujoco.mj_differentiatePos(m,velocity,2*eps,minus,plus);finite[:,k]=velocity
  d.qpos[:]=saved
  error=float(np.max(abs(finite-E)));assert error<1e-6,error
  cases[label]=dict(joint_order=list(JOINTS),follower_armature_matrix_kg_m2=reflected.tolist(),follower_armature_diagonal_kg_m2={name:float(reflected[i,i]) for i,name in enumerate(JOINTS)},
      total_reduced_mass_diagonal_kg_m2={name:float(full[i,i]) for i,name in enumerate(JOINTS)},
      follower_fraction_of_reduced_diagonal={name:float(reflected[i,i]/full[i,i]) for i,name in enumerate(JOINTS)},
      maximum_projection_finite_difference_error=error)
 rows['id_mode' if mode else 'forward_mode']=dict(follower_count=len(follower),follower_armature_values=np.unique(m.dof_armature[follower]).tolist(),equalities_disabled=bool(m.opt.disableflags & int(mujoco.mjtDisableBit.mjDSBL_EQUALITY)),poses=cases)
report=dict(scope='Kinematic reduction E.T M E with pelvis fixed; not trajectory torques or an OpenSim inertia match. ID flags/armatures inspected only, not adopted.',modes=rows)
(OUT/'FOLLOWER_INERTIA.json').write_text(json.dumps(report,indent=2))
for mode,r in rows.items():
 for pose,p in r['poses'].items():print(mode,pose,'knee_r',p['follower_armature_diagonal_kg_m2']['knee_angle_r'],'fraction',p['follower_fraction_of_reduced_diagonal']['knee_angle_r'],'knee_l',p['follower_armature_diagonal_kg_m2']['knee_angle_l'])
