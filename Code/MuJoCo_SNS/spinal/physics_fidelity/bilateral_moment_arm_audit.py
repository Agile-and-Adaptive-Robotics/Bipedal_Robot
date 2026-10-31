"""Bilateral matched-coordinate, equality-aware MuJoCo/OpenSim route audit."""
import os,sys,json,math
from pathlib import Path
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parents[1];ROOT=H.parents[2]
sys.path[:0]=[str(H),str(H/'functional_gait2392')]
import numpy as np,mujoco,opensim as osim,bsolve_ik as B
from mechanics import CoupledArms,JOINTS
OUT=Path(__file__).resolve().parent/'runs/audit_20261030';OUT.mkdir(parents=True,exist_ok=True)
SOURCE=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'
om=osim.Model(str(SOURCE));om.setUseVisualizer(False);s=om.initSystem();coords=om.getCoordinateSet();muscles=om.getMuscles()
m,d=B.build_air_model(False);ar=CoupledArms(m)
SIGN_SOURCE=H/'hunt_training/runs/ik_audit/fk_signs.json'
signs={r['coordinate']:r['sign'] for r in json.loads(SIGN_SOURCE.read_text())}
families=('bifemsh','vas_med','vas_int','vas_lat','soleus','tib_ant','rect_fem','semimem','semiten','bifemlh','med_gas','lat_gas')
axes=('hip_flexion','hip_adduction','hip_rotation','knee_angle','ankle_angle','subtalar_angle','mtp_angle')
rows={}
for pose,angles in [('neutral',{}),('bilateral_swing',dict(hip_flexion=30.,knee_angle=-60.,ankle_angle=10.))]:
 d.qpos[ar.addresses]=0.
 for name in JOINTS:
  value=math.radians(angles.get(name.rsplit('_',1)[0],0.))
  coords.get(name).setValue(s,value,False)
  d.qpos[m.joint(name).qposadr[0]]=value*signs[name]
 ar.project(d);mujoco.mj_forward(m,d);A=ar.update(d).copy();fd=ar.finite_difference(d)
 om.realizePosition(s)
 d.act[:]=0.;mujoco.mj_forward(m,d);passive_force=d.actuator_force.copy()
 d.act[:]=1.;mujoco.mj_forward(m,d);active_force=np.abs(d.actuator_force-passive_force)
 d.act[:]=0.;mujoco.mj_forward(m,d)
 entries=[];lengths=[];capacities=[]
 for family in families:
  for side in ('r','l'):
   name=family+'_'+side;aid=m.actuator(name).id;muscle=muscles.get(name)
   lengths.append(dict(muscle=name,mujoco_length_m=float(d.actuator_length[aid]),opensim_length_m=float(muscle.getLength(s))))
   gp=m.actuator_gainprm[aid];lr=m.actuator_lengthrange[aid];L0=(lr[1]-lr[0])/(gp[1]-gp[0]);Lt=lr[0]-L0*gp[0]
   capacities.append(dict(muscle=name,active_capacity_N=float(active_force[aid]),passive_force_N=float(abs(passive_force[aid])),fitted_Fmax_N=float(gp[2]),source_Fmax_N=float(muscle.getMaxIsometricForce()),fitted_fiber_optimum_m=float(L0),fitted_tendon_offset_m=float(Lt),normalized_fiber_length=float((d.actuator_length[aid]-Lt)/L0),gainprm=gp.tolist(),actuator_lengthrange=lr.tolist()))
   for axis in axes:
    joint=axis+'_'+side;k=JOINTS.index(joint)
    om_arm=float(muscle.computeMomentArm(s,coords.get(joint)));mj_arm=float(A[k,aid])*signs[joint];f=float(fd[k,aid])*signs[joint]
    entries.append(dict(muscle=name,joint=joint,mujoco_arm_m=mj_arm,mujoco_native_arm_m=float(A[k,aid]),coordinate_sign=signs[joint],opensim_arm_m=om_arm,finite_difference_arm_m=f,
      absolute_error_mm=1000*abs(mj_arm-om_arm),relative_error=abs(mj_arm-om_arm)/max(abs(om_arm),.001),
      significant_sign_disagreement=bool(abs(mj_arm)>.001 and abs(om_arm)>.001 and mj_arm*om_arm<0)))
 bilateral=[]
 for family in families:
  for axis in axes:
   r=next(x for x in entries if x['muscle']==family+'_r' and x['joint']==axis+'_r')
   l=next(x for x in entries if x['muscle']==family+'_l' and x['joint']==axis+'_l')
   bilateral.append(dict(muscle_family=family,joint_axis=axis,mujoco_right_m=r['mujoco_arm_m'],mujoco_left_m=l['mujoco_arm_m'],
    mujoco_bilateral_difference_mm=1000*abs(r['mujoco_arm_m']-l['mujoco_arm_m']),
    opensim_right_m=r['opensim_arm_m'],opensim_left_m=l['opensim_arm_m'],opensim_bilateral_difference_mm=1000*abs(r['opensim_arm_m']-l['opensim_arm_m'])))
 force_bilateral=[]
 for family in families:
  r=next(x for x in capacities if x['muscle']==family+'_r');l=next(x for x in capacities if x['muscle']==family+'_l')
  force_bilateral.append(dict(muscle_family=family,right_active_capacity_N=r['active_capacity_N'],left_active_capacity_N=l['active_capacity_N'],left_to_right_active_capacity_ratio=l['active_capacity_N']/max(r['active_capacity_N'],1e-9),right_normalized_fiber_length=r['normalized_fiber_length'],left_normalized_fiber_length=l['normalized_fiber_length']))
 rows[pose]=dict(angles_deg=angles,entries=entries,lengths=lengths,bilateral=bilateral,capacities=capacities,force_bilateral=force_bilateral,
  max_mujoco_finite_difference_error_m=float(np.max(abs(A-fd))),
  significant_sign_disagreements=[x for x in entries if x['significant_sign_disagreement']],
  sizeable_cross_engine_errors=[x for x in entries if x['absolute_error_mm']>5 and x['relative_error']>.2],
  sizeable_bilateral_mujoco_differences=[x for x in bilateral if x['mujoco_bilateral_difference_mm']>2])
report=dict(source_opensim=str(SOURCE),source_mujoco=str(B.runner.MODEL),coordinate_sign_source=str(SIGN_SOURCE),coordinate_signs=signs,opensim_version=osim.GetVersion(),mujoco_version=mujoco.__version__,
 scope='Neutral and identical bilateral anatomical hip30/knee-60/ankle10; all other selected coordinates zero. Both q and native moment arms mapped using validated FK signs. Moment arm is torque per positive muscle tension, -dL/dq, m/rad. Native MJ world x anterior,y left,z superior versus OpenSim x anterior,y superior,z right. No scalar Fmax normalization enters moment arms.',
 poses=rows)
(OUT/'BILATERAL_MOMENT_ARMS.json').write_text(json.dumps(report,indent=2))
for pose,r in rows.items():
 print(pose,'FDerror',r['max_mujoco_finite_difference_error_m'],'signerrors',len(r['significant_sign_disagreements']),'largebilateral',len(r['sizeable_bilateral_mujoco_differences']))
 for x in r['sizeable_bilateral_mujoco_differences']:print('BILATERAL',x)
 for x in r['sizeable_cross_engine_errors']:
  if x['joint'].startswith(('knee_angle','ankle_angle')):print('ENGINE',x)
