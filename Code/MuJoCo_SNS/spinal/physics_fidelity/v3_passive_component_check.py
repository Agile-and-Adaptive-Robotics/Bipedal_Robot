"""Separate source passive/active forces and knee-torque bias at low drive."""
import os,sys,json,math
from pathlib import Path
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parents[1];ROOT=H.parents[2];sys.path[:0]=[str(H),str(H/'functional_gait2392')]
import numpy as np,mujoco,opensim as osim,bsolve_ik as B
from mechanics import CoupledArms,JOINTS
HERE=Path(__file__).resolve().parent;OUT=HERE/'runs/audit_20261030'
spec=json.loads((HERE/'shared_force_law_functional_v3.json').read_text());signs={r['coordinate']:r['sign'] for r in json.loads((H/'hunt_training/runs/ik_audit/fk_signs.json').read_text())}
v2spec=json.loads((HERE/'shared_force_law_functional_v2.json').read_text())
om=osim.Model(str(ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/mjc/gait2392_simbody/gait2392_simbody.osim'));om.setUseVisualizer(False);s=om.initSystem();cs=om.getCoordinateSet();muscles=om.getMuscles()
m,d=B.build_air_model(False);ar=CoupledArms(m);families=('bifemsh','semimem','semiten','bifemlh','rect_fem','tfl');rows={}
for label,hip,knee in [('neutral',0.,0.),('bent_standing',7.,-10.),('swing',30.,-60.)]:
 for name in JOINTS:
  base=name.rsplit('_',1)[0];value=math.radians(hip if base=='hip_flexion' else knee if base=='knee_angle' else 0.)
  cs.get(name).setValue(s,value,False);d.qpos[m.joint(name).qposadr[0]]=value*signs[name]
 ar.project(d);mujoco.mj_forward(m,d);A=ar.update(d).copy()
 for i in range(muscles.getSize()):muscles.get(i).setActivation(s,.01)
 om.equilibrateMuscles(s);om.realizeDynamics(s)
 entries=[]
 for family in families:
  for side in ('r','l'):
   name=family+'_'+side;aid=m.actuator(name).id;L=float(d.actuator_length[aid]);ms=muscles.get(name)
   source_passive=float(ms.getPassiveFiberForce(s)*ms.getCosPennationAngle(s));source_active=float(ms.getActiveFiberForce(s)*ms.getCosPennationAngle(s))
   gp=m.actuator_gainprm[aid,:9];bp=m.actuator_biasprm[aid,:9];lr=m.actuator_lengthrange[aid]
   original_passive=-mujoco.mju_muscleBias(L,lr,1.,bp)
   v=spec['muscles'][family];gp2=np.asarray(v['gainprm'][:9]);bp2=np.asarray(v['biasprm'][:9]);lr2=np.asarray(v['lengthrange'])
   shared_passive=-mujoco.mju_muscleBias(L,lr2,1.,bp2)
   v2=v2spec['muscles'].get(family)
   v2_passive=-mujoco.mju_muscleBias(L,np.asarray(v2['lengthrange']),1.,np.asarray(v2['biasprm'][:9])) if v2 else original_passive
   k=JOINTS.index('knee_angle_'+side);rho=float(A[k,aid])
   entries.append(dict(muscle=name,source_passive_tension_N=source_passive,source_active_tension_at_activation_001_N=source_active,
    original_mujoco_passive_tension_N=float(original_passive),shared_mujoco_passive_tension_N=float(shared_passive),
    v2_mujoco_passive_tension_N=float(v2_passive),v2_mujoco_passive_knee_torque_Nm=float(rho*v2_passive),
    original_mujoco_passive_knee_torque_Nm=float(rho*original_passive),shared_mujoco_passive_knee_torque_Nm=float(rho*shared_passive),
    source_passive_knee_torque_Nm=float(ms.computeMomentArm(s,cs.get('knee_angle_'+side))*source_passive)))
 rows[label]=entries
report=dict(scope='Source components at activation0.01, not zero activation; rigid MuJoCo passive force is activation-independent, whereas compliant source fiber state can change with activation. Knee torques use each engine moment arm. No dynamics or causal gait claim.',poses=rows)
(OUT/'V3_PASSIVE_COMPONENT_CHECK.json').write_text(json.dumps(report,indent=2))
for pose,points in rows.items():
 for side in ('r','l'):
  p=[x for x in points if x['muscle'].endswith('_'+side) and x['muscle'].split('_')[0] in ('bifemsh','semimem','semiten','bifemlh')]
  print(pose,side,'hamstring_passive_knee_Nm source/original/v2/v3',*(sum(x[k] for x in p) for k in ('source_passive_knee_torque_Nm','original_mujoco_passive_knee_torque_Nm','v2_mujoco_passive_knee_torque_Nm','shared_mujoco_passive_knee_torque_Nm')))
