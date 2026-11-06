"""Independent OpenSim verification of generalized insertion-wrench torque."""
from pathlib import Path
import json
import numpy as np
import opensim as o
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/audit_20261105'
source=ROOT/'Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim'
o.Logger.setLevelString('error');m=o.Model(str(source));cases=[('hip_flexion_r','pelvis','femur_r',.3),('hip_adduction_r','pelvis','femur_r',.2),('hip_rotation_r','pelvis','femur_r',.2),('knee_angle_r','femur_r','tibia_r',-.8),('ankle_angle_r','tibia_r','talus_r',.2),('subtalar_angle_r','talus_r','calcn_r',.2),('mtp_angle_r','calcn_r','toes_r',.2)]
for name,parent,child,q in cases:
 a=o.PathActuator();a.setName('test_'+name);a.set_appliesForce(False)
 a.addNewPathPoint('origin',m.updBodySet().get(parent),o.Vec3(-.04,.08,.06))
 a.addNewPathPoint('insertion',m.updBodySet().get(child),o.Vec3(.03,-.12,.025));m.addForce(a)
m.finalizeConnections();s=m.initSystem()
def T(body):
 t=m.getBodySet().get(body).getTransformInGround(s);x=np.eye(4);x[:3,:3]=[[t.R().get(i,j) for j in range(3)] for i in range(3)];x[:3,3]=[t.p().get(i) for i in range(3)];return x
rows=[]
for name,parent,child,q in cases:
 for j in range(m.getCoordinateSet().getSize()):
  c=m.updCoordinateSet().get(j);c.setValue(s,c.getDefaultValue(),False)
 c=m.updCoordinateSet().get(name);c.setValue(s,q,False);m.realizePosition(s)
 rel=np.linalg.inv(T(parent))@T(child);p=np.array([-.04,.08,.06]);ins=np.array([.03,-.12,.025]);endpoint=rel[:3,:3]@ins+rel[:3,3]
 force=p-endpoint;force/=np.linalg.norm(force);forceChild=rel[:3,:3].T@force;moment=np.cross(ins,forceChild)
 h=1e-6;ts=[]
 for value in [q+h,q-h]:
  c.setValue(s,value,False);m.realizePosition(s);ts.append(np.linalg.inv(T(parent))@T(child))
 A=rel[:3,:3].T@((ts[0][:3,:3]-ts[1][:3,:3])/(2*h));omega=np.array([A[2,1]-A[1,2],A[0,2]-A[2,0],A[1,0]-A[0,1]])/2
 velocity=rel[:3,:3].T@((ts[0][:3,3]-ts[1][:3,3])/(2*h))
 predicted=float(moment@omega+forceChild@velocity)
 c.setValue(s,q,False);m.realizePosition(s);act=o.PathActuator.safeDownCast(m.getForceSet().get('test_'+name))
 reference=act.getGeometryPath().computeMomentArm(s,c)
 assert abs(predicted-reference)<1e-6,(name,predicted,reference)
 rows.append(dict(coordinate=name,coordinate_axis_in_child=omega.tolist(),translation_velocity_in_child=velocity.tolist(),predicted_moment_arm_m=predicted,opensim_moment_arm_m=reference,error_m=abs(predicted-reference)))
(OUT/'coordinate_wrench_opensim_verification.json').write_text(json.dumps(rows,indent=2))
print('PASS: generalized wrench equals native OpenSim moment arms for all seven right-leg coordinates. Max error:',max(r['error_m'] for r in rows))
