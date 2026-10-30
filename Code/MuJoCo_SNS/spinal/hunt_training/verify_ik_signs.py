"""Forward-kinematics sign adapter measured independently in both engines."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
H=Path(__file__).resolve().parent;S=H.parent;sys.path[:0]=[str(S),str(S/'functional_gait2392')]
import numpy as np,opensim as osim,mujoco,bsolve_ik as B
from scipy.spatial.transform import Rotation
from mechanics import CoupledArms,JOINTS
root=S.parents[2]/'Solid_Models/OpenSim/Gait2392_Robotbody'
om=osim.Model(str(root/'subject01_simbody.osim'));state=om.initSystem();cs=om.getCoordinateSet()
for i in range(cs.getSize()):
 c=cs.get(i)
 if c.getLocked(state):c.setLocked(state,False)
 c.setValue(state,0.)
om.realizePosition(state)
def rotation(body):
 r=om.getBodySet().get(body).getTransformInGround(state).R().asMat33()
 return np.array([[r.get(i,j) for j in range(3)] for i in range(3)])
mm,md=B.build_air_model(False);arms=CoupledArms(mm);md.qpos[arms.addresses]=0.;arms.project(md);mujoco.mj_kinematics(mm,md)
adapter=np.array([[1.,0.,0.],[0.,0.,-1.],[0.,1.,0.]])
rows=[]
for joint,jid,adr in zip(JOINTS,arms.jids,arms.addresses):
 body='torso' if joint.startswith('lumbar') else ('femur_' if joint.startswith('hip') else 'tibia_' if joint.startswith('knee') else 'talus_' if joint.startswith('ankle') else 'calcn_' if joint.startswith('subtalar') else 'toes_')+joint[-1]
 base_os=rotation(body);bid=mujoco.mj_name2id(mm,mujoco.mjtObj.mjOBJ_BODY,body);base_mj=md.xmat[bid].reshape(3,3).copy()
 cs.get(joint).setValue(state,.1);om.realizePosition(state);md.qpos[adr]=.1;arms.project(md);mujoco.mj_kinematics(mm,md)
 axis_os=Rotation.from_matrix(rotation(body)@base_os.T).as_rotvec()/.1
 axis_mj=adapter.T@Rotation.from_matrix(md.xmat[bid].reshape(3,3)@base_mj.T).as_rotvec()/.1
 dot=float(axis_os@axis_mj/(np.linalg.norm(axis_os)*np.linalg.norm(axis_mj)))
 rows.append(dict(coordinate=joint,source_axis=axis_os.tolist(),mujoco_axis_in_isb=axis_mj.tolist(),axis_dot=dot,sign=1 if dot>0 else -1))
 cs.get(joint).setValue(state,0.);om.realizePosition(state);md.qpos[adr]=0.;arms.project(md);mujoco.mj_kinematics(mm,md)
assert all(abs(r['axis_dot'])>.99 for r in rows),rows
out=H/'runs/ik_audit/fk_signs.json';out.write_text(json.dumps(rows,indent=2))
print(json.dumps(rows,indent=2))
