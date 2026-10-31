"""Read-only finite differences at the saved near-flat W2L108 pose."""
import os,sys,json
from pathlib import Path
os.environ.setdefault('CONDA_PREFIX',r'D:\Anaconda\envs\myo')
import numpy as np
import mujoco
B=Path(__file__).resolve().parent;V=B.parent/'virtual_walker_20261031';O=B/'runs/canonical_12II/ik_reference'
sys.path.insert(0,str(V));import native_kernel as K
S=V/'runs/W2L_VirtualWalker_Newer';meta=json.loads((S/'IMPORT_AUDIT.json').read_text())
m=mujoco.MjModel.from_xml_path(str(S/'canonical108_ground_baseline_v1/model.xml'));d=mujoco.MjData(m)
q=np.load(O/'flat_existing_rails.npz')['qpos'];d.qpos[:]=q;mujoco.mj_forward(m,d)
def ID(kind,s):return mujoco.mj_name2id(m,kind,s)
def body(s):return ID(mujoco.mjtObj.mjOBJ_BODY,s)
def joint(s):return ID(mujoco.mjtObj.mjOBJ_JOINT,s)
def coordinates(side):
 h=d.xanchor[joint('hip_'+side)].copy();k=d.xanchor[joint('knee_'+side)].copy();a=d.xanchor[joint('ankle_'+side)].copy()
 thigh=k-h;shank=a-k
 hp=np.arctan2(thigh[0],-thigh[1]);sp=np.arctan2(shank[0],-shank[1])
 F=d.xmat[body('foot_'+side)].reshape(3,3);T=d.xmat[body('toe_'+side)].reshape(3,3)
 fp=np.arctan2(F[1,0],F[0,0]);tp=np.arctan2(T[1,0],T[0,0])
 return np.array([hp,hp-sp,fp-sp,tp-fp]),dict(hip=h.tolist(),knee=k.tolist(),ankle=a.tolist(),foot_center=d.xpos[body('foot_'+side)].tolist(),toe_center=d.xpos[body('toe_'+side)].tolist())
labels=['hip_flexion_from_knee_anterior_to_hip','knee_flexion_from_thigh_minus_shank_anterior_tilt','ankle_dorsiflexion_from_foot_pitch_minus_shank_anterior_tilt','toe_extension_from_toe_minus_foot_pitch']
rows=[];eps=1e-5
for side in ['L','R']:
 d.qpos[:]=m.qpos0;mujoco.mj_forward(m,d);home,_=coordinates(side)
 for k,kind in enumerate(['hip','knee','ankle','toe']):
  j=joint(kind+'_'+side);qa=m.jnt_qposadr[j];va=m.jnt_dofadr[j]
  d.qpos[:]=q;mujoco.mj_forward(m,d);actual,points=coordinates(side)
  d.qpos[qa]+=eps;mujoco.mj_forward(m,d);plus,pplus=coordinates(side)
  d.qpos[qa]-=2*eps;mujoco.mj_forward(m,d);minus,pminus=coordinates(side)
  row=dict(joint=kind+'_'+side,native_q_deg=float(np.degrees(q[qa])),anatomical_coordinate=labels[k],home_offset_deg=float(np.degrees(home[k])),anatomical_at_reference_deg=float(np.degrees(actual[k])),positive_q_anatomical_derivative=float((plus[k]-minus[k])/(2*eps)),knee_position_derivative_m_per_rad=((np.array(pplus['knee'])-pminus['knee'])/(2*eps)).tolist(),ankle_position_derivative_m_per_rad=((np.array(pplus['ankle'])-pminus['ankle'])/(2*eps)).tolist(),toe_center_derivative_m_per_rad=((np.array(pplus['toe_center'])-pminus['toe_center'])/(2*eps)).tolist())
  rows.append(row)
d.qpos[:]=q;mujoco.mj_forward(m,d);controller=K.NativeController(meta['source_snapshot'],meta['muscles']);muscles=[]
for pre,target,co in controller.outputs:
 b=meta['muscles'][target];a=ID(mujoco.mjtObj.mjOBJ_ACTUATOR,b['name'].replace(' ','_'))
 actions=[]
 for row in rows:
  va=m.jnt_dofadr[joint(row['joint'])];v=float(d.actuator_moment[a,va])
  if abs(v)>1e-5:actions.append(dict(joint=row['joint'],moment_arm_m=v,physical_action='positive_'+row['anatomical_coordinate'] if v>0 else 'negative_'+row['anatomical_coordinate']))
 muscles.append(dict(MN=controller.names[pre],muscle=b['name'],adapter=co.tolist(),actions=actions))
out=dict(method='Central finite differences 1e-5rad of actual hip/knee/ankle anchors and foot/toe body orientations at exact saved qpos, plus actual pulling actuator Jacobian and source MN bindings.',coordinate_labels=labels,joints=rows,MN_routes=muscles)
(O/'FK_ACTION_AUDIT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
