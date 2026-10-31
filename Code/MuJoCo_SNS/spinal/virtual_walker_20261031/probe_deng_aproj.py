"""Read-only physical torque audit of the locally archived original rat model."""
import importlib.util,json
from pathlib import Path
import numpy as np
from import_virtual_walker import REPO,HERE
spec=importlib.util.spec_from_file_location('legacy_readonly',REPO/'Code/MuJoCo_SNS/spinal/w2l_mujoco/make_w2l_mjcf.py');G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G)
p=REPO/'Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination/runs/deng_source_audit_20261031/K_A synergy.aproj.txt';model=G.parse_aproj(str(p));body=model['by_id'];rows=[]
for side in ('LH','RH'):
 foot=next(b for b in model['bodies'] if b.name==side+'_Foot' and b.joint);joint=foot.joint;anchor=np.array(foot.p_world)+np.array(foot.R_world)@joint['pos_local_al'];axis=np.array(foot.R_world)@np.array(joint['R_local_al'])@np.array([-1.,0,0])
 def distal(a):
  while a is not None:
   if a is foot:return True
   a=a.parent
  return False
 for muscle in model['muscles']:
  if not muscle['name'].startswith(side+'_AnkleZ'):continue
  attachments=[body[i] for i in muscle['attach_ids']];torque=np.zeros(3)
  for a,b in zip(attachments[:-1],attachments[1:]):
   pa,pb=np.array(a.p_world),np.array(b.p_world);force=(pb-pa)/np.linalg.norm(pb-pa)
   if distal(a):torque+=np.cross(pa-anchor,force)
   if distal(b):torque+=np.cross(pb-anchor,-force)
  rows.append(dict(muscle=muscle['name'],muscle_id=muscle['id'],attachment_names=[a.name for a in attachments],attachment_parents=[a.parent.name for a in attachments],attachment_world_m=[a.p_world for a in attachments],ankle_world_position_m=anchor.tolist(),native_positive_joint_axis=axis.tolist(),world_torque_Nm_per_N=torque.tolist(),native_positive_joint_torque_per_N=float(torque@axis),foot_world_rotation=foot.R_world,foot_dims_xyz_m=[foot.lwh[0],foot.lwh[2],foot.lwh[1]]))
out=dict(source=str(p),frame='Original AnimatLab source world, before anatomical frame adaptation.',scope='Unit-tension physical torque only; PF/MN signals and anatomical dorsiflexion direction require separate identification.',muscles=rows);dest=HERE/'runs'/'DENG2019_TORQUE_PROBE.json';dest.write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
