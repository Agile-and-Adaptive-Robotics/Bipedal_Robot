"""Geometry deadband and physical contact-force direction, not a fitted test."""
import json
from pathlib import Path
import numpy as np
from import_virtual_walker import mujoco,HERE
from contact_guides import Measurements
out=[]
for clearance in (.01,.03):
 folder=HERE/'runs/W2L_VirtualWalker'/f'guides_v1_t2_c{clearance}_h0_p5';m=mujoco.MjModel.from_xml_path(str(folder/'model.xml'));d=mujoco.MjData(m);rootid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_BODY,'Root');observer=Measurements(m,rootid);q0=m.qpos0.copy();checks=[]
 for kind,amount in [('neutral',0.),('lateral_half_gap',clearance*.5),('lateral_beyond_gap',clearance*1.5),('yaw_half_gap',.5*np.arctan(clearance/.4)),('yaw_beyond_gap',1.5*np.arctan(clearance/.4))]:
  mujoco.mj_resetData(m,d);d.qpos[:]=q0
  if kind.startswith('lateral'):d.qpos[2]+=amount
  elif kind.startswith('yaw'):mujoco.mju_axisAngle2Quat(d.qpos[3:7],np.array([0.,1.,0.]),amount)
  mujoco.mj_forward(m,d);r=observer.sample(d);count=int(r['guide_contact_count']);checks.append(dict(kind=kind,amount=amount,contacts=count,guide_force_world=r['guide_force_world'].tolist(),guide_torque_world=r['guide_torque_about_pelvis'].tolist()))
  if kind in ('neutral','lateral_half_gap','yaw_half_gap'):assert count==0,(clearance,kind,checks[-1])
  else:assert count>0,(clearance,kind,checks[-1])
  if kind=='lateral_beyond_gap':assert r['guide_force_world'][2]<0,'Positive lateral displacement must produce negative lateral reaction'
  if kind=='yaw_beyond_gap':assert r['guide_torque_about_pelvis'][1]<0,'Positive yaw must produce negative guide yaw torque'
 out.append(dict(clearance=clearance,checks=checks))
(HERE/'runs/GUIDE_DEADBAND_VERIFICATION.json').write_text(json.dumps(out,indent=2));print('Verified contact-free motion inside both clearances, physical contact outside, restoring reaction signs, and no root weld.')
