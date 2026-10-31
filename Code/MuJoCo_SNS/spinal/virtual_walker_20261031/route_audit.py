import json
import numpy as np
from import_virtual_walker import build,mujoco,safe
def audit(newer=False):
 path=build(newer);m=mujoco.MjModel.from_xml_path(str(path));d=mujoco.MjData(m);mujoco.mj_forward(m,d);meta=json.loads((path.parent/'IMPORT_AUDIT.json').read_text());jn=['hip_L','knee_L','ankle_L','hip_R','knee_R','ankle_R'];rows=[]
 for bank in meta['muscles']:
  if bank['type']!='LinearHillMuscle':continue
  tid=mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_TENDON,'route_'+safe(bank['name']));arms={}
  for joint in jn:
   q=m.jnt_qposadr[mujoco.mj_name2id(m,mujoco.mjtObj.mjOBJ_JOINT,joint)];d.qpos[q]=1e-5;mujoco.mj_forward(m,d);plus=d.ten_length[tid];d.qpos[q]=-1e-5;mujoco.mj_forward(m,d);minus=d.ten_length[tid];d.qpos[q]=0;arms[joint]=float(-(plus-minus)/2e-5)
  rows.append(dict(muscle=bank['name'],moment_arms_m=arms,spanned_joints=[j for j,v in arms.items() if abs(v)>1e-4]))
 out=dict(method='Central tendon-length finite difference at assembled pose; source physical routes, not muscle-name inference.',muscles=rows,all_uniarticular=all(len(x['spanned_joints'])==1 for x in rows));(path.parent/'ROUTE_AUDIT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
if __name__=='__main__':audit()
