"""Measure foot heading and lateral ordering in the recorded pelvis frame."""
import json
from pathlib import Path
from run_neural_robot import rotation_error
import mujoco
import numpy as np
HERE=Path(__file__).resolve().parent
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'));data=mujoco.MjData(model)
trace=np.load(HERE/'runs/verified65_exact_fast.npz')
pelvis=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'pelvis')
feet=[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+s) for s in ('r','l')]
headings=[];lateral=[]
for q in trace['qpos'][trace['time']>=15][::5]:
    data.qpos[:]=q;mujoco.mj_kinematics(model,data)
    rotation=data.xmat[pelvis].reshape(3,3).T
    directions=[rotation@data.xmat[f].reshape(3,3)[:,0] for f in feet]
    headings.append([np.rad2deg(np.arctan2(v[2],v[0])) for v in directions])
    lateral.append([(rotation@(data.xpos[f]-data.xpos[pelvis]))[2] for f in feet])
headings=np.array(headings);lateral=np.array(lateral)
result=dict(run='verified65_exact_fast',window_start_s=15.,
    heading_definition='atan2(foot-local anterior axis right-component, anterior-component), expressed in pelvis ISB frame',
    minimum_foot_heading_deg=headings.min(axis=0).tolist(),maximum_foot_heading_deg=headings.max(axis=0).tolist(),
    median_foot_heading_deg=np.median(headings,axis=0).tolist(),
    minimum_right_minus_left_ankle_lateral_spacing_m=float(np.min(lateral[:,0]-lateral[:,1])),
    ankle_order_reversed_fraction=float(np.mean(lateral[:,0]<lateral[:,1])),
    scope='Body-axis heading and ankle ordering only; foot meshes can overlap before ankle order reverses')
(HERE/'foot_heading_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
