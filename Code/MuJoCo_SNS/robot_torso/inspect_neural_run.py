import json,sys
from pathlib import Path
import numpy as np
import mujoco
from dynamic_routes import DynamicRoutes
from run_neural_robot import rotation_error
HERE=Path(__file__).resolve().parent
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'))
config=json.loads((HERE/'robotbody_routes.json').read_text());dynamic=DynamicRoutes(model,config)
log=np.load(sys.argv[1])
for index in (0,min(10,len(log['time'])-1),-1):
 data=mujoco.MjData(model);data.qpos[:]=log['qpos'][index];data.qvel[:]=log['qvel'][index];data.act[:]=log['activation'][index]
 dynamic.update(data)
 force=-data.actuator_force;predicted=dynamic.moment_arms[model.actuator_trnid[:,0]].T@force
 actual=data.qfrc_actuator+dynamic.extra_length_jacobian[model.actuator_trnid[:,0]].T@data.actuator_force
 print('T',log['time'][index],'ROOT',data.qpos[:3],'TILT',np.rad2deg(rotation_error(data.qpos[3:7])), 'KNEE',log['knee_coordinate'][index],
       'force consistency',np.max(np.abs(actual-predicted)), 'contacts',data.ncon)
 for side in ('r','l'):
  body=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+side)
  print('FOOT',side,data.xpos[body], 'HIP',np.rad2deg(rotation_error(data.qpos[model.jnt_qposadr[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint')]:][:4])))
 print('DESIRED',np.round(log['desired_activation'][index],2),'ACTUAL',np.round(data.act,2))
