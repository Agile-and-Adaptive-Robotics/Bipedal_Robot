"""Identify actual abducting leverage among the selected robot actuators."""
import json
from run_neural_robot import *
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'));data=mujoco.MjData(model)
dynamic=DynamicRoutes(model,json.loads((HERE/'robotbody_routes.json').read_text()))
dynamic.update(data)
rows=[]
for side,sign in [('r',-1.),('l',1.)]:
    joint=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'thigh_'+side+'_joint')
    dof=model.jnt_dofadr[joint]
    arms=sign*dynamic.moment_arms[model.actuator_trnid[:,0],dof]
    capacities=-model.actuator_gainprm[:,0]*arms
    entries=[dict(actuator=mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i),
                  abduction_arm_m=float(arms[i]),capacity_Nm=float(capacities[i]))
             for i in range(model.nu) if capacities[i]>.05]
    entries.sort(key=lambda x:x['capacity_Nm'],reverse=True)
    rows.append(dict(side=side,abduction_capacity_Nm=float(capacities[capacities>0].sum()),
                     adduction_capacity_Nm=float(-capacities[capacities<0].sum()),actuators=entries))
(HERE/'abduction_capacity_probe.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows,indent=2))
