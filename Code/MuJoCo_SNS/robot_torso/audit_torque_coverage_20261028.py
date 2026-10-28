"""Check actual CAD tendon authority on the knee-loop tangent, not name labels."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
HERE=Path(__file__).resolve().parent
import numpy as np
import mujoco
from dynamic_routes import DynamicRoutes
from robot_kinematics import anatomical_knee_tangent
model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_abduction_diagnostic_v2.xml'));data=mujoco.MjData(model)
routes=json.loads((HERE/'robotbody_routes.json').read_text());dynamic=DynamicRoutes(model,routes)
trace=np.load(HERE/'runs/abd_defense_exact.npz')
acts=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
records=[]
for t in (0.,5.,15.,24.,35.):
    row=int(np.argmin(abs(trace['time']-t)));data.qpos[:]=trace['qpos'][row];dynamic.update(data)
    tangent=anatomical_knee_tangent(model,data)
    matrix=tangent.T@dynamic.moment_arms[model.actuator_trnid[:,0],6:].T
    records.append(np.linalg.norm(matrix,axis=0))
records=np.array(records);inactive=[a for a,norm in zip(acts,records.max(axis=0)) if norm<1e-7]
result=dict(model=str(HERE/'ah001_abduction_diagnostic_v2.xml'),sample_times_s=[0,5,15,24,35],
    actuators=acts,maximum_sampled_free_joint_moment_arm_norm_m=records.max(axis=0).tolist(),
    inactive_on_free_joint_tangent=inactive,count_inactive=len(inactive),
    scope='Five recorded poses; body internal forces can still exist, but these routes supply no resolved free-joint torque at the sampled poses. Not a whole-RoM proof.')
path=HERE.parent/'spinal/reports_20261028_family_audit/robot_torque_coverage.json';path.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
