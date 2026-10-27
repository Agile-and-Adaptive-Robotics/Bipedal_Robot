"""Verify the corrected robot frame, joints, route map, and reference capacities."""
import json
from pathlib import Path
import mujoco
import numpy as np
from dynamic_routes import DynamicRoutes
HERE=Path(__file__).resolve().parent
def main():
    model=mujoco.MjModel.from_xml_path(str(HERE/'ah001_neural.xml'))
    config=json.loads((HERE/'robotbody_routes.json').read_text())
    dynamic=DynamicRoutes(model,config);data=mujoco.MjData(model);dynamic.update(data)
    assert np.allclose(model.opt.gravity,[0,-9.81,0])
    ground=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_GEOM,'ground')
    assert np.allclose(data.geom_xmat[ground].reshape(3,3)[:,2],[0,1,0],atol=1e-12)
    lumbar=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'lumbar_ball')
    assert model.jnt_type[lumbar]==int(mujoco.mjtJoint.mjJNT_BALL)
    assert abs(model.body_mass.sum()-75)<1e-9
    assert model.nu==54 and model.ntendon==54
    for side in ('r','l'):
        names=[mujoco.mj_id2name(model,mujoco.mjtObj.mjOBJ_ACTUATOR,i) for i in range(model.nu)]
        assert sum(name.endswith('_'+side) for name in names)==27
        foot=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'foot_'+side)
        assert model.body_jntnum[foot]==0
        es=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_ACTUATOR,'ercspn_'+side)
        assert -model.actuator_gainprm[es,0]==2500
    rows=np.flatnonzero(data.efc_type==int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
    jac=data.efc_J.reshape(data.nefc,model.nv)[rows]
    rank=int(np.linalg.matrix_rank(jac,tol=1e-8))
    assert rank==4 and model.nv-rank==17
    assert np.max(np.abs(data.efc_pos[rows]))<1e-8
    assert data.ncon>0
    data.act[:]=1
    dynamic.update(data)
    assert np.allclose(data.actuator_force,model.actuator_gainprm[:,0])
    moving_extra=float(np.max(np.abs(dynamic.extra_length_jacobian)))
    assert moving_extra>1e-4
    output=dict(frame='x anterior, y superior, z right',mass_kg=float(model.body_mass.sum()),
        actuators=model.nu,tendons=model.ntendon,actuators_per_side=27,raw_velocity_dofs=model.nv,
        knee_loop_constraint_rank=rank,independent_dofs=model.nv-rank,lumbar_joint='ball',
        initial_ground_contacts=data.ncon,maximum_moving_guide_jacobian_term_m=moving_extra,
        status='Frame, construction, reference capacity, contact, and dynamic-path checks pass; control performance is assessed separately')
    (HERE/'ah001_verification.json').write_text(json.dumps(output,indent=2))
    print(json.dumps(output,indent=2))
if __name__=='__main__':main()
