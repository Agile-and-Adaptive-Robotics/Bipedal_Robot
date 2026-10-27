"""Mechanical checks and CAD preview; these do not constitute a walking test."""
import json
from pathlib import Path

import mujoco
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent


def main():
    model = mujoco.MjModel.from_xml_path(str(HERE/'ah001_with_torso.xml'))
    data = mujoco.MjData(model)
    mujoco.mj_forward(model, data)
    equality_rows = np.flatnonzero(data.efc_type == int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
    jacobian = data.efc_J.reshape(data.nefc,model.nv)[equality_rows]
    rank = int(np.linalg.matrix_rank(jacobian, tol=1e-8))
    assert model.neq == 2 and rank == 4
    assert model.nv-rank == 15  # free root 6, hips 6, knees 2, lumbar 1
    assert abs(model.body_mass.sum()-75) < 1e-9
    assert np.max(np.abs(data.efc_pos[equality_rows])) < 1e-8
    lumbar = mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'lumbar_flexion')
    address = model.jnt_qposadr[lumbar]
    dof_lumbar = model.jnt_dofadr[lumbar]
    analytic_moment_arms = -data.ten_J[:,dof_lumbar].copy()
    step_angle = 1e-6
    data.qpos[address] += step_angle
    mujoco.mj_forward(model,data)
    plus=data.ten_length.copy()
    data.qpos[address] -= 2*step_angle
    mujoco.mj_forward(model,data)
    minus=data.ten_length.copy()
    data.qpos[address] += step_angle
    mujoco.mj_forward(model,data)
    moment_arm_error=float(np.max(np.abs(analytic_moment_arms+(plus-minus)/(2*step_angle))))
    assert moment_arm_error < 1e-8
    torso_body=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'lumbar_torso')
    whole_com=np.sum(model.body_mass[:,None]*data.xipos,axis=0)/75
    whole_com_from_hinge=(whole_com-data.xpos[torso_body]).tolist()
    camera = mujoco.MjvCamera()
    mujoco.mjv_defaultCamera(camera)
    camera.lookat[:] = [-.06, 0, .64]
    camera.distance = 2.05
    camera.azimuth = 140
    camera.elevation = -12
    renderer = mujoco.Renderer(model, height=1000, width=1200)
    renderer.update_scene(data,camera)
    Image.fromarray(renderer.render()).save(HERE/'ah001_preview.png')
    del renderer
    back = mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_JOINT,'back_r_joint')
    dof = model.jnt_dofadr[back]
    maximum_error = 0.
    for step in range(1000):
        # Explicit diagnostic support: cancel gravity and damp every DOF. Apply
        # a small knee torque to check mobility and loop closure numerically.
        data.qfrc_applied[:] = data.qfrc_bias - .05*data.qvel
        data.qfrc_applied[dof] += .005*np.sin(2*np.pi*step*model.opt.timestep)
        mujoco.mj_step(model,data)
        rows = np.flatnonzero(data.efc_type == int(mujoco.mjtConstraint.mjCNSTR_EQUALITY))
        maximum_error = max(maximum_error,float(np.max(np.abs(data.efc_pos[rows]))))
        assert np.isfinite(data.qpos).all() and np.isfinite(data.qvel).all()
        assert not any(w.number for w in data.warning)
    assert maximum_error < .001
    output = dict(total_mass_kg=float(model.body_mass.sum()), nv=model.nv,
                  loop_constraint_rank=rank, independent_dofs=model.nv-rank,
                  supported_mechanical_test_seconds=float(data.time),
                  maximum_loop_error_m=maximum_error,
                  knee_diagnostic_angle_change_rad=float(data.qpos[model.jnt_qposadr[back]]),
                  tendon_moment_arm_finite_difference_error_m=moment_arm_error,
                  whole_model_com_from_hinge_neutral_m=whole_com_from_hinge,
                  scope='Supported mechanical exercise only; no SNS control or walking claim')
    (HERE/'ah001_verification.json').write_text(json.dumps(output,indent=2))
    print(json.dumps(output,indent=2))


if __name__ == '__main__':
    main()
