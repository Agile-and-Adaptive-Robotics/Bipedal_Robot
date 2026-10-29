"""Pelvis-frame ankle/knee ordering from recorded MuJoCo states.

Two centimetres of body-origin separation is an engineering screening gate,
not a full foot-edge clearance or self-collision test and not a literature value.
"""
import numpy as np
import mujoco

def assess_frontal(path):
    import bsolve_ik as B
    model, data = B.build_air_model(False)
    ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, name)
           for name in ('pelvis', 'calcn_r', 'calcn_l', 'tibia_r', 'tibia_l')]
    mujoco.mj_kinematics(model, data)
    hips = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            for name in ('hip_flexion_r', 'hip_flexion_l')]
    separation = data.xanchor[hips[0]] - data.xanchor[hips[1]]
    right_local = data.xmat[ids[0]].reshape(3, 3).T @ (separation / np.linalg.norm(separation))
    # This converted model preserves OpenSim's local z-right pelvis frame.
    # Its 90-degree body rotation maps local +z to approximately world -y.
    # Applying that rotation to a world -y vector would rotate it twice.
    if not np.allclose(right_local, [0., 0., 1.], atol=1e-3):
        raise ValueError('Hip anchors do not verify the expected OpenSim pelvis frame.')
    with np.load(path, allow_pickle=True) as z:
        q = z['qfull']
        keep = np.flatnonzero((z['t'] >= 5.) & (z['t'] < 13.))[::10]
        if not len(keep) or not np.all(np.isfinite(q[keep])):
            return dict(frontal_gate=False, frontal_incomplete=True)
        ankle = []; knee = []
        for i in keep:
            data.qpos[:] = q[i]
            mujoco.mj_kinematics(model, data)
            right = data.xmat[ids[0]].reshape(3, 3) @ right_local
            ankle.append(float((data.xpos[ids[1]] - data.xpos[ids[2]]) @ right))
            knee.append(float((data.xpos[ids[3]] - data.xpos[ids[4]]) @ right))
    return dict(frontal_gate=bool(min(ankle) > .02 and min(knee) > .02),
        minimum_ankle_width_m=min(ankle), minimum_knee_width_m=min(knee),
        crossed_ankle_fraction=float(np.mean(np.asarray(ankle) < 0.)),
        pelvis_local_right_axis=right_local.tolist(),
        anatomical_axis_verified_from_hip_anchors=True,
        frontal_scope='Pelvis-frame calcaneus/tibia origins; not foot-edge clearance or full self-collision.')
