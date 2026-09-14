"""Probe: pelvis translation mapping - foot heights/sides for tz sign."""
import numpy as np
import mujoco

import bsolve_ik as B

t_ik, ik_names, ik_vals = B.read_mot(B.IK_MOT)
ik = {n: ik_vals[:, i] for i, n in enumerate(ik_names)}
model, data = B.build_air_model()
tal_r = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "talus_r")
tal_l = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "talus_l")

for tz_sign in (1.0, -1.0):
    print(f"--- pelvis_tz_os * {tz_sign:+.0f} ---")
    for k in (0, 30, 60, 90, 120):
        mujoco.mj_resetDataKeyframe(model, data, 0)
        for name in B.ROT_DRIVERS:
            jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            data.qpos[model.jnt_qposadr[jid]] = \
                np.deg2rad(ik[name][k]) * 1.0
        for name, val in (("pelvis_tx", ik["pelvis_tx"][k]),
                          ("pelvis_ty", ik["pelvis_ty"][k]),
                          ("pelvis_tz", tz_sign * ik["pelvis_tz"][k])):
            jid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, name)
            data.qpos[model.jnt_qposadr[jid]] = val
        B.apply_eq_followers(model, data)
        mujoco.mj_forward(model, data)
        print(f"  t={t_ik[k]:.2f} talus_r z {data.xpos[tal_r][2]:+0.3f} "
              f"y {data.xpos[tal_r][1]:+0.3f} | talus_l z "
              f"{data.xpos[tal_l][2]:+0.3f} y {data.xpos[tal_l][1]:+0.3f} | "
              f"pelvis os ty {ik['pelvis_ty'][k]:.2f} tz "
              f"{ik['pelvis_tz'][k]:+.3f}")
