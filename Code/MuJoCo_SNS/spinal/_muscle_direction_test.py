"""Definitive muscle-direction test: activate ONE knee muscle at a time in
the real simulation (real actuators, not moment matrices) and record where
the knee_angle goes and which way the heel travels.

Geometry ground truth (from _knee_sign_test.py): model faces +X;
NEGATIVE knee_angle = human flexion (heel up & backward);
POSITIVE knee_angle = impossible forward bend.
"""
import mujoco
import numpy as np
import runner as R

model = mujoco.MjModel.from_xml_path(str(R.MODEL))
data = mujoco.MjData(model)
mujoco.mj_resetDataKeyframe(model, data, 0)
kp = R.capture_pose(model)

TESTS = ["vas_lat_r", "rect_fem_r", "semimem_r", "bifemsh_r",
         "med_gas_r", "none"]
for muscle in TESTS:
    m = R.apply_harness(model, data, kxy=1500.0, leg_damping=2.0)
    d = mujoco.MjData(m)
    R.seed_pose(m, d, kp)
    mujoco.mj_forward(m, d)
    aid = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_ACTUATOR, muscle) \
        if muscle != "none" else -1
    kadr = int(m.joint("knee_angle_r").qposadr[0])
    knee_hist = []
    for k in range(500):  # 1 s at dt=2 ms
        if aid >= 0:
            d.ctrl[aid] = 1.0
        d.qvel[:] *= 0.98          # light global damping for readability
        if k < 100:
            d.qvel[:] = 0.0        # settle at keyframe first
            mujoco.mj_forward(m, d)
        else:
            mujoco.mj_step(m, d)
        knee_hist.append(np.degrees(d.qpos[kadr]))
    kh = np.array(knee_hist[100:])
    print(f"{muscle:12s} knee: start {kh[0]:+7.1f}  min {kh.min():+7.1f}  "
          f"max {kh.max():+7.1f}  end {kh[-1]:+7.1f}  "
          f"-> {'FLEXION (-)' if kh.min() < kh[0] - 1 else 'EXTENSION/other (+)' if kh.max() > kh[0] + 1 else 'no motion'}")
