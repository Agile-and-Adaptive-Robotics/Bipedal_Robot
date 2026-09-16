"""Ground truth for bridge test c (and reference data for tests a/b).

Runs the sensor-patched simbridge model in Python mujoco 2.3.7 with the
SAME initial state as the patched blockset (keyframe 0) and a CONSTANT
ctrl = 1.0 on vas_med_r (actuator index 28 of 92), stepping with mj_step
for 1 s (200 x 5 ms). Saves sensordata per step to logs/gt_const.npz.
"""
from pathlib import Path

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
MODEL = (HERE.parents[4] / "Solid_Models" / "OpenSim" / "Gait2392_Robotbody"
         / "mjc" / "gait2392_simbody" / "gait2392_simbody_cvt3_simbridge.xml")

VAS_MED_R = 28          # actuator index, of 92
NSTEP = 80              # 0.4 s at 5 ms: the uncontrolled model (only one
                        # muscle fired) goes solver-unstable at ~0.445 s,
                        # so stay inside the clean window for comparisons

m = mujoco.MjModel.from_xml_path(str(MODEL))
d = mujoco.MjData(m)
mujoco.mj_resetDataKeyframe(m, d, 0)      # match the patched block init
d.ctrl[:] = 0.0
d.ctrl[VAS_MED_R] = 1.0

sens = np.zeros((NSTEP, m.nsensordata))
knee = np.zeros(NSTEP)
kadr = m.jnt_qposadr[mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_JOINT, "knee_angle_r")]
for k in range(NSTEP):
    mujoco.mj_step(m, d)
    sens[k] = d.sensordata
    knee[k] = d.qpos[kadr]
    if not np.isfinite(sens[k]).all():
        raise SystemExit(f"non-finite sensordata at step {k}")

np.savez(HERE.parent / "logs" / "gt_const.npz", sens=sens, knee=knee,
         names=["knee_r_pos", "knee_r_vel", "vas_med_r_len",
                "vas_med_r_vel", "vas_med_r_frc"])
from scipy.io import savemat  # local import: scipy only needed for the .mat
savemat(HERE.parent / "logs" / "gt_const.mat", {"sens": sens, "knee": knee})
print(f"ran {NSTEP} steps clean")
print(f"knee_r_pos: start {np.degrees(sens[0, 0]):+.2f} deg  "
      f"end {np.degrees(sens[-1, 0]):+.2f} deg")
print(f"vas_med_r_frc: first {sens[0, 4]:+.1f} N  "
      f"max |frc| {np.abs(sens[:, 4]).max():+.1f} N  end {sens[-1, 4]:+.1f} N")
