"""Static hold capacity probe (2026-09-30): can the M1 body's muscles
support body weight AT ALL? Full activation (ctrl=1 on all 12 muscles)
from the spawn keyframe, 3 s, harness releases after 1 s. If the pelvis
collapses even at max effort, the rigid-tendon/no-Kse-Kpe port lacks the
passive load path Li's LinearHill muscles provide, and ground walking
needs the elastic muscle model - not tuning.
Usage: python probe_hold.py [hold_s=1.0]
"""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import numpy as np
import mujoco

HERE = os.path.dirname(os.path.abspath(__file__))
MJCF = os.path.join(HERE, "w2l_mjcf_fixed.xml")
HOLD = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
MODE = (sys.argv[2] if len(sys.argv) > 2 else "all").lower()  # all|ext

m = mujoco.MjModel.from_xml_path(MJCF)
d = mujoco.MjData(m)
mujoco.mj_resetDataKeyframe(m, d, 0)
mujoco.mj_forward(m, d)

_ctrl = np.zeros(m.nu)
for i in range(m.nu):
    nm = mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
    if MODE == "all" or "_ext" in nm:
        _ctrl[i] = 1.0
root_bid = int(m.jnt_bodyid[0])
mass = float(np.sum(m.body_mass))
wgt = mass * (-m.opt.gravity[2])
z0 = float(d.qpos[2])
K = 4.0 * wgt / 0.01
C = 2.0 * float(np.sqrt(K * mass))
print(f"mass {mass:.2f} kg  weight {wgt:.1f} N  spawn z {z0:.3f} m  "
      f"timestep {m.opt.timestep}")
print("muscles:", [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
                   for i in range(m.nu)])
print("sum Fmax:", float(np.sum(m.actuator_gainprm[:, 2])), "N")

n = int(3.0 / m.opt.timestep)
zmin_after = 1e9
for i in range(n):
    t = i * m.opt.timestep
    d.ctrl[:] = _ctrl                   # probe activation pattern
    if t < HOLD:
        Fz = wgt + K * (z0 - d.qpos[2]) - C * d.qvel[2]
        d.xfrc_applied[root_bid, 2] = max(Fz, 0.0)
    else:
        d.xfrc_applied[root_bid, :] = 0.0
        zmin_after = min(zmin_after, float(d.qpos[2]))
    mujoco.mj_step(m, d)
    if not np.all(np.isfinite(d.qpos)):
        print(f"NaN at t={t:.2f}")
        break
    if i % 500 == 0:
        kn = [float(np.degrees(d.qpos[m.jnt_qposadr[j]]))
              for j in range(m.njnt) if m.jnt_type[j] == 3]
        print(f"t={t:4.1f}  z={d.qpos[2]:.3f}  hinges(deg)="
              + " ".join(f"{v:+6.1f}" for v in kn))
print(f"\nmin pelvis z AFTER release: {zmin_after:.3f} m "
      f"(spawn {z0:.3f}) -> "
      f"{'HOLDS at max effort' if zmin_after > 0.8 * z0 else 'COLLAPSES even at max effort'}")
