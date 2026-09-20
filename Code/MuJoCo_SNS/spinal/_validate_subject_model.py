"""Post-conversion validation gate for a new subject01 MJCF.

Usage:
    set AARL_MODEL=<path to new xml>   (or pass the path as argv[1])
    python _validate_subject_model.py

Checks (GO/NO-GO before any tuning work touches the new model):
 1. loads in mujoco 2.3.7
 2. 92 actuators, and the key muscle names bsolve/runner expect
 3. standing keyframe: pelvis height ~1.02 m (the SCALED subject's
    standing height; the old unscaled model sits at 0.95 m)
 4. lowest foot body-origin z at the keyframe ~on the floor (old model:
    +0.014 m); leg-drop pelvis_z - foot_z recorded for scale comparison
 5. 1000 steps at 2 ms with zero ctrl: finite, no explosion
Prints PASS/FAIL per check + a one-line verdict; exits nonzero on fail.
"""
import io
import os
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import mujoco

model_path = Path(sys.argv[1] if len(sys.argv) > 1
                  else os.environ.get("AARL_MODEL", ""))
if not model_path.is_file():
    print("FAIL: no model path (argv[1] or AARL_MODEL)")
    sys.exit(2)

m = mujoco.MjModel.from_xml_path(str(model_path))
d = mujoco.MjData(m)
mujoco.mj_resetDataKeyframe(m, d, 0) if m.nkey else mujoco.mj_resetData(m, d)
mujoco.mj_forward(m, d)

fails = []

# 2. actuators
acts = [mujoco.mj_id2name(m, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(m.nu)]
print(f"actuators: {m.nu} (expect 92)")
if m.nu != 92:
    fails.append("actuator count")
need = ("vas_lat_r", "soleus_r", "rect_fem_r", "semimem_r", "tib_ant_l")
missing = [n for n in need if n not in acts]
if missing:
    fails.append(f"missing muscles {missing}")

def body_z(nm):
    b = mujoco.mj_name2id(m, mujoco.mjtObj.mjOBJ_BODY, nm)
    return float(d.xpos[b][2]) if b >= 0 else None

# 3+4. standing geometry (uses the keyframe pose as-is)
pel = body_z("pelvis")
foot_zs = {nm: body_z(nm) for nm in
           ("toe_r", "calcn_r", "toe_l", "calcn_l") if body_z(nm) is not None}
foot_low = min(foot_zs.values()) if foot_zs else None
print(f"standing pelvis z: {pel:.3f} m (scaled subject ~1.02; "
      f"old unscaled model 0.95)")
if pel is None or not 0.98 <= pel <= 1.10:
    fails.append("pelvis standing height out of scaled-subject range")
if foot_low is not None:
    print(f"lowest foot origin z at keyframe: {foot_low:+.3f} m "
          f"(old model +0.014)")
    if not -0.02 <= foot_low <= 0.10:
        fails.append("foot origin z at keyframe out of range")
    print(f"leg drop (pelvis - foot): {pel - foot_low:.3f} m "
          f"(old model 0.936)")

# 5. zero-ctrl stability
bad = False
for _ in range(1000):
    d.ctrl[:] = 0.0
    mujoco.mj_step(m, d)
    if not np.isfinite(d.qpos).all() or not np.isfinite(d.qvel).all():
        bad = True
        break
print(f"1000-step zero-ctrl: {'EXPLODED' if bad else 'finite'} "
      f"(falling over is fine; NaN is not)")
if bad:
    fails.append("nonfinite dynamics at zero ctrl")

if fails:
    print("VERDICT: FAIL -", "; ".join(fails))
    sys.exit(1)
print("VERDICT: PASS - model is loadable, complete, scaled, and stable "
      "enough for the harness/integration audit")
print("NEXT: point AARL_MODEL at it and re-run audit_ik_ground.py - the "
      "GRF-contact foot height should move from +0.087 m to ~0.")
