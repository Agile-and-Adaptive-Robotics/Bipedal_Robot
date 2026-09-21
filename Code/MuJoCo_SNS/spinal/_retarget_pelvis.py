"""Pelvis re-target (Ben's "run our own IK on gait2392", 2026-09-18).

The scaled subject's IK joint angles are valid on any adult skeleton,
but the subject's legs are ~7-8% longer than stock gait2392
(audit_ik_ground.py: standing pelvis 1.02 m vs 0.95 m), so replaying
the IK verbatim leaves the loaded foot ~8.7 cm above our MJCF floor.
Keep the measured JOINT angles; re-solve only the pelvis HEIGHT per
frame so the GRF-loaded foot sits on the floor.

Closed form: qpos[pelvis_ty] is a world-vertical slide, so lifting the
pelvis by d lifts both feet by exactly d (verified numerically per
run).  With per-foot contact weights w from the measured vertical GRF,
d = (w_r*(z0 - z_r) + w_l*(z0 - z_l)) / (w_r + w_l), where z0 is the
foot-origin height of the model's own standing keyframe ("on the
floor" for this model).  Frames with no meaningful load hold the
previous d (regularization; walking always has one loaded foot).

Writes bsolve_retarget.npz = bsolve_out.npz with corrected qpos plus
'retarget_dty', and prints before/after contact-window foot heights.
Pure post-process: read-only on the model, no bsolve rerun needed.
"""
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import mujoco

import bsolve_ik as B
import runner

SRC = Path("bsolve_out.npz")
DST = Path("bsolve_retarget.npz")

d = np.load(SRC, allow_pickle=True)
qpos = np.array(d["qpos"], float)
t = np.asarray(d["t"], float)
T, nq = qpos.shape

model, data = B.build_air_model()          # kinematics only
adr_ty = model.jnt_qposadr[
    mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_JOINT, "pelvis_ty")]

foot_ids = [mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, nm)
            for nm in ("toes_r", "calcn_r", "toes_l", "calcn_l")]
foot_ids = [i for i in foot_ids if i >= 0]
assert len(foot_ids) == 4, "foot bodies missing"

def foot_z_low(row=None):
    if row is not None:
        data.qpos[:] = row
    mujoco.mj_kinematics(model, data)
    return min(float(data.xpos[b][2]) for b in foot_ids)

def foot_z_sides():
    mujoco.mj_kinematics(model, data)
    zr = min(float(data.xpos[b][2]) for b in foot_ids[:2])
    zl = min(float(data.xpos[b][2]) for b in foot_ids[2:])
    return zr, zl

# --- sanity: dz/dty must be exactly 1 (world-vertical slide) ----------
data.qpos[:] = qpos[0]
z0r, z0l = foot_z_sides()
data.qpos[adr_ty] += 0.01
z1r, z1l = foot_z_sides()
print(f"slide probe: dz_r {z1r - z0r:+.6f}, dz_l {z1l - z0l:+.6f} "
      f"(expect +0.01)")
assert abs((z1r - z0r) - 0.01) < 1e-6 and abs((z1l - z0l) - 0.01) < 1e-6, \
    "pelvis_ty is not a world-vertical slide on this model"
print("slide check: dz/dty = 1 exactly for both feet")

# --- floor reference = the model's own standing keyframe --------------
mujoco.mj_resetDataKeyframe(model, data, 0)
mujoco.mj_kinematics(model, data)
z_stand = min(float(data.xpos[b][2]) for b in foot_ids)
print(f"standing-keyframe foot-origin z (floor ref): {z_stand:+.4f} m")

# --- per-foot contact weights from the measured vertical GRF ----------
t_g, g_names, g_vals = B.read_mot(B.GRF_MOT)
Fr, Pr, Tr = B.grf_at(t, t_g, g_names, g_vals, pref="")
Fl, Pl, Tl = B.grf_at(t, t_g, g_names, g_vals, pref="1_")
Fvr, Fvl = Fr[:, 2], Fl[:, 2]
w_r = np.clip(Fvr / 150.0, 0.0, 1.0)
w_l = np.clip(Fvl / 150.0, 0.0, 1.0)

# --- closed-form correction -------------------------------------------
qpos_new = qpos.copy()
dty = np.zeros(T)
prev = 0.0
for k in range(T):
    data.qpos[:] = qpos[k]
    zr, zl = foot_z_sides()
    ws = w_r[k] + w_l[k]
    if ws < 0.05:
        dd = prev
    else:
        dd = (w_r[k] * (z_stand - zr) + w_l[k] * (z_stand - zl)) / ws
    dty[k] = dd
    prev = dd
    qpos_new[k, adr_ty] = qpos[k, adr_ty] + dd

print(f"dty: mean {dty.mean():+.4f} m  range [{dty.min():+.4f}, "
      f"{dty.max():+.4f}]  (subject legs longer than stock -> negative)")

# --- validation: contact-window foot height before/after ---------------
contact_r = Fvr > 50.0
in_stride = (t >= t[0]) & (t <= t[-1])
def window_stats(rows):
    zc, zs = [], []
    for k in range(T):
        data.qpos[:] = rows[k]
        zr, zl = foot_z_sides()
        if in_stride[k]:
            (zc if contact_r[k] else zs).append(
                min(zr, zl) if contact_r[k] else zr)
    return (float(np.median(zc)) if zc else np.nan,
            float(np.percentile(zc, 10)) if zc else np.nan,
            float(np.median(zs)) if zs else np.nan)

bc, bc10, bs = window_stats(qpos)
ac, ac10, asw = window_stats(qpos_new)
print(f"right-foot z during GRF-contact: BEFORE median {bc:+.3f} "
      f"(p10 {bc10:+.3f})  ->  AFTER {ac:+.3f} (p10 {ac10:+.3f}); "
      f"floor ref {z_stand:+.3f}")
print(f"right-foot z during swing:       BEFORE median {bs:+.3f}  ->  "
      f"AFTER {asw:+.3f} (clearance above floor ref "
      f"{asw - z_stand:+.3f} m)")

ok = abs(ac - z_stand) < 0.02 and abs(ac10 - z_stand) < 0.03
print("VERDICT:", "PASS - loaded foot now sits on the model floor"
      if ok else "CHECK - residual contact error remains")

out = {k: d[k] for k in d.files}
out["qpos"] = qpos_new
out["retarget_dty"] = dty
out["retarget_note"] = np.array(
    "2026-09-18 pelvis re-target: measured joint angles kept, pelvis "
    "height re-solved on gait2392 so GRF-loaded foot sits on the model "
    "floor (closed-form dty, vGRF-weighted; see _retarget_pelvis.py)")
np.savez(DST, **out)
print(f"written: {DST}")
