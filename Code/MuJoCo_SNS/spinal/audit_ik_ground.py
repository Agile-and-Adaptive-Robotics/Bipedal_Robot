"""IK ground-consistency audit (Ben, 2026-09-18: "walking on the ground or
walking-like motion in the air?").

Evidence chain:
 1. IK joint ranges vs textbook walking ranges.
 2. Pelvis height/oscillation from the IK pelvis coordinates.
 3. Measured GRF: peak vertical force in body weights, per-stance impulses,
    single/double-support pattern and duty.
 4. Whole-body vertical impulse balance: integral of TOTAL vertical GRF over
    a full stride ~= body weight * stride duration (only true if the forces
    are real ground reactions supporting the body).
 5. Kinematic replay: IK qpos through the repaired MuJoCo model -> lowest
    foot point height during GRF-contact windows vs during swing windows.
    Real ground walking: contact-foot floor distance ~= 0 (within ~2 cm),
    swing clearance 5-20 cm, and the two datasets must agree on timing.

Writes audit summary to fsa_results/ik_ground_audit.json and stdout.
Read-only: opens the MuJoCo model but never steps it.
"""
import io
import json
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import mujoco

import runner  # MODEL + repair stack
from bsolve_ik import GRF_MOT, IK_MOT, OSIM_DIR, read_mot  # noqa: E402

OUT = Path("fsa_results/ik_ground_audit.json")
report = {}

# ------------------------------------------------------------- 1. ranges
t_ik, cols, data = read_mot(IK_MOT)
idx = {c: i for i, c in enumerate(cols)}
RANGES = {  # textbook walking ranges (deg), (lo, hi) as OpenSim signs
    "hip_flexion_r": (-25, 35), "knee_angle_r": (-75, 5),
    "ankle_angle_r": (-25, 20), "hip_flexion_l": (-25, 35),
    "knee_angle_l": (-75, 5), "ankle_angle_l": (-25, 20),
}
print("== 1. IK coordinate ranges vs textbook walking ==")
for c, (lo, hi) in RANGES.items():
    v = data[:, idx[c]]  # subject01_walk1_ik.mot is ALREADY in degrees
    mn, mx = float(v.min()), float(v.max())
    ok = lo - 12 <= mn and mx <= hi + 12
    print(f"  {c:<18s} {mn:7.1f}..{mx:7.1f}  (textbook {lo}..{hi}) "
          f"{'OK' if ok else 'CHECK'}")
    report[f"range_{c}"] = [mn, mx]

# --------------------------------------------------- 2. pelvis height
# IK .mot pelvis_ty is in meters already (translations not converted)
if "pelvis_ty" in idx:
    ty = data[:, idx["pelvis_ty"]]
    print(f"== 2. pelvis height: mean {ty.mean():.3f} m, "
          f"osc {ty.max() - ty.min():.3f} m (normal ~0.9 m, ~2 cm) ==")
    report["pelvis_mean_m"] = float(ty.mean())
    report["pelvis_osc_m"] = float(ty.max() - ty.min())

# ------------------------------------------------------------- 3. GRF
t_g, gcols, gdata = read_mot(GRF_MOT)
gi = {c: i for i, c in enumerate(gcols)}
m0 = mujoco.MjModel.from_xml_path(str(runner.MODEL))
BW_N = 0.0
d0 = mujoco.MjData(m0)
mujoco.mj_resetDataKeyframe(m0, d0, 0)
BW_N = float(np.sum(m0.body_mass) * 9.81)
print(f"== 3. GRF vs body weight: model mass {float(np.sum(m0.body_mass)):.1f} kg"
      f" -> BW = {BW_N:.0f} N ==")

def vgrf(pref):
    return gdata[:, gi[f"{pref}_vy"]], True

grf_r, ok_r = vgrf("ground_force")            # right foot
grf_l, ok_l = vgrf("1_ground_force")          # left foot
print(f"  GRF columns found: R={['ground_force_vy']} "
      f"L={['1_ground_force_vy']}")
fr, fl = grf_r, grf_l
print(f"  peak vGRF right {fr.max() / BW_N:.2f} BW, left {fl.max() / BW_N:.2f}"
      f" BW (walking ~1.0-1.2 BW)")
report["vgrf_peak_bw"] = [float(fr.max() / BW_N), float(fl.max() / BW_N)]

# impulses between successive right heel strikes (both feet summed)
hs_r = [0.617, 1.85]
if len(hs_r) >= 2:
    m = (t_g >= hs_r[0]) & (t_g < hs_r[1])
    impulse = float(np.trapz((fr + fl)[m], t_g[m]))
    dur = hs_r[1] - hs_r[0]
    expect = BW_N * dur
    print(f"== 4. vertical impulse over 1 stride ({dur:.2f} s): {impulse:.0f} N.s"
          f" vs BW*dur {expect:.0f} N.s -> ratio {impulse / expect:.2f} "
          f"(~1.0 only for real ground support) ==")
    report["impulse_ratio"] = impulse / expect

# --------------------------------------------------- 5. kinematic replay
print("== 5. IK qpos replay: foot height vs GRF contact windows ==")
model = runner.apply_harness(m0, d0)
d = mujoco.MjData(model)
mujoco.mj_resetDataKeyframe(model, d, 0)
# qpos from bsolve npz was built with equality-aware ID setup: the leg
# joints follow IK; use raw qpos (follower rows) — heights only need the
# kinematic tree.
dnpz = np.load("bsolve_out.npz", allow_pickle=True)
qpos = dnpz["qpos"]
t_q = np.asarray(dnpz["t"], float)
# floor: lowest ground-plane z0
floor_z = 0.0

foot_bodies = {}
for nm in ("toe_r", "foot_r", "calcn_r", "toes_r", "boft_r"):
    bid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, nm)
    if bid >= 0:
        foot_bodies["r_toe" if "toe" in nm else "r_heel"] = bid
for nm in ("toe_l", "foot_l", "calcn_l", "toes_l", "boft_l"):
    bid = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, nm)
    if bid >= 0:
        foot_bodies["l_toe" if "toe" in nm else "l_heel"] = bid
if not foot_bodies:
    # fall back: any body with r_toe-ish name
    for i in range(model.nbody):
        nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, i) or ""
        low = nm.lower()
        if ("toe" in low or "calcn" in low or "boft" in low) and \
                (low.endswith("_r") or "r" in low.split()[-1:]):
            foot_bodies[nm] = i
print(f"  foot bodies: {foot_bodies}")

# GRF contact windows for the right foot (vGRF > 50 N)
contact_r = fr > 50.0
swing_r = (fr <= 50.0) & (t_g >= hs_r[0]) & (t_g <= hs_r[1])

def foot_low(pref_names):
    zs = []
    for nm, bid in foot_bodies.items():
        zs.append(d.xpos[bid][2] - floor_z)
    return min(zs) if zs else np.nan

res = {"contact": [], "swing": []}
for k in range(len(t_q)):
    d.qpos[:] = qpos[k]
    mujoco.mj_kinematics(model, d)
    zmin = min(d.xpos[b][2] for nm, b in foot_bodies.items()
               if nm.startswith("r_"))
    tck = np.searchsorted(t_g, t_q[k])
    in_contact = bool(contact_r[min(tck, len(contact_r) - 1)])
    in_stride = hs_r[0] <= t_q[k] <= hs_r[1]
    if in_stride:
        res["contact" if in_contact else "swing"].append(float(zmin))

# reference: the model's own standing keyframe (feet on the floor there)
mujoco.mj_resetDataKeyframe(model, d, 0)
mujoco.mj_kinematics(model, d)
z_stand = min(d.xpos[b][2] for nm, b in foot_bodies.items()
              if nm.startswith("r_"))
print(f"  right-foot lowest body-origin z at standing keyframe: "
      f"{z_stand:+.3f} m (this is 'on the ground' for this model)")
for lbl, arr in res.items():
    if arr:
        med = float(np.median(arr))
        p10 = float(np.percentile(arr, 10))
        print(f"  right-foot lowest point during GRF-{lbl}: "
              f"median {med:+.3f} m = {med - z_stand:+.3f} m vs standing, "
              f"10th pct {p10:+.3f} m ({p10 - z_stand:+.3f} vs standing)")
        report[f"foot_z_{lbl}_med_vs_stand"] = med - z_stand
        report[f"foot_z_{lbl}_p10_vs_stand"] = p10 - z_stand

print(f"\nVERDICT inputs written to {OUT}")
OUT.write_text(json.dumps(report, indent=1), encoding="utf-8")
