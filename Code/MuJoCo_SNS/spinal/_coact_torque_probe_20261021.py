"""Static torque-balance probe at FULL co-contraction (2026-10-21).

Set every muscle activation to 1.0 at the standing keyframe pose and
compute the net actuator torque about each leg joint (actuator_moment.T
@ actuator_force from mj_forward). This answers, without any dynamics,
WHICH joint carries a net torque at max co-contraction and in which
direction - the agonist/antagonist force-mismatch measurement behind
the co-contraction free-stand failure. Also reports the per-joint
antagonist-group force sums (PF vs DF etc.) from the same call."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np
import mujoco
import runner as R
from muscle_map import classify

model = mujoco.MjModel.from_xml_path(str(R.MODEL))
data = mujoco.MjData(model)
mujoco.mj_resetDataKeyframe(model, data, 0)
R.apply_start_pose(model, data)

acts = [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_ACTUATOR, i)
        for i in range(model.nu)]
aid = {a: i for i, a in enumerate(acts)}

data.ctrl[:] = 1.0
data.act[:] = 1.0
mujoco.mj_forward(model, data)
tau_mus = data.actuator_moment.T @ data.actuator_force   # [nv]
F = np.abs(data.actuator_force)

print("FULL-ACTIVATION STATIC TORQUE BALANCE (standing pose, no rig)")
print(f"{'joint':22s} {'net muscle torque':>18s}  dominant side")
grav = data.qfrc_bias
for jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
           "hip_flexion_l", "knee_angle_l", "ankle_angle_l",
           "lumbar_extension"):
    j = model.joint(jn)
    dof = int(model.jnt_dofadr[j.id])
    t = float(tau_mus[dof])
    print(f"{jn:22s} {t:+18.1f} N·m   (gravity bias {float(grav[dof]):+.1f})")

print("\nper-group summed |force| at full activation (N), both legs:")
groups = {}
for a in acts:
    mi = classify(a)
    if mi is None:
        continue
    groups.setdefault(mi.groups[0], 0.0)
    groups[mi.groups[0]] += float(F[aid[a]])
for g, v in sorted(groups.items(), key=lambda kv: -kv[1]):
    print(f"  {g:12s} {v:8.0f}")

print("\nankle detail R (act, group, force N, moment-arm torque N·m):")
rows = []
for a in acts:
    mi = classify(a)
    if mi is None or mi.side != "r" or a[-1] != "r":
        continue
    if mi.groups[0] in ("ankle_pf", "ankle_df"):
        dof = int(model.jnt_dofadr[model.joint("ankle_angle_r").id])
        t = float(data.actuator_moment[aid[a], dof]
                  * data.actuator_force[aid[a]])
        rows.append((mi.groups[0], a, float(F[aid[a]]), t))
pf = sum(r[3] for r in rows if r[0] == "ankle_pf")
df = sum(r[3] for r in rows if r[0] == "ankle_df")
for r in sorted(rows):
    print(f"  {r[0]:9s} {r[1]:14s} F={r[2]:7.0f}  tau={r[3]:+7.1f}")
print(f"  -> net ankle_r torque at full co-activation: PF {pf:+.1f} + "
      f"DF {df:+.1f} = {pf + df:+.1f} N·m "
      f"({'plantarflexor-dominant' if pf + df < 0 else 'dorsiflexor-dominant'})")
