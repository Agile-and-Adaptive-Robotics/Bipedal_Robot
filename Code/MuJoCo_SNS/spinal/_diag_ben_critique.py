"""Both-leg diagnosis of the s3b winner run (Ben's critique 2026-09-20):
left vs right activity, hip DC offset, ankle pattern, and - via a
re-run with contact logging - per-foot ground contact and pelvis height
vs the standing keyframe. Ground truth for what to fix."""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import _curriculum as C
import runner as R
import params as P

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(3, {**mul, **st3})

# instrument: patch the stance-feedback block's inputs by running the
# runner normally but capturing contact via monkeypatched mujoco calls
# is heavy; instead run and read npz + qfull foot heights afterward.
R.main(["--eval", "--drive", repr(st3["drive"])])

import mujoco

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro, qf = z["t"], z["q"], z["neuro"], z["qfull"]
names = [str(s) for s in z["neuro_names"]]
kj = [str(s) for s in z["key_joints"]]
m = (t >= 5.0)

print("=== joint ranges t>=5 s, RIGHT vs LEFT (deg, OpenSim sign) ===")
for j in ("hip_flexion_r", "hip_flexion_l", "knee_angle_r",
          "knee_angle_l", "ankle_angle_r", "ankle_angle_l"):
    i = kj.index(j)
    v = q[m, i]
    print(f"{j:15s} {v.min():7.1f}..{v.max():6.1f}  mean {v.mean():6.1f} "
          f" amp {v.max() - v.min():5.1f}")

print("\n=== neural, RIGHT vs LEFT (mV) ===")
for ch in ("RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
           "PF_E1_l", "PF_F1_l"):
    if ch in names:
        v = neuro[m, names.index(ch)]
        print(f"{ch:9s} {v.min():6.2f}..{v.max():6.2f}  swing "
              f"{v.max() - v.min():5.2f}")

# interleg antiphase quality
e_r = neuro[m, names.index("RG_E_r")]
e_l = neuro[m, names.index("RG_E_l")]
c = np.corrcoef(e_r - e_r.mean(), e_l - e_l.mean())[0, 1]
print(f"\nRG_E R/L corr (antiphase wants < 0): {c:+.2f}")

# foot heights + contact from the model, replaying qfull
model = mujoco.MjModel.from_xml_path(str(R.MODEL))
data = mujoco.MjData(model)
feet = {}
for b in range(model.nbody):
    nm = mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_BODY, b) or ""
    for side in ("r", "l"):
        if nm.endswith("_" + side) and any(
                f in nm for f in ("calcn", "toes")):
            feet.setdefault(side, []).append(b)
print("\nfoot bodies:", {s: [mujoco.mj_id2name(
    model, mujoco.mjtObj.mjOBJ_BODY, b) for b in v] for s, v in
    feet.items()})

nrow = int(m.sum())
ht = {s: np.zeros(nrow) for s in feet}
load = {s: np.zeros(nrow) for s in feet}
ti = np.flatnonzero(m)
for k, step in enumerate(ti):
    data.qpos[:] = qf[step, :model.nq]
    data.qvel[:] = 0.0
    mujoco.mj_forward(model, data)
    for s, bodies in feet.items():
        ht[s][k] = min(data.xpos[b, 2] for b in bodies)
        fsum = 0.0
        for ci in range(data.ncon):
            con = data.contact[ci]
            for g in (con.geom1, con.geom2):
                if model.geom_bodyid[g] in bodies:
                    mujoco.mj_contactForce(model, data, ci,
                                           np.zeros(6))
                    fsum += 1
                    break
        load[s][k] = fsum
for s in ("r", "l"):
    frac = float((load[s] > 0).mean())
    print(f"foot_{s}: height min {ht[s].min():.3f} mean {ht[s].mean():.3f}"
          f" max {ht[s].max():.3f} m | contact frames {frac * 100:.0f}%")
