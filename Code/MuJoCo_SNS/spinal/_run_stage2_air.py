"""Airwalk report (2026-09-20): run the stage-2 (afferented air) winner
through the runner with the full input schedule and report the gait-cycle
metrics from the npz with CORRECT units (q channels are already degrees).

Usage: python _run_stage2_air.py [time_s] [--gif]
"""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import _curriculum as C
import runner as R

HERE = Path(__file__).parent
t_end = "16"
args_given = [a for a in sys.argv[1:] if not a.startswith("--gif")]
if args_given:
    t_end = args_given[0]

mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st2 = json.loads((HERE / "curriculum_stage2.json")
                 .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(2, {**mul, **st2})
R.main(["--no-ground", "--time", t_end, "--drive", repr(st2["drive"])])

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(s) for s in z["neuro_names"]]
kj = [str(s) for s in z["key_joints"]]
m = (t >= 5.0) & (t <= t[-1])
print(f"\n=== stage-2 AIR window t={t[m][0]:.1f}..{t[m][-1]:.1f} s ===")
nan = not (np.all(np.isfinite(q[m])) and np.all(np.isfinite(neuro[m])))
print(f"nan: {nan}")
for j in ("knee_angle_r", "hip_flexion_r", "ankle_angle_r",
          "knee_angle_l", "hip_flexion_l"):
    i = kj.index(j)
    print(f"{j:14s} {q[m, i].min():7.1f} .. {q[m, i].max():7.1f} deg "
          f"(amp {q[m, i].max() - q[m, i].min():5.1f})")
for ch in ("RG_E_r", "RG_F_r", "RG_E_l", "RG_F_l",
           "PF_E1_r", "PF_F1_r"):
    if ch in names:
        v = neuro[m, names.index(ch)]
        print(f"{ch:14s} {v.min():6.2f} .. {v.max():6.2f} mV")
rge = neuro[m, names.index("RG_E_r")]
thr = 0.5 * max(rge.max(), 1e-9)
on = (rge > thr).astype(int)
rises = int(np.sum(np.diff(on) == 1))
dur = t[m][-1] - t[m][0]
print(f"RG_E_r bursts: {rises} in {dur:.1f} s -> {rises / dur:.2f} Hz")
rgl = neuro[m, names.index("RG_E_l")]
x = np.corrcoef(rge - rge.mean(), rgl - rgl.mean())[0, 1]
print(f"L/R RG_E correlation (want < 0): {x:.2f}")
