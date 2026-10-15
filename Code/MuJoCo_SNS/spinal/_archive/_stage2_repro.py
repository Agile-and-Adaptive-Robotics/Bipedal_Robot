"""Exact reproduction of the curriculum stage-2 objective for the stored
winner: same flags (--no-ground --time 14 --drive repr), same window
(t 5..17), same metric (3*rises + 0.5*(-knee_min), RG_E_r = neuro[:,2]).
Usage: python _stage2_repro.py
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
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st2 = json.loads((HERE / "curriculum_stage2.json")
                 .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(2, {**mul, **st2})
R.main(["--no-ground", "--time", "14", "--drive", repr(st2["drive"])])

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
m = (t >= 5.0) & (t <= 17.0)
print("finite:", bool(np.all(np.isfinite(q[m])) and np.all(np.isfinite(neuro[m]))))
rge = neuro[m, 2]
thr = 0.5 * max(rge.max(), 1e-9)
on = (rge > thr).astype(int)
rises = int(np.sum(np.diff(on) == 1))
knee = q[m, 4]
score = 3.0 * rises + 0.5 * (-float(knee.min()))
print(f"RG_E_r window max {rge.max():.2f} mV, rises {rises}, "
      f"knee_min {float(knee.min()):.1f} deg -> objective {score:.3f} "
      f"(study recorded 36.906)")
