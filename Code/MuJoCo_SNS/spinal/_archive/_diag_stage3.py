"""Diagnose the stage-3 ground failure: run the seeded config in eval
mode and capture the metrics + per-0.5s forensics."""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import _curriculum as C
import runner as R

p = json.loads(open("curriculum_stage3.json", encoding="utf-8").read())["params"]
mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(3, {**mul, **p})
m = R.main(["--eval", "--drive", repr(p["drive"])])
print("\n=== eval metrics ===")
for k, v in m.items():
    print(f"  {k}: {v}")
z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(s) for s in z["neuro_names"]]
w = (t >= 5.0) & (t <= t[-1])
if w.any():
    for ch in ("RG_E_r", "RG_F_r", "PF_E1_r", "PF_F1_r"):
        if ch in names:
            v = neuro[w, names.index(ch)]
            print(f"{ch}: {v.min():.2f}..{v.max():.2f} mV")
    # q channels are ALREADY degrees of the named KEY_JOINTS (runner fix
    # 2026-09-18) - no np.degrees here; that double conversion printed
    # knee "-3250 deg" in the 09-18 log.
    knee = q[w, 4]
    hip = q[w, 3]
    print(f"knee_r: {knee.min():.1f}..{knee.max():.1f} deg, "
          f"hip_r: {hip.min():.1f}..{hip.max():.1f} deg")
