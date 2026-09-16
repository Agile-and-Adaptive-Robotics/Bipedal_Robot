"""Metrics for the laminated curriculum final run: kine_score, duty,
cadence, joint ranges from ground_curriculum_final.npz (or argv[1])."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import kine_ref
import params as P

npz = sys.argv[1] if len(sys.argv) > 1 else "ground_curriculum_final.npz"
z = np.load(npz, allow_pickle=True)
n_done = int(np.sum(z["t"] > 0)) + 1
t = z["t"][:n_done]
q = z["q"][:n_done]
neuro = z["neuro"][:n_done]
com = z["com"][:n_done]
k = kine_ref.compare(t, q, neuro, P.SCHEDULE["walk"][0])
print(f"== final-run metrics from {npz}")
if k:
    walk_dur = P.SCHEDULE["walk"][1] - P.SCHEDULE["walk"][0]
    print(f"kine_score {k['kine_score']:.3f}  duty {k['duty']:.3f}  "
          f"cadence {k['n_cycles'] / walk_dur:.3f} Hz  "
          f"n_cycles {k['n_cycles']}")
    for j in ("hip", "knee", "ankle"):
        print(f"{j:6s} {k[j + '_min']:7.1f} .. {k[j + '_max']:7.1f} deg  "
              f"(range {k[j + '_max'] - k[j + '_min']:5.1f}, "
              f"rmse {k['rmse_' + j]:.1f})")
    print(f"knee_min {k['knee_min']:.1f}")
jn = [str(a) for a in z["key_joints"]]
for i, name in enumerate(jn):
    print(f"{name:18s} {q[:, i].min():7.1f} .. {q[:, i].max():7.1f} deg")
comz = com[:, 2]
print(f"COM z min {np.nanmin(comz):.3f} (kz), tilt max "
      f"{np.nanmax(q[:, 0]):.1f} deg")
