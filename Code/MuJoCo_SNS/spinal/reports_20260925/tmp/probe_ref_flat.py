"""Probe: is the LEFT OpenSim reference cycle flat because IK coverage starts late?

Checks (read-only):
  1. IK file time support vs GRF loading onsets (right/left).
  2. ref['l'] / ref['r'] joint traces: plateau detection (repeated first value).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))  # spinal/
import kine_ref  # noqa: E402

t_ik, names, vals = kine_ref._read_mot(kine_ref.IK_MOT)
t_g, gn, gv = kine_ref._read_mot(kine_ref.GRF_MOT)
vy_r = gv[:, gn.index("ground_force_vy") - 1]
vy_l = gv[:, gn.index("1_ground_force_vy") - 1]
on_r = kine_ref._loading_onsets(t_g, vy_r)
on_l = kine_ref._loading_onsets(t_g, vy_l)

print(f"IK   time support : {t_ik[0]:.4f} .. {t_ik[-1]:.4f} s  ({len(t_ik)} rows)")
print(f"GRF  time support : {t_g[0]:.4f} .. {t_g[-1]:.4f} s  ({len(t_g)} rows)")
print(f"right loading onsets ({len(on_r)}): {np.round(on_r[:5], 4)} ...")
print(f"left  loading onsets ({len(on_l)}): {np.round(on_l[:5], 4)} ...")
print()
for side, onsets in (("r", on_r), ("l", on_l)):
    t0, t1 = onsets[0], onsets[1]
    cover0 = 100.0 * (max(t_ik[0], t0) - t0) / (t1 - t0)
    cover1 = 100.0 * (min(t_ik[-1], t1 - 1e-9) - t0) / (t1 - t0)
    print(f"[{side}] cycle onsets {t0:.4f} -> {t1:.4f}  (T={t1 - t0:.3f} s)")
    print(f"[{side}] IK-covered phase span: {cover0:.1f}% .. {cover1:.1f}%"
          f"   <-- outside this, np.interp HOLDS the edge value")
    if side == "l":
        print(f"[l] IK starts {'BEFORE' if t_ik[0] <= t0 else 'AFTER'} left cycle start"
              f" by {abs(t_ik[0] - t0):.4f} s"
              f" = {100 * abs(t_ik[0] - t0) / (t1 - t0):.1f}% of the cycle")
print()
ref = kine_ref.ref_cached()
for side in ("r", "l"):
    for j in ("hip", "knee", "ankle"):
        y = ref[side][j]
        n_head = len(y) - np.flatnonzero(np.diff(y) != 0).min() if np.any(np.diff(y) != 0) else len(y)
        head_flat = np.all(y[:50] == y[0])
        print(f"ref[{side}][{j:5s}]  y[0]={y[0]:8.3f}  flat through first 50%? {head_flat}"
              f"   first {n_head} samples identical"
              f"   range={np.ptp(y):6.2f} deg")
