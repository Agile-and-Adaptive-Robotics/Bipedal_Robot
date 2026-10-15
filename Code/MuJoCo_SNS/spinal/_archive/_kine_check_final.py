"""Final consistency check: compare() on the on-disk winner npz."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import kine_ref as KR

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro, contact = z["t"], z["q"], z["neuro"], z["contact"]
for ws in (5.0, 2.0):
    k = KR.compare(t, q, neuro, ws, ref=KR.ref_cached(), contact=contact)
    if k is None:
        print(f"walk_start {ws}: None")
        continue
    print(f"walk_start {ws}: kine_score {k['kine_score']:.4f} "
          f"duty_r {k['duty_r']:.4f} duty_l {k.get('duty_l', float('nan')):.4f} "
          f"n_cyc_r {k.get('n_cycles_r')} ds {k.get('ds', float('nan')):.3f}")
