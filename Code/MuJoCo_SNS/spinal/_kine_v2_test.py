"""Direct kine_ref v2 test on the current spinal_run.npz: traceback +
both-leg breakdown of the score."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import kine_ref

z = np.load("spinal_run.npz", allow_pickle=True)
print("npz keys:", sorted(z.files))
t, q, neuro = z["t"], z["q"], z["neuro"]
contact = z["contact"] if "contact" in z.files else None
if contact is not None:
    print(f"contact N: r max {contact[:, 0].max():.1f} mean "
          f"{contact[:, 0].mean():.1f} | l max {contact[:, 1].max():.1f} "
          f"mean {contact[:, 1].mean():.1f}")
    m = t >= 5.0
    print(f"contact>20N frac: r {(contact[m, 0] > 20).mean():.3f} "
          f"l {(contact[m, 1] > 20).mean():.3f}")

try:
    k = kine_ref.compare(t, q, neuro, 5.0,
                         ref=kine_ref.ref_cached(), contact=contact)
    if k is None:
        print("compare returned None (no cycles either leg)")
    else:
        for key in sorted(k):
            v = k[key]
            if isinstance(v, float):
                print(f"  {key:22s} {v:9.3f}")
            else:
                print(f"  {key:22s} {v}")
except Exception:
    import traceback
    traceback.print_exc()
