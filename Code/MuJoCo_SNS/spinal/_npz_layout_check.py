"""One-off: verify the spinal_run.npz channel layout (named q vs qfull)
so the curriculum knee-index question is settled against the artifact the
curriculum actually reads."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

z = np.load("spinal_run.npz", allow_pickle=True)
print("keys:", sorted(z.files))
q = z["q"]
print("q shape:", q.shape)
kj = [str(s) for s in z["key_joints"]] if "key_joints" in z.files else None
print("key_joints:", kj)
if "qfull" in z.files:
    print("qfull shape:", z["qfull"].shape)
if kj and q.ndim == 2 and q.shape[1] == len(kj):
    print("q[4] =", kj[4], "| range:", float(q[:, 4].min()), "..",
          float(q[:, 4].max()))
print("t range:", float(z["t"][0]), "..", float(z["t"][-1]))
