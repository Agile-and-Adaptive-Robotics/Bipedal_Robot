"""Peak walking joint-torque DEMAND from the bsolve ID residuals
(tau_fit: per-frame fit rows = net muscle moment required at each joint,
subject01_walk1 + measured GRF, computed the right way in bsolve_ik).
"""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

z = np.load("bsolve_out.npz", allow_pickle=True)
names = [n for n in z["fit_joints"]]
tau = z["tau_fit"]          # [T, nrows]
t = z["t"]
print(f"frames {len(t)}, fit rows: {names}")
mass = 78.5
for jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
           "lumbar_extension", "hip_flexion_l", "knee_angle_l",
           "ankle_angle_l"):
    if jn in names:
        i = names.index(jn)
        col = tau[:, i]
        pk = int(np.argmax(np.abs(col)))
        print(f"  {jn:20s} peak |tau| {abs(col[pk]):6.1f} N*m "
              f"({abs(col[pk]) / mass:.2f} N*m/kg) at t={t[pk]:.2f} s "
              f"(sign {col[pk]:+.1f})")
