"""Where in the run did stepping actually happen? 4-second sliding windows
over spinal_run.npz: knee range, RG_E/RG_F burst counts."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(s) for s in z["neuro_names"]]
kj = [str(s) for s in z["key_joints"]]
ik, ih = kj.index("knee_angle_r"), kj.index("hip_flexion_r")
ie, ifl = names.index("RG_E_r"), names.index("RG_F_r")
if "cfg" in z.files:
    print("cfg:", str(z["cfg"]))
w = 4.0
print(f"{'t0':>5} {'knee_r range':>14} {'hip_r range':>13} "
      f"{'E bursts':>9} {'F bursts':>9}")
for t0 in np.arange(0.0, t[-1] - w, 2.0):
    m = (t >= t0) & (t < t0 + w)
    if not m.any():
        continue
    kr = q[m, ik]
    hr = q[m, ih]
    e = neuro[m, ie]
    f = neuro[m, ifl]

    def bursts(sig):
        thr = 0.5 * max(sig.max(), 1e-9)
        if sig.max() < 1.0:
            return 0
        on = (sig > thr).astype(int)
        return int(np.sum(np.diff(on) == 1))
    print(f"{t0:5.0f} {kr.min():6.1f}..{kr.max():6.1f} "
          f"{hr.min():5.1f}..{hr.max():5.1f} {bursts(e):9d} {bursts(f):9d}")
