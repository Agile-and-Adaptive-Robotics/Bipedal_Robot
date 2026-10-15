"""Check the RG collapse Ben circled in hindlimb_style_nap_air.png:
when do RG_E/RG_F stop alternating, relative to the DRIVE schedule
(walk ends 15 s, ramp_down 15-17 s, stand2 17-22 s)?"""
import io
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

z = np.load("nap_air_walk.npz", allow_pickle=True)
t, neuro = z["t"], z["neuro"]
names = [str(n) for n in z["neuro_names"]]
i_drive = names.index("DRIVE")
i_rge, i_rgf = names.index("RG_E_r"), names.index("RG_F_r")

print(" t[s]  DRIVE   RG_E_r  RG_F_r")
for tt in np.arange(10.0, 22.01, 0.5):
    k = int(np.argmin(np.abs(t - tt)))
    print(f"{tt:5.1f}  {neuro[k, i_drive]:6.2f}  {neuro[k, i_rge]:6.2f}  "
          f"{neuro[k, i_rgf]:6.2f}")

# last time the differential swings > 2 mV (alternation alive)
d = neuro[:, i_rge] - neuro[:, i_rgf]
alive = t[np.abs(np.diff(np.sign(np.round(d, 1)))) > 0]
swings = t[np.abs(d) > 2.0]
print(f"\nlast |RG_E-RG_F| > 2 mV at t = {swings.max():.2f} s "
      f"(walk window ends 15.0, ramp_down ends 17.0)")
print(f"final state: RG_E {neuro[-1, i_rge]:.2f} mV, "
      f"RG_F {neuro[-1, i_rgf]:.2f} mV (rest = 0)")
