"""Why can't the left leg lift? (t54 winner npz, walk window t>=2)
Left vs right: RG, PF cells, MN activations of the swing-critical
muscles (knee flexors semimem/bifemsh, hip flexors psoas/sar/rect_fem,
ankle DF tib_ant) vs the stance extensors (vas_lat, soleus, glut_max).
If left flexor MN drive IS present and strong, the latch is mechanical/
extensor co-contraction; if absent, it's wiring/gain."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro, act, contact = z["t"], z["q"], z["neuro"], z["act"], \
    z["contact"]
names = [str(s) for s in z["neuro_names"]]
acts = [str(s) for s in z["key_acts"]]
kj = [str(s) for s in z["key_joints"]]
m = t >= 2.0

print("=== neural (mV): right vs left ===")
for ch in ("RG_E_r", "RG_F_r", "PF_E1_r", "PF_F1_r", "PF_F2_r",
           "RG_E_l", "RG_F_l", "PF_E1_l", "PF_F1_l", "PF_F2_l"):
    if ch in names:
        v = neuro[m, names.index(ch)]
        print(f"{ch:9s} {v.min():6.2f}..{v.max():6.2f}  swing "
              f"{v.max() - v.min():5.2f}")

print("\n=== MN activations: swing-critical vs extensors ===")
groups = [("knee_flex", ("semimem_r", "semem?" )), ]
WATCH = ("semimem", "bifemsh", "psoas", "sar", "rect_fem", "tib_ant",
         "vas_lat", "soleus", "med_gas", "glut_max2")
for base in WATCH:
    for s in ("r", "l"):
        a = f"{base}_{s}"
        if a in acts:
            v = act[m, acts.index(a)]
            print(f"{a:12s} mean {v.mean():.3f}  max {v.max():.3f}")

print("\n=== joints + contact (deg / N) ===")
for j in ("knee_angle_l", "hip_flexion_l", "ankle_angle_l",
          "knee_angle_r"):
    i = kj.index(j)
    print(f"{j:14s} {q[m, i].min():7.1f}..{q[m, i].max():6.1f}  "
          f"mean {q[m, i].mean():6.1f}")
print(f"contact r: mean {contact[m, 0].mean():6.1f} N  min "
      f"{contact[m, 0].min():6.1f}  max {contact[m, 0].max():6.1f}")
print(f"contact l: mean {contact[m, 1].mean():6.1f} N  min "
      f"{contact[m, 1].min():6.1f}  max {contact[m, 1].max():6.1f}")
