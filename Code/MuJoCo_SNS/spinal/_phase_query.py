import numpy as np

d = np.load("spinal_run.npz", allow_pickle=True)
t, act, neuro = d["t"], d["act"], d["neuro"]
names = list(d["key_acts"])
walk = (t > 6.0) & (t < 14.5)
rgE, rgF = neuro[:, 2], neuro[:, 3]
e_hi, f_hi = 0.6 * np.max(rgE[walk]), 0.6 * np.max(rgF[walk])
st = walk & (rgE > e_hi) & (rgF < 0.5 * f_hi)
sw = walk & (rgF > f_hi) & (rgE < 0.5 * e_hi)
for a in ("rect_fem_r", "sar_r", "grac_r", "psoas_r", "semimem_r",
          "glut_med1_r", "tib_ant_r"):
    i = names.index(a)
    print(f"{a:12s} stance {act[st, i].mean():5.2f}  swing "
          f"{act[sw, i].mean():5.2f}  peak {act[walk, i].max():5.2f}")
