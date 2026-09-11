"""Phase-aligned activation averages from the last spinal_run.npz: mean
activation of key muscles during RG-E (stance) vs RG-F (swing) windows of
the right leg. Shows exactly who co-contracts when the knee should swing."""
import numpy as np

d = np.load("spinal_run.npz", allow_pickle=True)
t, act, neuro = d["t"], d["act"], d["neuro"]
key_acts = list(d["key_acts"])
names = ("DRIVE", "POSTURE", "RG_E_r", "RG_F_r", "PF_E2_r", "PF_F1_r",
         "RG_E_l", "RG_F_l")
walk = (t > 6.0) & (t < 14.5)
rgE, rgF = neuro[:, 2], neuro[:, 3]
# adaptive thresholds: RG amplitude scales with DRIVE
e_hi = 0.6 * np.max(rgE[walk])
f_hi = 0.6 * np.max(rgF[walk])
stance_m = walk & (rgE > e_hi) & (rgF < 0.5 * f_hi)
swing_m = walk & (rgF > f_hi) & (rgE < 0.5 * e_hi)
print(f"stance samples {stance_m.sum()}, swing samples {swing_m.sum()}")
print(f"{'muscle':14s} {'stance':>7s} {'swing':>7s}")
for a in ("vas_lat_r", "rect_fem_r", "semimem_r", "med_gas_r", "soleus_r",
          "tib_ant_r", "psoas_r", "glut_max2_r"):
    i = key_acts.index(a)
    print(f"{a:14s} {act[stance_m, i].mean():7.2f} {act[swing_m, i].mean():7.2f}")
knee = d["q"][:, 1]
print(f"\nknee_r mean stance {np.degrees(np.nanmean(np.where(stance_m, np.radians(knee), np.nan))):+.2f} deg"
      if False else
      f"knee_r mean stance {knee[stance_m].mean():+.2f} deg, swing {knee[swing_m].mean():+.2f} deg")
