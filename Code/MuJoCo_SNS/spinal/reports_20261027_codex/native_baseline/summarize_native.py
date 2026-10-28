"""Summarize the native baseline trajectory (eval_spinal_native.sto.sto)."""
import json

import numpy as np

PATH = "eval_spinal_native.sto.sto"
with open(PATH) as fh:
    lines = [fh.readline() for _ in range(7)]
names = lines[6].strip().split("\t")
data = np.loadtxt(PATH, skiprows=7)
col = {n: i for i, n in enumerate(names)}

t = data[:, col["time"]]
tx = data[:, col["pelvis_tx"]]
ty = data[:, col["pelvis_ty"]]
tilt = np.degrees(data[:, col["pelvis_tilt"]])

summary = {
    "rows": int(len(t)),
    "dt_storage": float(t[1] - t[0]),
    "duration_s": float(t[-1]),
    "pelvis_tx_final_m": float(tx[-1]),
    "pelvis_tx_range_m": [float(tx.min()), float(tx.max())],
    "mean_speed_mps_free_1s": float(np.polyfit(t[t >= 1.0], tx[t >= 1.0], 1)[0]),
    "pelvis_ty_range_m": [float(ty.min()), float(ty.max())],
    "pelvis_tilt_deg_range": [float(tilt.min()), float(tilt.max())],
    "phase_left_range": [float(data[:, col["phase_left"]].min()),
                         float(data[:, col["phase_left"]].max())],
    "phase_right_range": [float(data[:, col["phase_right"]].min()),
                          float(data[:, col["phase_right"]].max())],
}

for leg, pref in (("L", "leg0_l"), ("R", "leg1_r")):
    gy = data[:, col[pref + ".grf_norm_y"]]
    gz = data.get if False else None
    summary[f"grf_norm_y_{leg}"] = {
        "max": float(gy.max()),
        "stance_fraction_gt_0.1": float(np.mean(gy > 0.1)),
    }

# knee/ankle/hip ranges (rad) for pose record
for jn in ("hip_flexion_r", "knee_angle_r", "ankle_angle_r",
           "hip_flexion_l", "knee_angle_l", "ankle_angle_l"):
    v = data[:, col[jn]]
    summary[f"{jn}_deg_range"] = [float(np.degrees(v.min())),
                                  float(np.degrees(v.max()))]

# a couple of muscle excitation stats (regression-target availability)
for m in ("soleus_r", "tib_ant_r", "vasti_l"):
    ex = data[:, col[m + ".excitation"]]
    summary[f"{m}.excitation_range"] = [float(ex.min()), float(ex.max())]

with open("native_baseline_summary.json", "w") as fh:
    json.dump(summary, fh, indent=1)
print(json.dumps(summary, indent=1))
