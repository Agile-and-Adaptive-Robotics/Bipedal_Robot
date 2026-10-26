# Comparison plot: intact vs ablated from the two saved npz runs.
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\results"
fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
for arm, col in (("intact", "#0072B2"), ("ablated", "#D55E00")):
    d = np.load(R + f"\\ablation_{arm}_20261026.npz", allow_pickle=True)
    keys = list(d.keys())
    t = d["t"] if "t" in d else d[keys[0]]
    # find the pelvis-height and root-forward channels by name if present
    def pick(*cands):
        for c in cands:
            for k in keys:
                if c in k:
                    return d[k]
        return None
    z = pick("pelz", "root_z", "zpos", "height")
    x = pick("root_x", "xpos", "fwd", "rootx")
    if z is not None:
        axes[0].plot(t, z, col, label=arm, lw=1.4)
    if x is not None:
        axes[1].plot(t, x, col, label=arm, lw=1.4)
axes[0].set_ylabel("pelvis height (m)"); axes[0].legend(); axes[0].grid(alpha=.3)
axes[1].set_ylabel("forward position (m)"); axes[1].set_xlabel("time (s)")
axes[1].grid(alpha=.3)
fig.suptitle("Feedback ablation, frozen tuning: intact (blue) vs contact+Ib ablated (orange)")
fig.tight_layout()
fig.savefig(R + r"\pictures\ablation_pair_20261026.png", dpi=220)
print("wrote pictures\\ablation_pair_20261026.png; npz keys:", keys[:8])
