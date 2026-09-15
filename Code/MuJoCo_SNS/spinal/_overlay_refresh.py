"""Regenerate the OpenSim-overlay gait-cycle figure from a v9 capture
(npz + best json), replacing opensim_overlay_gait_cycles.png.
Usage: python _overlay_refresh.py <capture.npz>
"""
import io
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = Path(__file__).parent
npz = sys.argv[1] if len(sys.argv) > 1 else "v9_best_trial36.npz"
z = np.load(HERE / npz, allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
jj = [str(a) for a in z["key_joints"]]
import kine_ref
mean, n_cyc, duty = kine_ref.sim_cycles(t, q, neuro, 2.0)

ref = kine_ref.ref_cached()
GRID = kine_ref.GRID
fig, ax = plt.subplots(3, 1, figsize=(7.0, 8.6), sharex=True)
fig.subplots_adjust(hspace=0.12, left=0.10, right=0.98, top=0.94,
                    bottom=0.07)
cols = dict(hip="#0072b2", knee="#d55e00", ankle="#009e73")
names = dict(hip="hip [deg]", knee="knee [deg]", ankle="ankle [deg]")
for i, j in enumerate(("hip", "knee", "ankle")):
    for c in (mean[j], ):
        pass
    # individual cycles are not retained in the npz path; plot mean only
    ax[i].plot(GRID, mean[j], color=cols[j], lw=2.2,
               label=f"SNS sim mean ({n_cyc} cyc)")
    ax[i].plot(GRID, ref[j], "--", color="black", lw=1.8,
               label="OpenSim IK")
    ax[i].set_ylabel(names[j], fontsize=9)
    ax[i].grid(True, alpha=0.25, lw=0.4)
    ax[i].tick_params(labelsize=8)
    ax[i].legend(fontsize=7.5, loc="upper right", framealpha=0.85)
ax[-1].set_xlabel("gait cycle [%]", fontsize=9)
fig.suptitle(f"cycle-normalized gait cycles, right leg ({npz}; "
             f"0 = RG-E rise ~ stance start)", fontsize=10)
out = HERE / "figures" / "opensim_overlay_gait_cycles.png"
fig.savefig(out, dpi=170)
print(f"saved {out} (n_cycles {n_cyc}, duty {duty:.2f})")
