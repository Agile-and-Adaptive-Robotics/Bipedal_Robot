"""Paper-style figure in the layout of SNS-Toolbox figure_hindlimb:
neural activity (RG half-centers, PF cells, MN output) over joint angles
for one recorded run. Usage: _figure_hindlimb_style.py <run.npz> <tag>
"""
import io
import os
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

HERE = Path(__file__).parent
npz, tag = sys.argv[1], sys.argv[2]
z = np.load(HERE / npz, allow_pickle=True)
t = z["t"]
q = z["q"]            # degrees (KEY_JOINTS order in z["key_joints"])
act = z["act"]        # KEY_ACTS order in z["key_acts"]
neuro = z["neuro"]
acts = [str(a) for a in z["key_acts"]]
joints = [str(a) for a in z["key_joints"]]
neuro_names = [str(a) for a in z["neuro_names"]]

ai = {a: i for i, a in enumerate(acts)}
ji = {a: i for i, a in enumerate(joints)}
ni = {a: i for i, a in enumerate(neuro_names)}
# Window: show only the commanded walk phase.  Air runs discard the first
# second after DRIVE reaches its plateau so the panel emphasizes the settled
# limit cycle without including the intentional ramp-down/holding phase.
kind = "ground" if "ground" in npz else "air"
w0, w1 = (5.0, 15.0) if kind == "ground" else (6.0, 15.0)
m = (t >= w0) & (t <= w1)

fig, ax = plt.subplots(5, 1, figsize=(9.0, 8.2), sharex=True,
                       gridspec_kw=dict(height_ratios=[1, 1, 1, 1, 1]))
fig.subplots_adjust(hspace=0.13, left=0.09, right=0.985, top=0.95,
                    bottom=0.07)

ax[0].plot(t[m], neuro[m, ni["RG_E_r"]], color="#0072b2", lw=1.2,
           label="RG-E (ext)")
ax[0].plot(t[m], neuro[m, ni["RG_F_r"]], color="#d55e00", lw=1.2,
           label="RG-F (flx)")
ax[0].set_ylabel("RG\npotential [mV]", fontsize=8)
ax[1].plot(t[m], neuro[m, ni["PF_E1_r"]], color="#0072b2", lw=0.9,
           label="PF-E1")
ax[1].plot(t[m], neuro[m, ni["PF_E2_r"]], color="#56b4e9", lw=0.9,
           label="PF-E2")
ax[1].plot(t[m], neuro[m, ni["PF_F1_r"]], color="#d55e00", lw=0.9,
           label="PF-F1")
ax[1].plot(t[m], neuro[m, ni["PF_F2_r"]], color="#e69f00", lw=0.9,
           label="PF-F2")
ax[1].set_ylabel("PF\npotential [mV]", fontsize=8)
ax[2].plot(t[m], act[m, ai["vas_lat_r"]], color="#009e73", lw=1.2,
           label="MN knee-ext (vas_lat, activation)")
ax[2].plot(t[m], act[m, ai["semimem_r"]], color="#cc79a7", lw=1.2,
           label="MN knee-flx (semimem, activation)")
ax[2].set_ylabel("knee MN\ndrive [0-1]", fontsize=8)
ax[3].plot(t[m], q[m, ji["knee_angle_r"]], color="black", lw=1.3)
ax[3].set_ylabel("knee angle\n[deg]", fontsize=8)
ax[4].plot(t[m], q[m, ji["hip_flexion_r"]], color="black", lw=1.3)
ax[4].set_ylabel("hip angle\n[deg]", fontsize=8)
ax[4].set_xlabel("time [s]", fontsize=8)

for a in ax:
    a.grid(True, alpha=0.25, lw=0.4)
    a.tick_params(labelsize=7)
    a.set_xlim(w0, w1)
for a in ax[:3]:
    a.legend(fontsize=6.2, ncol=4, loc="upper right", framealpha=0.85)
fig.suptitle(f"gait2392 spinal SNS - right leg ({tag}, RoM limits + "
             f"Renshaw 0.5)", fontsize=9.5)
def save_replace(path, **kwargs):
    """Render beside the target, then replace it after transient preview locks."""
    tmp = path.with_name(path.stem + ".rendering" + path.suffix)
    fig.savefig(tmp, **kwargs)
    for attempt in range(12):
        try:
            os.replace(tmp, path)
            return
        except OSError:
            if attempt == 11:
                raise
            time.sleep(0.25)


save_replace(HERE / "figures" / f"hindlimb_style_{tag}.png", dpi=200)
save_replace(HERE / "figures" / f"hindlimb_style_{tag}.pdf")
print(f"saved figures/hindlimb_style_{tag}.png/.pdf")
