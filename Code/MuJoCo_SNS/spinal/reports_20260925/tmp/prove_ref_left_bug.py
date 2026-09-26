"""Show the left-reference bug: current (edge-filled) vs corrected (IK-covered cycle).

Read-only vs kine_ref; writes one PNG into reports_20260925/figs/.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import kine_ref  # noqa: E402

t_ik, names, vals = kine_ref._read_mot(kine_ref.IK_MOT)
t_g, gn, gv = kine_ref._read_mot(kine_ref.GRF_MOT)
vy_l = gv[:, gn.index("1_ground_force_vy") - 1]
on_l = kine_ref._loading_onsets(t_g, vy_l)
col = {n: i for i, n in enumerate(names[1:])}

ref = kine_ref.ref_cached()

# Corrected left cycle: first onset pair FULLY inside IK support
ok = [(a, b) for a, b in zip(on_l[:-1], on_l[1:]) if a >= t_ik[0] and b <= t_ik[-1]]
t0, t1 = ok[0]
m = (t_ik >= t0) & (t_ik < t1)
fix = {j.split("_")[0]: np.interp(
    kine_ref.GRID, (t_ik[m] - t0) / (t1 - t0) * 100.0, vals[m, col[j]])
    for j in ("hip_flexion_l", "knee_angle_l", "ankle_angle_l")}

fig, axes = plt.subplots(3, 1, figsize=(7.5, 8.5), sharex=True)
for ax, j in zip(axes, ("hip", "knee", "ankle")):
    ax.plot(kine_ref.GRID, ref["l"][j], "k--", lw=1.6,
            label="as-built ref l (frozen 0-40%: IK starts 0.5 s, cycle starts 0.005 s)")
    ax.plot(kine_ref.GRID, fix[j], color="tab:red", lw=1.8,
            label=f"corrected: left cycle {t0:.3f}-{t1:.3f} s (fully IK-covered)")
    ax.set_ylabel(f"{j} (deg, +flex/DF)")
    ax.grid(alpha=0.3)
axes[0].legend(fontsize=8, loc="lower right")
axes[2].set_xlabel("left gait cycle (%)")
axes[0].set_title("subject01_walk1 LEFT reference cycle - the np.interp edge-fill bug Ben spotted\n"
                  "IK file covers 0.50-2.50 s; the as-built left cycle starts at 0.005 s "
                  "=> first 39.6% is a held constant")
fig.tight_layout()
out = Path(__file__).resolve().parents[1] / "figs" / "ref_left_bug_proof.png"
fig.savefig(out, dpi=130)
print("wrote", out)
print(f"corrected left cycle spans {t0:.4f}-{t1:.4f} s; "
      f"hip range {np.ptp(fix['hip']):.1f} deg (as-built {np.ptp(ref['l']['hip']):.1f})")
