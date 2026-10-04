"""GOAL 4 figure — ground-gait comparison, non-spiking s3k winner vs the
UNTUNED spiking mirror (both at the s3k winner parameters, 16 s ground runs).

Data (existing artifacts, no re-run):
  spinal/spinal_run_spkbase_s3k.npz       (kine_score -161.56754173676563)
  spinal/spinal_run_spkbase_spiking.npz   (kine_score -236.3496680668847)

Panels:
  (A) knee angle, right (dashed orange) + left (solid indigo), non-spiking
      s3k: left leg frozen near extension.
  (B) knee angle, spiking mirror: both legs cycle (8 cycles right,
      9 cycles left; cycle knee minimum -79.8 deg).
  (C) hip angle, non-spiking.  (D) hip angle, spiking mirror.

Output: Figures/30-results/spiking_mirror_ground_gait.{pdf,png} + _alt.txt
"""
import io
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from goal4_fig_style import apply_style, INDIGO, ORANGE, save_fig  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

SPINAL = HERE.parents[1]
FIGDIR = (SPINAL.parents[2] / "Documentation" / "Reports and Papers"
          / "Dissertation" / "Figures" / "30-results")

apply_style()

JN = {j: i for i, j in enumerate(
    np.load(SPINAL / "spinal_run_spkbase_s3k.npz", allow_pickle=True)["key_joints"])}


def knee_hip(npz):
    d = np.load(npz, allow_pickle=True)
    q, t = d["q"], d["t"]
    return (t, q[:, JN["knee_angle_r"]], q[:, JN["knee_angle_l"]],
            q[:, JN["hip_flexion_r"]], q[:, JN["hip_flexion_l"]])


(t, kr_ns, kl_ns, hr_ns, hl_ns) = knee_hip(SPINAL / "spinal_run_spkbase_s3k.npz")
(t2, kr_sp, kl_sp, hr_sp, hl_sp) = knee_hip(
    SPINAL / "spinal_run_spkbase_spiking.npz")
assert np.allclose(t, t2)

fig, ax = plt.subplots(2, 2, figsize=(7.5, 6.0), sharex=True)

pairs = [
    (ax[0, 0], kr_ns, kl_ns, "Knee angle (deg)", "(A) Non-spiking s3k winner"),
    (ax[0, 1], kr_sp, kl_sp, "Knee angle (deg)", "(B) Spiking mirror, untuned"),
    (ax[1, 0], hr_ns, hl_ns, "Hip flexion angle (deg)", "(C) Non-spiking s3k winner"),
    (ax[1, 1], hr_sp, hl_sp, "Hip flexion angle (deg)", "(D) Spiking mirror, untuned"),
]
for a, r, l, ylab, title in pairs:
    a.plot(t, l, color=INDIGO, lw=1.1, label="left")
    a.plot(t, r, color=ORANGE, lw=1.1, ls="--", label="right")
    a.set_ylabel(ylab)
    a.set_title(title)
    a.grid(True, lw=0.4, alpha=0.5)
for a in ax[1, :]:
    a.set_xlabel("Time (s)")
ax[0, 0].legend(loc="upper right", frameon=False)
fig.tight_layout()

stem = str(FIGDIR / "spiking_mirror_ground_gait")
save_fig(fig, stem)

alt = """Alt text for spiking_mirror_ground_gait.pdf

Four panels of joint-angle traces over a 16-second ground walk of the
converted MuJoCo robot, comparing the non-spiking spinal-network winner
(left column) with the untuned spiking-neuron mirror run at the same winner
parameters (right column). Top row: knee angle in degrees. In panel A
(non-spiking), the dashed orange right-leg trace sweeps repeatedly between
roughly minus 100 and plus 8 degrees while the solid indigo left-leg trace
stays nearly flat between minus 42 and plus 15 degrees: the left leg is
frozen and only the right leg steps. In panel B (spiking mirror), both the
indigo left trace and the dashed orange right trace sweep with deep
flexion, the left leg completing nine cycles and the right leg eight.
Bottom row: hip flexion angle in degrees. In panel C (non-spiking), the
right hip oscillates between about minus 18 and plus 44 degrees while the
left hip moves only a few degrees. In panel D (spiking mirror), both hips
oscillate between roughly minus 35 and plus 43 degrees. The figure shows
that the untuned spiking mirror, although it scores worse overall, produces
a bilateral gait that the non-spiking winner cannot.
"""
(FIGDIR / "spiking_mirror_ground_gait_alt.txt").write_text(alt,
                                                          encoding="utf-8")
print("wrote", FIGDIR / "spiking_mirror_ground_gait_alt.txt")
