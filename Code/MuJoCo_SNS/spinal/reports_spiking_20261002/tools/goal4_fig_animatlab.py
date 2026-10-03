"""GOAL 4 figure — AnimatLab spiking conversion, honest outcome.

Data (existing chart .txt byproducts of the 2026-10-02 campaign, no re-run):
  models/W2L_modern/L Angles.txt           graded baseline (oscillates,
                                           RG period 0.443 s)
  models/W2L_modern_spiking/L Angles.txt   spiking conversion v1 (rhythm
                                           dies; 2 RG spikes in 5 s)
  models/bisect/W2L_neuronsonly/L Angles.txt
                                           thresholds-only ablation (spiking
                                           thresholds, graded synapses kept:
                                           oscillation preserved)

Angles are radians in the chart files; converted to degrees here (same
conversion as tools/analyze.py). Panels: rows = the three variants;
columns = left hip (deg) and left knee (deg).

Output: Figures/93-AppendixB/animatlab_spiking_conversion.{pdf,png} +
_alt.txt (directory created if absent).
"""
import io
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from goal4_fig_style import apply_style, INDIGO, ORANGE, PINK, save_fig  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

REP = HERE.parent            # reports_spiking_20261002
SPINAL = REP.parent          # spinal
ROOT = SPINAL.parents[2]     # repo root (parents: MuJoCo_SNS -> Code -> root)
FIGDIR = (ROOT / "Documentation" / "Reports and Papers"
          / "Dissertation" / "Figures" / "93-AppendixB")
FIGDIR.mkdir(parents=True, exist_ok=True)


def read_angles(path):
    arr = np.loadtxt(path, skiprows=1)
    t = arr[:, 1]
    names = open(path, encoding="utf-8", errors="replace").readline().split()
    cols = {n: i for i, n in enumerate(names)}
    return t, np.degrees(arr[:, cols["hip_L"]]), np.degrees(arr[:, cols["knee_L"]])


rows = [
    ("Graded baseline (oscillates)", "W2L_modern", INDIGO),
    ("Spiking conversion (rhythm dies)", "W2L_modern_spiking", ORANGE),
    # the ablation's charts sit directly in bisect/ (the successful
    # bisect2 re-run; ranges verified against
    # metrics_bisect_W2L_neuronsonly.json: hip -15.1..23.2, knee -0.8..47.9)
    ("Thresholds only, graded synapses (oscillates)", "bisect", PINK),
]

apply_style()
fig, ax = plt.subplots(3, 2, figsize=(7.5, 7.0), sharex=True)
for r, (title, sub, color) in enumerate(rows):
    t, hip, knee = read_angles(REP / "models" / sub / "L Angles.txt")
    for c, (sig, ylab) in enumerate([(hip, "Left hip angle (deg)"),
                                     (knee, "Left knee angle (deg)")]):
        a = ax[r, c]
        a.plot(t, sig, color=color, lw=1.1)
        a.set_title(f"({'ABC'[r]}{1 + c}) {title}" if c == 0 else title)
        a.set_ylabel(ylab, fontsize=10)
        a.grid(True, lw=0.4, alpha=0.5)
        if r == 2:
            a.set_xlabel("Time (s)")
fig.tight_layout()

stem = str(FIGDIR / "animatlab_spiking_conversion")
save_fig(fig, stem)

alt = """Alt text for animatlab_spiking_conversion.pdf

Six panels in three rows and two columns showing left-leg joint angles in
degrees over five seconds of the AnimatLab walker family. The left column
is the left hip angle, and the right column is the left knee angle. Row A,
in indigo, is the graded non-spiking baseline: the hip oscillates smoothly
between about minus 15 and plus 23 degrees and the knee between about minus
1 and plus 60 degrees at a rhythm-generator period of 0.443 seconds. Row B,
in orange, is the direct spiking conversion with thresholds set to the
spiking regime and all chemical synapses converted to per-spike synapses:
both traces are nearly flat, with the hip pinned near 6 degrees and the
knee near 11 to 12 degrees; the rhythm dies. Row C, in pink, is the
thresholds-only ablation, in which the neurons carry spiking-regime
thresholds but the graded synapses are kept: the oscillation returns with
essentially baseline hip and knee ranges. The figure demonstrates that the
synapse conversion, not the spiking thresholds, is what breaks the rhythm.
"""
(FIGDIR / "animatlab_spiking_conversion_alt.txt").write_text(alt,
                                                            encoding="utf-8")
print("wrote", FIGDIR / "animatlab_spiking_conversion_alt.txt")
