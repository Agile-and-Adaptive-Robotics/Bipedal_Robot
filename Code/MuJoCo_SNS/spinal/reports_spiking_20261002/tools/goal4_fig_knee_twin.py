"""GOAL 4 figure — Simulink knee-reflex twin: committed non-spiking
KneeReflexDemo vs KneeReflexDemo_Spiking at the matched operating point.

Data: goal4_knee_traces.mat (exported this session by
goal4_export_knee_traces.m; printed metrics matched goal2_knee_spiking.mat:
baseline rise 15.79 deg, settle 43.5 (41.6-44.9); twin rise 15.79 deg,
settle 44.4 (41.9-46.7)).

Panels:
  (A) knee angle, both models overlaid over the full 5 s.
  (B) last 2 s window: the spiking twin's settle band is slightly hotter
      (41.9-46.7 vs 41.6-44.9 deg).

Output: Figures/30-results/simulink_knee_spiking_twin.{pdf,png} + _alt.txt
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
import scipy.io as sio  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

REP = HERE.parent
SPINAL = REP.parent
ROOT = SPINAL.parents[2]
FIGDIR = (ROOT / "Documentation" / "Reports and Papers" / "Dissertation"
          / "Figures" / "30-results")

d = sio.loadmat(REP / "goal4_knee_traces.mat", squeeze_me=True,
                struct_as_record=False)
base, spk = d["base"], d["spk"]
tb = np.asarray(base.t); thb = np.degrees(np.asarray(base.th))
ts = np.asarray(spk.t); ths = np.degrees(np.asarray(spk.th))

apply_style()
fig, ax = plt.subplots(2, 1, figsize=(7.5, 5.2))

ax[0].plot(tb, thb, color=INDIGO, lw=1.2, label="Non-spiking demo")
ax[0].plot(ts, ths, color=ORANGE, lw=1.2, ls="--", label="Spiking twin")
ax[0].set_ylabel("Knee angle (deg)")
ax[0].set_title("(A) Full 5 s run at the matched operating point")
ax[0].legend(loc="lower right", frameon=False)
ax[0].grid(True, lw=0.4, alpha=0.5)

w = ts > ts[-1] - 2.0
ax[1].plot(ts[w], ths[w], color=ORANGE, lw=1.2, ls="--", label="Spiking twin")
wb = tb > tb[-1] - 2.0
ax[1].plot(tb[wb], thb[wb], color=INDIGO, lw=1.2, label="Non-spiking demo")
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Knee angle (deg)")
ax[1].set_title("(B) Settle window (last 2 s)")
ax[1].legend(loc="lower right", frameon=False)
ax[1].grid(True, lw=0.4, alpha=0.5)

fig.tight_layout()
stem = str(FIGDIR / "simulink_knee_spiking_twin")
save_fig(fig, stem)

alt = """Alt text for simulink_knee_spiking_twin.pdf

Two panels of knee-angle traces in degrees from the Simulink knee-reflex
demonstration. Panel A shows the full five-second run: the solid indigo
trace is the committed non-spiking demonstration and the dashed orange
trace is its spiking twin at the matched operating point; the two traces
rise together from zero to about 16 degrees at 0.05 seconds, overshoot to
about 47 degrees, and settle near 44 degrees, lying on top of one another
for the whole run. Panel B zooms into the final two seconds: both traces
alternated around their settle level, with the non-spiking demo oscillating
between about 41.6 and 44.9 degrees and the spiking twin between about 41.9
and 46.7 degrees, a slightly wider and faster ripple. The figure shows that
replacing the sensory interneuron layer with spiking cells leaves the
plant-side knee response essentially unchanged at the matched operating
point.
"""
(FIGDIR / "simulink_knee_spiking_twin_alt.txt").write_text(alt,
                                                          encoding="utf-8")
print("wrote", FIGDIR / "simulink_knee_spiking_twin_alt.txt")
