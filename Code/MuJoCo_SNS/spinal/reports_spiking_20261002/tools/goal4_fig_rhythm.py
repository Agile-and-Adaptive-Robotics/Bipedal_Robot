"""GOAL 4 figure — constant-DRIVE rhythm: non-spiking vs spiking RG.

Data: goal4_gate2_traces.npz, dumped by tools/goal4_dump_gate2_traces.py
(a verbatim replication of the GATE 2 recipe; valid only if that script's
printed periods match logs/gate2_rhythm2.log: non-spiking 0.399 s,
spiking 0.621 s).

Panels:
  (A) non-spiking: V[RG_E_r] - V[RG_F_r] (mV) — analog alternation,
      period 0.399 +/- 0.164 s.
  (B) spiking mirror: cumulative spike counts of the RAW RG_E_r / RG_F_r
      cells — alternating burst staircases, envelope period 0.621 s.

Output: Figures/30-results/spiking_mirror_rhythm.{pdf,png} + _alt.txt
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

REP = HERE.parent
SPINAL = REP.parent
ROOT = SPINAL.parents[2]
FIGDIR = (ROOT / "Documentation" / "Reports and Papers" / "Dissertation"
          / "Figures" / "30-results")

d = np.load(REP / "goal4_gate2_traces.npz")
t, v_ns, ce, cf = d["t"], d["v_ns"], d["cnt_e"], d["cnt_f"]

apply_style()
fig, ax = plt.subplots(2, 1, figsize=(7.5, 5.2), sharex=True)

ax[0].plot(t, v_ns, color=INDIGO, lw=1.0)
ax[0].set_ylabel("RG ext minus RG flx readout (mV)")
ax[0].set_title("(A) Non-spiking network: analog half-center alternation "
                "(period 0.399 s)")
ax[0].grid(True, lw=0.4, alpha=0.5)

ax[1].plot(t, ce, color=ORANGE, lw=1.2, label="RG ext, right")
ax[1].plot(t, cf, color=INDIGO, lw=1.2, ls="--", label="RG flx, right")
ax[1].set_xlabel("Time (s)")
ax[1].set_ylabel("Cumulative spike count")
ax[1].set_title("(B) Spiking mirror: alternating burst staircases "
                "(envelope period 0.621 s)")
ax[1].legend(loc="upper left", frameon=False)
ax[1].grid(True, lw=0.4, alpha=0.5)

fig.tight_layout()
stem = str(FIGDIR / "spiking_mirror_rhythm")
save_fig(fig, stem)

alt = """Alt text for spiking_mirror_rhythm.pdf

Two panels over twenty seconds of network-only simulation under constant
drive. Panel A shows the difference between the right rhythm-generator
extensor and flexor readout voltages in millivolts for the non-spiking
analog network: an approximately square wave alternating between about
minus 3 and plus 5 millivolts with an average period of 0.399 seconds,
irregular from cycle to cycle. Panel B shows the cumulative spike counts
of the raw right extensor and flexor rhythm-generator cells of the spiking
mirror as two staircase traces: the dashed indigo flexor staircase and the
solid orange extensor staircase rise in alternating bursts, each burst
adding about ten spikes, with the burst-to-burst alternation regular at an
envelope period of 0.621 seconds. The spiking mirror therefore
self-sustains the same alternating rhythm but at a 55.6 percent longer
period than the analog network at this drive.
"""
(FIGDIR / "spiking_mirror_rhythm_alt.txt").write_text(alt, encoding="utf-8")
print("wrote", FIGDIR / "spiking_mirror_rhythm_alt.txt")
