"""Reporter verification probes (2026-10-01), laptop myoconv env.

Probe A: build_w2l_split_net.build(comm=1.0) census — expect 87 neurons / 164 synapses
         (the census fix both the MuJoCo-port and BPA tracks rely on).
Probe B: full W2L net under tonic+kickoff ONLY (no heel contact trains) — expect a
         latch (zero RG-E bursts), the finding that forced the contact-train protocol
         in the SNS_W2L_CPG validation. Mirrors campaigns/20260930/w2l_numpy_ref_20260930.py
         with the contact-train lines removed.
"""
from __future__ import annotations

import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

SPINAL = Path(__file__).resolve().parents[1] / "Code" / "MuJoCo_SNS" / "spinal"
sys.path.insert(0, str(SPINAL / "w2l_mujoco"))
sys.path.insert(0, str(SPINAL / "w2l_cpg"))

# --- Probe A: split-net census -------------------------------------------------
from build_w2l_split_net import build as build_split  # noqa: E402

net = build_split(comm=1.0)
print(f"PROBE_A SPLIT_CENSUS neurons={len(net.idx)} synapses={net.n_synapses}")

# --- Probe B: tonic + kickoff only, no contact trains ---------------------------
from build_w2l_net import DT, build  # noqa: E402

DURATION, SKIP = 12.0, 2.0
full = build()
nsteps = int(round(DURATION / DT))
u = full.make_inputs()
i_s1 = full.input_index("Stimulus_1")
i_s2 = full.input_index("Stimulus_2 (10nA)")
i_te = [full.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
i_tf = [full.input_index("TONIC " + w) for w in ("L RG flx", "R RG flx")]
idx = full.idx
lE = np.zeros(nsteps)
rE = np.zeros(nsteps)
for k in range(nsteps):
    t = k * DT
    u[:] = 0.0
    for i in i_te:
        u[i] = 2.0
    for i in i_tf:
        u[i] = 3.0
    if t < 0.01:
        u[i_s1] = 10.0
        u[i_s2] = 10.0
    V = full.step(u)
    lE[k] = V[idx["L RG ext"]]
    rE[k] = V[idx["R RG ext"]]

win = np.arange(nsteps) * DT >= SKIP


def bursts(sig):
    on = sig > 0.5 * sig.max()
    return np.flatnonzero(on[1:] & ~on[:-1]) + 1


sL, sR = bursts(lE[win]), bursts(rE[win])
print(f"PROBE_B TONIC_ONLY window {SKIP:.0f}-{DURATION:.0f}s "
      f"L RG-E max {lE[win].max():.3f} mV bursts {len(sL)} | "
      f"R RG-E max {rE[win].max():.3f} mV bursts {len(sR)}")
print("PROBE_B VERDICT:", "LATCH (no rhythm)" if (len(sL) < 4 and len(sR) < 4) else "OSCILLATES")
