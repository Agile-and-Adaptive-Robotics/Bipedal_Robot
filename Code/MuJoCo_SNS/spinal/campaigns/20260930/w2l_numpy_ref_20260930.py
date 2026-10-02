"""Numpy reference for the SNS_W2L_CPG Simulink core validation (2026-10-01).

Runs the FULL transcribed W2L bilateral CPG (spinal/w2l_cpg/build_w2l_net.py,
95 neurons / 208 synapses / 10 ports) under the SAME input protocol the
Simulink representative core receives:
  - constant tonic drive on all four RG half-centers (E 2.0 nA, F 3.0 nA -
    smoke_w2l.py defaults, stack DRIVE convention),
  - the verbatim .aproj antiphase kickoff pair: Stimulus_1 (-> L RG ext) and
    Stimulus_2 (-> R RG flx), 10 nA for t in [0, 0.01 s],
  - alternating heel-contact trains (smoke_w2l.py defaults: 1.0 Hz, 0.30
    active fraction, 20 nA adapter current): L heel active phase [0, 0.3),
    R heel [0.5, 0.8).
    PROVEN NECESSARY 2026-10-01: under tonic+kickoff only the full net LATCHES
    (R RG-E max -1.154 mV over 2-12 s, zero bursts) - the W2L architecture is
    contact-driven (its AnimatLab template has no tonic-drive wiring at all;
    the TONIC ports are this package's runner convention). The Simulink core
    receives the SAME heel trains through the identical graded contact-neuron
    construction (tau 0.04, output synapses g=1.0 onto the ipsilateral
    extensor layers).
Duration 12 s at dt = 2 ms (build_w2l_net.DT), fixed-tau_h numpy backend
(SNS_NumpyFixedTau, exactly the compiled class the smoke uses).

Analysis = smoke_w2l.py / check_rhythm.py conventions: drop the first 2 s,
burst starts = rising edges of (V > 0.5 * window max) on each RG-E trace,
period = mean diff * dt, antiphase = Pearson r of the two RG-E traces.

Writes w2l_numpy_ref.json + w2l_numpy_ref.npz next to this file.
Exit 0 iff >= 4 RG-E bursts on BOTH sides and r < -0.5.
"""
from __future__ import annotations

import io
import json
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = Path(__file__).parent
SPINAL = HERE.parent.parent
sys.path.insert(0, str(SPINAL / "w2l_cpg"))

from build_w2l_net import DT, build  # noqa: E402

DURATION = 12.0
SKIP = 2.0
TONIC_E, TONIC_F = 2.0, 3.0
CONTACT_CURRENT = 20.0     # nA (AnimatLab adapter Gain C=20)
CONTACT_FREQ = 1.0         # Hz
CONTACT_WIDTH = 0.30       # active fraction of each heel pulse


def main() -> int:
    net = build()
    print(f"full W2L net: {len(net.idx)} neurons, {net.n_synapses} synapses, "
          f"{len(net.inputs)} input ports")
    nsteps = int(round(DURATION / DT))
    u = net.make_inputs()
    i_s1 = net.input_index("Stimulus_1")
    i_s2 = net.input_index("Stimulus_2 (10nA)")
    i_te = [net.input_index("TONIC " + w) for w in ("L RG ext", "R RG ext")]
    i_tf = [net.input_index("TONIC " + w) for w in ("L RG flx", "R RG flx")]
    i_hl = net.input_index("L heel contact")
    i_hr = net.input_index("R heel contact")

    watch = ["L RG ext", "L RG flx", "R RG ext", "R RG flx"]
    V_all = np.zeros((nsteps, len(net.idx)))
    period_s = 1.0 / CONTACT_FREQ
    for k in range(nsteps):
        t = k * DT
        u[:] = 0.0
        for i in i_te:
            u[i] = TONIC_E
        for i in i_tf:
            u[i] = TONIC_F
        if t < 0.01:
            u[i_s1] = 10.0
            u[i_s2] = 10.0
        phase = (t % period_s) / period_s
        if phase < CONTACT_WIDTH:
            u[i_hl] = CONTACT_CURRENT
        if 0.5 <= phase < 0.5 + CONTACT_WIDTH:
            u[i_hr] = CONTACT_CURRENT
        V = net.step(u)
        if not np.isfinite(V).all():
            print(f"W2L_NUMPY_REF FAIL non-finite V at t={t:.3f} s")
            return 1
        V_all[k] = V

    win = np.arange(nsteps) * DT >= SKIP
    idx = net.idx
    lE, rE = V_all[win, idx["L RG ext"]], V_all[win, idx["R RG ext"]]

    def bursts(sig):
        on = sig > 0.5 * sig.max()
        return np.flatnonzero(on[1:] & ~on[:-1]) + 1

    sL, sR = bursts(lE), bursts(rE)
    nL, nR = len(sL), len(sR)
    r = float(np.corrcoef(lE, rE)[0, 1])
    perL = float(np.diff(sL).mean() * DT) if nL >= 2 else float("nan")
    perR = float(np.diff(sR).mean() * DT) if nR >= 2 else float("nan")

    print(f"protocol: tonic E {TONIC_E} / F {TONIC_F} nA + kickoff pair + "
          f"alternating {CONTACT_FREQ:g} Hz heel trains "
          f"({CONTACT_CURRENT:g} nA, width {CONTACT_WIDTH}), "
          f"{DURATION:.0f} s @ dt={DT}")
    print(f"window {SKIP:.0f}-{DURATION:.0f} s | L RG E max {lE.max():.3f} mV, "
          f"R RG E max {rE.max():.3f} mV")
    print(f"bursts: L={nL} (period {perL:.3f} s)  R={nR} (period {perR:.3f} s)")
    print(f"left-right RG-E antiphase r = {r:.3f}")

    res = dict(protocol="tonic E 2 / F 3 nA + kickoff, no contact",
               duration_s=DURATION, dt_s=DT, skip_s=SKIP,
               period_L_s=perL, period_R_s=perR, bursts_L=int(nL),
               bursts_R=int(nR), antiphase_r=r,
               lE_max_mV=float(lE.max()), rE_max_mV=float(rE.max()),
               neurons=len(net.idx), synapses=int(net.n_synapses))
    (HERE / "w2l_numpy_ref.json").write_text(
        json.dumps(res, indent=2), encoding="utf-8")
    np.savez_compressed(HERE / "w2l_numpy_ref.npz", V=V_all,
                        t=np.arange(nsteps) * DT,
                        names=np.array(sorted(net.idx, key=lambda n: net.idx[n])))
    print(f"saved {HERE / 'w2l_numpy_ref.json'} + .npz")

    if nL < 4 or nR < 4:
        print(f"W2L_NUMPY_REF FAIL too few RG-E bursts (L={nL}, R={nR})")
        return 1
    if not (r < -0.5):
        print(f"W2L_NUMPY_REF FAIL antiphase r={r:.3f} not < -0.5")
        return 1
    print(f"W2L_NUMPY_REF PASS period_L={perL:.3f} period_R={perR:.3f} "
          f"r={r:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
