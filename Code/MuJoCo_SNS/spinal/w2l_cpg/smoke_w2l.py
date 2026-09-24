"""Smoke test for the transcribed W2L bilateral contact-driven CPG.

Builds the network (build_w2l_net.build, numpy backend), then integrates
>= 20 s at the repo's dt = 2 ms with:
  (a) a CONSTANT MILD DRIVE  - tonic current on all four RG half-centers
      through the TONIC ports (stack DRIVE convention; F tonic > E tonic -
      the F tonic is the escape engine, equal tonics co-latch);
  (b) ALTERNATING L/R HEEL-CONTACT PULSE TRAINS at ~1 Hz through the
      contact-neuron ports, at the AnimatLab adapter gain C = 20 nA
      (SESSION_NOTES build_contact.pl).
plus the verbatim antiphase kickoff: BOTH Stimulus_1 (into L RG ext) and
Stimulus_2 (into R RG flx) are 10 nA for t in [0, 0.01 s] in the AnimatLab
.aproj.

VERDICT (exit 0 only on PASS):
  W2L_SMOKE PASS period=<s> antiphase_r=<-1..1> rg_e_bursts_L=<n> rg_e_bursts_R=<n>
when BOTH RG-E half-centers burst rhythmically (>= 4 bursts each) in
antiphase (Pearson r < -0.5). Otherwise "W2L_SMOKE FAIL <reason>", exit 1.

Burst/antiphase conventions follow spinal/check_rhythm.py (threshold at half
the window max; np.corrcoef of the two RG-E traces).

Usage:
    python smoke_w2l.py [--tonic-e N] [--tonic-f N] [--freq HZ] [--width S]
                        [--duration S] [--contact-current N] [--c1 N]
                        [--v3 N] [--contact-syn N] [--plot]
"""
from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, str(Path(__file__).parent))

from build_w2l_net import DT, build  # noqa: E402

T_END_DEFAULT = 24.0        # s, >= 20 per the smoke contract
ANALYSIS_SKIP = 2.0         # s, drop the kickoff transient from the verdict


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tonic-e", type=float, default=2.0, dest="tonic_e",
                    help="nA, constant mild drive on each RG EXTENSOR "
                         "half-center (stack DRIVE convention)")
    ap.add_argument("--tonic-f", type=float, default=3.0, dest="tonic_f",
                    help="nA, constant mild drive on each RG FLEXOR "
                         "half-center. Kept ABOVE tonic-e: the F tonic is "
                         "the escape engine that cuts the E plateau (tonic-f "
                         "2.5 left the E bursts too duty-wide for the "
                         "antiphase gate; equal tonics co-latch)")
    ap.add_argument("--freq", type=float, default=1.0,
                    help="Hz of the alternating heel-contact trains")
    ap.add_argument("--width", type=float, default=0.30,
                    help="s, active fraction-width of each heel pulse")
    ap.add_argument("--contact-current", type=float, default=20.0,
                    help="nA per active contact (AnimatLab adapter Gain C=20)")
    ap.add_argument("--contact-syn", type=float, default=None,
                    dest="contact_syn",
                    help="override contact-neuron OUTPUT synapse gain "
                         "(W2L_GAINS['contact']; the 20 nA adapter current "
                         "is separate)")
    ap.add_argument("--c1", type=float, default=None,
                    help="override c1 commissural inhibition gain "
                         "(W2L_GAINS['c1_inh'])")
    ap.add_argument("--v3", type=float, default=None,
                    help="override V3 commissural excitation gain "
                         "(W2L_GAINS['v3_weak'])")
    ap.add_argument("--duration", type=float, default=T_END_DEFAULT)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    gains = {}
    if args.contact_syn is not None:
        gains["contact"] = args.contact_syn
    if args.c1 is not None:
        gains["c1_inh"] = args.c1
    if args.v3 is not None:
        gains["v3_weak"] = args.v3
    net = build(gains=gains)
    print(f"neurons: {len(net.idx)}  synapses: {net.n_synapses}  "
          f"inputs: {len(net.inputs)}")

    dt = DT
    nsteps = int(round(args.duration / dt))
    u = net.make_inputs()
    i_s1 = net.input_index("Stimulus_1")
    i_s2 = net.input_index("Stimulus_2 (10nA)")
    i_tonic = {"E": [net.input_index("TONIC " + w) for w in
                     ("L RG ext", "R RG ext")],
               "F": [net.input_index("TONIC " + w) for w in
                     ("L RG flx", "R RG flx")]}
    i_heel = {"L": net.input_index("L heel contact"),
              "R": net.input_index("R heel contact")}

    watch = ["L RG ext", "L RG flx", "R RG ext", "R RG flx"]
    cols = {w: np.zeros(nsteps) for w in watch}

    kick_until = 0.01  # s - the .aproj stimulus window (verbatim)
    period_s = 1.0 / args.freq

    for k in range(nsteps):
        t = k * dt
        u[:] = 0.0
        # (a) constant mild drive (stack DRIVE convention: both half-centers
        # of a side carry tonic) + the verbatim 10 ms 10 nA kickoff pair
        for i in i_tonic["E"]:
            u[i] = args.tonic_e
        for i in i_tonic["F"]:
            u[i] = args.tonic_f
        u[i_s1] = 10.0 if t < kick_until else 0.0
        u[i_s2] = 10.0 if t < kick_until else 0.0
        # (b) alternating heel-contact trains at ~1 Hz, C = 20 nA
        phase = (t % period_s) / period_s
        u[i_heel["L"]] = args.contact_current if phase < args.width else 0.0
        u[i_heel["R"]] = args.contact_current if (phase >= 0.5 and
                                                  phase < 0.5 + args.width) \
            else 0.0
        V = net.step(u)
        if not np.isfinite(V).all():
            print(f"W2L_SMOKE FAIL non-finite membrane potential at "
                  f"t={t:.3f} s")
            return 1
        for w in watch:
            cols[w][k] = V[net.idx[w]]

    # ---------------- verdict analysis (check_rhythm.py conventions) --------
    win = np.arange(nsteps) * dt >= ANALYSIS_SKIP
    lE = cols["L RG ext"][win]
    rE = cols["R RG ext"][win]
    if lE.max() < 1e-6 or rE.max() < 1e-6:
        print("W2L_SMOKE FAIL an RG-E trace is flat "
              f"(L max {lE.max():.3g} mV, R max {rE.max():.3g} mV)")
        return 1

    def bursts(sig):
        on = sig > 0.5 * sig.max()
        starts = np.flatnonzero(on[1:] & ~on[:-1]) + 1
        return starts

    sL, sR = bursts(lE), bursts(rE)
    nL, nR = len(sL), len(sR)
    r = float(np.corrcoef(lE, rE)[0, 1])
    per_L = float(np.diff(sL).mean() * dt) if nL >= 2 else float("nan")
    per_R = float(np.diff(sR).mean() * dt) if nR >= 2 else float("nan")

    print(f"window {ANALYSIS_SKIP:.0f}-{args.duration:.0f} s | "
          f"L RG E max {lE.max():.2f} mV, R RG E max {rE.max():.2f} mV")
    print(f"bursts: L={nL} (period {per_L:.3f} s)  R={nR} "
          f"(period {per_R:.3f} s)")
    print(f"left-right RG-E correlation (want < -0.5 for antiphase): {r:.3f}")
    print(f"F-side check: L RG F max {cols['L RG flx'][win].max():.2f} mV, "
          f"R RG F max {cols['R RG flx'][win].max():.2f} mV")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(4, 1, figsize=(9, 6), sharex=True)
            t = np.arange(nsteps) * dt
            for a, w in zip(ax, watch):
                a.plot(t, cols[w], lw=0.8)
                a.set_ylabel(w, fontsize=7)
                a.grid(True, alpha=0.3)
            ax[-1].set_xlabel("t [s]")
            fig.suptitle("W2L smoke: tonic + alternating heel trains "
                         f"@{args.freq:g} Hz")
            fig.tight_layout()
            out = Path(__file__).parent / "smoke_w2l.png"
            fig.savefig(out, dpi=130)
            print(f"plot: {out}")
        except Exception as e:  # plotting must never flip the verdict
            print(f"(plot skipped: {e})")

    if nL < 4 or nR < 4:
        print(f"W2L_SMOKE FAIL too few RG-E bursts (L={nL}, R={nR}, want >= 4)")
        return 1
    if not (r < -0.5):
        print(f"W2L_SMOKE FAIL antiphase correlation r={r:.3f} not < -0.5")
        return 1

    print(f"W2L_SMOKE PASS period={per_L:.3f} antiphase_r={r:.3f} "
          f"rg_e_bursts_L={nL} rg_e_bursts_R={nR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
