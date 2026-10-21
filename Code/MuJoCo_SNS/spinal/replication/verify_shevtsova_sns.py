"""TEMPLATE-MODEL VERIFICATION: Shevtsova et al. 2026 lumbar circuit in
SNS-TOOLBOX (Ben directive 2026-10-20), LUMBAR-ONLY biped reduction.

Same circuit, same claims, same gates as verify_shevtsova.py (the
standalone numpy ODE that PASSES); the circuit under test is
shevtsova_sns.py - the 22-population build in sns_toolbox primitives
(Network / NonSpikingNeuron / NonSpikingNeuronWithPersistentSodiumChannel
/ NonSpikingSynapse + the repo backend patch family), driven by input
ports for alpha, integrated at the paper's dt = 1 ms. Mapping deviation
flags S-1..S-5 are documented in shevtsova_sns.py's docstring; the
headline one (S-1): the repo's fixed-tau NaP reduction holds her trot
band at alpha=0.55 (tau 0.30 s -> 4.45 Hz vs the standalone's 4.51) but
LATCHES at the co-active fixed point (~-42 mV both RG centers) at
alpha=1.05 for every fixed tau inside her stated 150-400 ms band, so the
default build restores her exact eq-8 tau_h bell through the same
one-line backend-patch mechanism (her published constants, no free
knobs).

Protocol: the paper's Fig-9-style stepped-alpha sweep (identical to
verify_shevtsova.py): deterministic (sigma 0), alpha in {0.3, 0.55, 0.7,
0.85, 1.0, 1.05}, 20 s each, first 5 s dropped, fresh initial condition
per alpha. Analysis = the paper's own definitions (midextension periods,
f(V_F)>=0.1 flexion criterion, normalized LR phase), imported from
shevtsova_standalone so the two verifiers measure identically.

LITERATURE CLAIMS TESTED (same quotes/gates as verify_shevtsova.py):
  C1  Trot at low alpha (alpha=0.55): both sides oscillating, freq in
      [3,6] Hz, LR phase median in [0.40,0.60].
  C2  Frequency rises with alpha; freq(1.05) >= 6.5 Hz (paper's sync
      onset ">6.5 Hz"); extensor fraction decreases 0.55 -> 1.05.
  C3  Hind alternation -> synchronization shift is a FULL-quadruped
      behavior (LPN pathways excluded from this reduction by design) -
      phase trajectory REPORTED, out of scope (rules _comment).

Verdict: SHEVTSOVA_SNS PASS iff C1 and C2 pass.

Re-run:  D:/Anaconda/envs/myo/python.exe verify_shevtsova_sns.py
         (optional --tau=fixed for the S-1 fixed-tau comparison build)
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                     # noqa: E402

from shevtsova_sns import build_circuit, simulate  # noqa: E402
from shevtsova_standalone import analyze           # noqa: E402

ALPHAS = (0.3, 0.55, 0.7, 0.85, 1.0, 1.05)


def main() -> int:
    tau_mode = "fixed" if "--tau=fixed" in sys.argv else "bell"
    c = build_circuit(tau_mode=tau_mode)
    print(f"build: {c['n_pop']} populations, {c['n_syn']} synapses, "
          f"{len(c['ports'])} alpha-drive input ports, dt 1 ms, "
          f"tau_mode={tau_mode}"
          + (" (fixed stand-in 0.300 s inside her stated 150-400 ms band)"
             if tau_mode == "fixed" else
             " (her eq-8 bell; flag S-1 in shevtsova_sns.py)"))
    rows = []
    for a in ALPHAS:
        t, out, _ = simulate(a, t_end_s=20.0, circuit=c)
        r = analyze(t, out, c)
        rows.append((a, r))
        print(f"alpha={a:4.2f}: freq L {r['l']['freq']:5.2f} R {r['r']['freq']:5.2f} Hz | "
              f"LR phase {r['phase_med']:5.3f} (n={r['phase_n']:3d}) | "
              f"ext_frac L {r['l']['ext_frac']:.2f} R {r['r']['ext_frac']:.2f} | "
              f"ext episodes {r['l']['n_ext']}/{r['r']['n_ext']}")

    by = {a: r for a, r in rows}
    r_lo, r_hi = by[0.55], by[1.05]

    c1_osc = np.isfinite(r_lo["l"]["freq"]) and np.isfinite(r_lo["r"]["freq"])
    c1_freq = c1_osc and 3.0 <= min(r_lo["l"]["freq"], r_lo["r"]["freq"]) \
        and max(r_lo["l"]["freq"], r_lo["r"]["freq"]) <= 6.0
    c1_ph = 0.40 <= r_lo["phase_med"] <= 0.60
    c2_rise = (c1_osc and np.isfinite(r_hi["l"]["freq"])
               and r_hi["l"]["freq"] > r_lo["l"]["freq"])
    c2_band = c2_rise and min(r_hi["l"]["freq"], r_hi["r"]["freq"]) >= 6.5
    c2_ext = (np.isfinite(r_lo["l"]["ext_frac"])
              and np.isfinite(r_hi["l"]["ext_frac"])
              and r_hi["l"]["ext_frac"] < r_lo["l"]["ext_frac"])

    checks = [
        ("C1a both sides oscillate at alpha=0.55", c1_osc),
        (f"C1b trot band 3-6 Hz at 0.55 "
         f"(got L {r_lo['l']['freq']:.2f} / R {r_lo['r']['freq']:.2f})", c1_freq),
        (f"C1c hind LR alternation ~0.5 at 0.55 (got {r_lo['phase_med']:.3f})",
         c1_ph),
        (f"C2a frequency rises 0.55->1.05 "
         f"({r_lo['l']['freq']:.2f} -> {r_hi['l']['freq']:.2f} Hz)", c2_rise),
        (f"C2b freq(1.05) >= 6.5 Hz (paper's sync onset; "
         f"got {min(r_hi['l']['freq'], r_hi['r']['freq']):.2f})", c2_band),
        (f"C2c extensor fraction shortens "
         f"({r_lo['l']['ext_frac']:.2f} -> {r_hi['l']['ext_frac']:.2f})",
         c2_ext),
    ]
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    print(f"  [SCOPE] C3 LR phase trajectory over alpha: "
          + ", ".join(f"{a:.2f}->{by[a]['phase_med']:.3f}" for a in ALPHAS)
          + " - full trot->bound transition needs the quadruped LPN "
            "pathways excluded from this reduction by design (documented)")
    verdict = all(ok for _, ok in checks[:4]) and c2_band and c2_ext
    print("SHEVTSOVA_SNS " + ("PASS" if verdict else "FAIL")
          + " (C1+C2; C3 recorded as out-of-scope for the lumbar-only "
            "reduction)")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
