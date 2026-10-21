"""TEMPLATE-MODEL VERIFICATION: Shevtsova et al. 2026 lumbar circuit
(eLife 14:RP107480), LUMBAR-ONLY biped reduction (see
shevtsova_standalone.py for the parameter provenance + flags F-A..F-E).

Protocol: the paper's Fig-9-style stepped-alpha sweep, deterministic
(sigma 0 = the paper's stable-trajectory runs at sigmaNoise 5 fA),
alpha in {0.3, 0.55, 0.7, 0.85, 1.0, 1.05}, 20 s each, first 5 s
dropped. Analysis = the paper's own definitions (midextension periods,
f(V_F)>=0.1 flexion criterion, normalized LR phase).

LITERATURE CLAIMS TESTED (quotes: Fig 3/4 text):
  C1  Trot at low alpha: "trot was usually expressed (alpha=0.55)";
      "trot ... expressed at oscillation frequencies from about 3 to
      6 Hz"; hind left-right "alternation (~0.5 phase difference)".
      GATE: at alpha=0.55 - both sides oscillating, freq in [3,6] Hz,
      LR phase median in [0.40,0.60].
  C2  Frequency rises with alpha "from ~1 to ~7 Hz" (Fig 9), "driven
      primarily by shortening of the extensor ... phases".
      GATE: freq(1.05) > freq(0.55); freq(1.05) >= 6.5 Hz (paper:
      synchronization seen ">6.5 Hz in model"); extensor fraction
      decreases from alpha 0.55 -> 1.05.
  C3  Hind left-right "gradually shifted from alternation (~0.5) to
      synchronization (~0/1)" with alpha. IN THE LUMBAR-ONLY REDUCTION
      this is EXPECTED NOT TO FULLY TRANSITION: the paper's gait
      changes arise in the FULL quadruped circuit (the diagonal and
      homolateral LPN pathways, dV0D/aV3/dSh2/dLPNi, are excluded from
      this reduction by design - rules _comment). GATE (scoped): the
      measured phase trajectory is REPORTED; the PASS criterion is only
      that the reduced circuit does not leave the alternation band it
      is entitled to hold (paper's own lumbar isolation behavior); the
      full transition is recorded as OUT OF SCOPE for the reduction.

Verdict: SHEVTSOVA_TEMPLATE PASS iff C1 and C2 pass.
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np                     # noqa: E402

from shevtsova_standalone import analyze, build_circuit, simulate

ALPHAS = (0.3, 0.55, 0.7, 0.85, 1.0, 1.05)


def main() -> int:
    c = build_circuit()
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
    print("SHEVTSOVA_TEMPLATE " + ("PASS" if verdict else "FAIL")
          + " (C1+C2; C3 recorded as out-of-scope for the lumbar-only "
            "reduction)")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
