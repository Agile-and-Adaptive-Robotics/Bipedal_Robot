"""TEMPLATE-MODEL VERIFICATION: Deng 2019 two-layer CPG (Nourse 2023
Fig 6A, Tables A4-A7) - the literature-direct circuit of the replication
program (Ben 2026-10-20: "run the template models... until you can repeat
the results from literature").

WHAT THE LITERATURE/RECORD CLAIMS (replication targets):
  1. ONE 10 nA / 1 ms pulse into RG_HC_ext kicks the circuit into a
     SELF-SUSTAINED alternation that persists with NO drive thereafter
     (Deng 2019 latch/oscillation regime; Nourse 2023 Fig 6A).
  2. Cross-platform period agreement: Simulink 1.938 s vs numpy 1.94 s
     (DESIGN.md 2026-09-14 record, deng_cpg_ode.py vs
     SNS_Deng_Library.slx) - tolerance here +-2%.
  3. RG ext/flx traces anticorrelated, Pearson r about -0.86 (record).
  4. The toolbox voltage-dependent tau_h VARIANT quenches (the tau_h
     porting rule; checked as a negative control in the 20,000 s main,
     not here).

Gate: DENG_TEMPLATE PASS|FAIL with the three numbers printed.
"""
from __future__ import annotations

import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
SPINAL = os.path.dirname(HERE)
sys.path.insert(0, SPINAL)

import numpy as np                     # noqa: E402
from scipy.signal import find_peaks    # noqa: E402

from deng_cpg_ode import NAMES, simulate  # noqa: E402

T_END_MS = 60_000.0   # ms = 60 s (>> the ~2 s cycle: persistence is the claim)
DT_MS = 0.1


def main() -> int:
    out, _ = simulate("fixed", t_end=T_END_MS, dt=DT_MS)
    t = np.arange(out.shape[0]) * DT_MS / 1000.0
    ext = out[:, NAMES.index("RG_HC_ext")]
    flx = out[:, NAMES.index("RG_HC_flx")]
    w = t > 30.0                        # late window: self-sustainment
    ext_l, flx_l = ext[w], flx[w]

    pk, _ = find_peaks(ext_l, prominence=1.0)
    period_s = float(np.median(np.diff(t[w][pk]))) * 1000.0 if len(pk) >= 3 \
        else float("nan")
    r = float(np.corrcoef(ext_l, flx_l)[0, 1]) if ext_l.std() > 1e-9 else 0.0
    swing = float(ext_l.max() - ext_l.min())

    ok_persist = swing > 5.0 and len(pk) >= 3     # mV-scale late oscillation
    ok_period = abs(period_s - 1940.0) / 1940.0 < 0.02 if ok_persist else False
    ok_anti = r < -0.70                            # record -0.86, margin

    print(f"late window [30,60] s: HC_ext swing "
          f"{swing:.2f} mV, {len(pk)} peaks")
    print(f"period = {period_s:.1f} ms (record: Simulink 1938 / numpy 1940)")
    print(f"RG ext/flx Pearson r = {r:.3f} (record -0.86)")
    checks = [("self-sustained (no drive after the single pulse)", ok_persist),
              ("period within 2% of the 1.94 s record", ok_period),
              ("ext/flx anticorrelation r < -0.70", ok_anti)]
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    verdict = all(ok for _, ok in checks)
    print("DENG_TEMPLATE " + ("PASS" if verdict else "FAIL"))
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
