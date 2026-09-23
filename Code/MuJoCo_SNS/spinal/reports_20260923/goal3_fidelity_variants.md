# Goal-3 fidelity variants - measured (2026-09-23)

## (A) Error-controlled substepping (tol sweep)

| tol (deg/sample) | refined frac | RMS vs 0.5ms (deg) | max vs 0.5ms (deg) | wall (s) |
|---|---|---|---|---|
| 0.01 | 5.0% | 0.467 | 12.753 | 3.6 |
| 0.05 | 0.1% | 0.312 | 8.338 | 3.3 |
| 0.2 | 0.0% | 0.251 | 6.797 | 3.2 |

(fixed 2 ms wall 1.5 s, fixed 0.5 ms wall 4.4 s; fixed 2ms-vs-0.5ms RMS 0.439 deg / max 11.897 deg)

## (B) Contact damping variants (impedance framework)

| variant | RMS vs baseline (deg) | max vs baseline (deg) | peak foot force (N) | max pen (mm) |
|---|---|---|---|---|
| baseline solref=[0.02, 1], solimp mid 0.95 | 0 | 0 | 1001 | -17.07 |
| less viscous: tc 0.01, zeta 0.7 | 3.203 | 46.006 | 3237 | -17.07 |
| non-linear depth: solimp mid 0.90, width 0.003 | 5.680 | 60.927 | 882 | -17.07 |

Naive direct-(k, b) via negative solref was also attempted and is UNSTABLE at the 2 ms timestep (NaN by t=0.012 s; it bypasses the solimp impedance scaling): a concrete demonstration of the goal-3 report's 'partially native' verdict. Exact Hunt-Crossley fitting and any production adoption remain open (Ben's ruling).

(C) SEE not re-implemented - see module docstring.
