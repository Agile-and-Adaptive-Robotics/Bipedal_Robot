# Goal-5 pilot v2: production s3k config vs two references

Seed = curriculum_stage3.json (s3k trial 34) run through the exact stage-3
machinery (_curriculum.set_stage + runner --eval), with kine_ref.REF_CACHE
swapped per reference. kine lower = better. The study's non-sampled keys
(rg_adapt, kx, e2_pf/f1_df/f1_kf/e2_adapt/post_kneext/post_hipext) were
reconstructed NEUTRALLY (current defaults / factors 1.0) - the s3k-era
code's own defaults were not archived.

| drive | vs subject01 | vs Falisse | Falisse margin |
|---|---|---|---|
| 1.766 (x0.85) | -241.58 | -228.00 | +13.6 |
| 2.077 (x1.00) | -240.91 | -222.47 | +18.4 |
| 2.389 (x1.15) | -239.44 | -224.88 | +14.6 |

## Findings

1. The walker WALKS through this harness (scores far from the -25
   no-rhythm sentinel) - the reference-swap path is valid and internally
   consistent: identical config and code, only the reference differs.
2. **The current s3k gait matches Falisse's predicted walking BETTER than
   subject01 at every drive tested** (margin +13.6 to +18.4). With
   Falisse's faster, stiffer cycle (T 1.113 s, knee_min -59.2 deg, hip
   range 51.6 deg vs subject01's 1.233 s / -69.7 deg / 43.3 deg), this
   says the s3k-tuned walker is morphologically closer to a stiff
   predictive gait than to the experimental reference it was tuned
   against - the subject01 misfit concentrates exactly in the known gaps
   (knee flexion depth, hip range).
3. Drive response is weak at +-15%: subject01 is flat (-241.6 -> -239.4),
   Falisse peaks at the seed drive. No strong cadence-tension gradient at
   this scale; the tension signal lives in WHICH reference is fit better,
   not in the drive derivative.

## Caveats

- Absolute scores are NOT comparable to the recorded study value -159.45:
  that number is the curriculum study's combined objective (kine + kz /
  tilt terms), while these are raw kine_score from the same 16 s eval.
- Reproduction fidelity is limited by the neutral reconstruction of the
  non-sampled keys (above); the cross-reference COMPARISONS are unaffected
  (same reconstruction on both sides).
- Companion finding (pilot v1): the v10-era `--fitted --best10` winner NO
  LONGER WALKS under current physics (100% double-support, zero cycles at
  any drive 2.6-3.7; goal5_pilot_results.csv + %TEMP% npz diagnostics).
  All v8b/v9/v10 winners predate the RoM-limit / contact-surgery era -
  treat their recorded scores as historical until re-verified.

## How to reproduce

python gait_lib_pilot2.py  (myo env, ~14 min; writes the csv + this table)
