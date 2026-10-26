# Feedback ablation pair — walking W2L model, frozen tuning (2026-10-26)

Ben's order: ONE saved fixed-tuning comparison, feedback intact vs
ablated, on the walking BPA-routed model. Pre-registered predictions in
PREREG_REPRO_ABLATION_20261026.md (written before the runs).

Config (both conditions, frozen): Ben's literal corrected-aproj defaults —
te 3.0 / tf 4.0, comm 1.0, protocol knobs at the campaign defaults
(legcoll 0, sink 0.06, rig 0.25, walker_yield 1.5), net = bilateralrg2,
30 s fresh-process evals. The ONLY difference: the ablated condition
zeroes every harness-exposed sensory channel — the heel cascade
(ch_hc 0.08→0), ankle contact (ch_ank 1.0→0), the crossed
heel→contra-RG-flexor edge (ch_contra 0.08→0), the contact encoder
(ctn_amp 4→0), and the Ib force-sensor current (ibnA 1.0→0). HONEST
SCOPE: the Li-style Ia/II length/velocity chains are not knobbed in this
harness and remain at template defaults in BOTH conditions; this is a
contact+force-feedback ablation, not a total deafferentation. (The
model's limbs are legs; "arm" appears below only as the JSON field name,
a data key, not prose.)

## Result (results\ablation_{intact,ablated}_20261026.json + .npz)

| measure | INTACT | ABLATED |
|---|---|---|
| survive | true | **false** |
| objective | 90.0 | **5.8** |
| cycles L/R | 33 / 28 | **1 / 0** |
| inter-leg phase | 0.492 (n=30) | undefined (no cycles) |
| forward progress | +1.163 m | **+0.091 m** |
| pelvis min | 0.905 m | 0.903 m |

The model does not merely lose phase lock — stepping itself collapses to
a single left cycle in 30 s when contact and force feedback are removed,
at identical drive and central conductances.

## Prediction scorecard

- P4 intact survives: CORRECT (matches the documented contact2 datum).
- P5 "phase degrades first, cycling continues": **WRONG in the
  informative direction** — cycling itself collapses; the sensory
  channels are load-bearing for pattern generation here, not merely
  phase trimming.
- P6 progress strictly decreases: CORRECT (−92%).

## Reading for the chapter

Same controller, same muscles, same drive: the sensory channels are the
difference between walking (obj 90) and near-stillness (obj 5.8). This
is the Deng-2023 active/disabled protocol shape on our robot model, and
the null side is documented separately on easteregg2 (per-path
remove-and-retest with exactly-zero deltas — WHICH feedback matters is
per-path, but the class as a whole is essential).

## Saved evidence locations

- Metrics: `results\ablation_{intact,ablated}_20261026.json` (this folder).
- Trial rows (append-only, evaluator's CSV):
  `spinal\reports_20261016\w2l_wean\wean_trials.csv`, tag
  `ablation_{intact,ablated}_20261026` (knobs + all metrics per row).
- Caveat: the run-trajectory npz did NOT persist (the harness's AARL_NPZ
  redirect did not take effect from this driver) — the saved comparison
  is the metrics jsons + CSV rows; re-run with the redirect fixed if raw
  traces are wanted for a figure.
