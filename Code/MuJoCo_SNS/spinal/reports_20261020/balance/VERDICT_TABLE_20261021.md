# Standing-balance + perturbation verdict table (2026-10-21, final)

All numbers first-hand from the artifacts in `spinal\` and
`spinal\reports_20261020\balance\` (sources in brackets). "Rig share" =
share of the ~706 N model weight NOT carried by the feet, measured from
contact forces post-warmup (`_s3_rigshare_20261021.json`,
`_s3p_syn6_rigshare_20261021.json`). Perturbation battery = 7 evals
(quiet + 4 pelvis pushes 60 N/150 ms at t=4 s + platform translate 6 cm
+ platform tilt 6 deg), DRIVE 0, 8 s (`_curriculum_s3p.py`).

| model | verdict | quiet stand | perturbation battery (games score) | rig share (load on feet) | Ib ablation |
|---|---|---|---|---|---|
| HUMAN main (s3k stock) | LEANS ON RIG | 40/40 trials stand; winner = incumbent (vest 0, rig 1.0), score 17.2945, fresh BIT-EXACT [curr_s3_balance; `_s3_verify_stock`] | 20 trials (19/20 zero-fall); winner = trial 0 incumbent + vest_prop 0.375, score −46.08, 0 falls/7, mean push sway 9.0 cm, tilt 28.8 deg; fresh BIT-EXACT (−46.0837) [curr_s3p_perturb; `_s3p_verify_stock_20261021.json`] | rig carries 49.8% (354 N on feet), tilt max 28.9 deg, contact duty 1.0/1.0 | no fall, but R knee crumples into 49.9 deg SUSTAINED flexion (intact −9.3 deg); foot load 50%→36% — Ib load feedback is what holds knee extension [`_s3_ibablate_stock`; rigshare] |
| w2lvar | LEANS ON RIG | 48 trials (47 stand); winner trial 35 score 25.6134 DB, fresh 25.3098 (+1.2%, fresh-process deterministic) [curr_w2lvar_s3_balance; `_s3_verify_w2l_bal`] | 18 trials (18/18 zero-fall); winner = trial 0 = the FIXED incumbent (vest 0, rig 1.0) at −17.49 — the vest-tuned quiet winner did NOT transfer to perturbations; fresh BIT-EXACT (−17.4920) [curr_w2lvar_s3p_perturb; `_s3p_verify_w2lvar_20261021.json`] | rig carries 52.8% (333 N on feet), tilt 25.3 deg [rigshare] | barely changes (tilt 25.3→19.7 deg, foot load 47→42%) [`_s3_ibablate_w2l_bal`] |
| syn6 | LEANS ON RIG (quiet); s3p "wean" pose = SAGS ON RIG | 48 trials (43 stand); winner trial 11 score 37.4630, fresh BIT-EXACT [curr_syn6_s3_balance; `_s3_verify_syn6_bal`] | 18 trials (17/18 zero-fall); winner trial 15 (−54.05 vs incumbent −68.42) adopted vest_prop 0.89 + stance width 9.2 deg + knee bias 4.7 deg + coact 0.23 + rig_scale 0.054, 0 falls/7, push sway 7.7 cm, tilt 32.1 deg; fresh BIT-EXACT (−54.0504) — BUT rig-share of THAT pose: feet carry only 29.8%, COM sags to 0.575 m; the low sway is a hang on the slack rig, NOT a free stand [curr_syn6_s3p_perturb; `_s3p_verify_syn6_20261021.json`; `_s3p_syn6_rigshare_20261021.json`] | quiet winner: rig carries 51.7% (341 N) [`_s3_rigshare`] | EXACTLY identical to intact (every metric) — the winner's extensor-Ib gains are all default 0; the syn6 quiet stand never used Ib feedback at all [`_s3_ibablate_syn6_bal`] |
| W2L robot (Hill, afferented) | PARTIAL: rigged stand only | rigged (0.25–1.0) stands all scales, sway 0.1–1.0 cm, tilt ≤21.4 deg, harness ~49–50% weight, heel duty 0.00; UNRIGGED step regime FALLS t=1.18 s; latch regime "survives" rig-0 only by hanging (sway 48 cm range, harness 50%) [w2l_mujoco\balance_w2l_20261020.log; balance_w2l_ibablate_20261021.log] | push battery at rig 0.25: 0 falls, post-event sway 2.2–4.8 cm (harness-held) [same logs] | harness carries ~50% | ibnA=0: rigged stand barely changes; unrigged fall 1.18→1.88 s (falls either way) [`balance_w2l_ibablate_20261021.log`] |
| BPA robot | LEANS ON RIG (structural) | retained-harness quiet stand (support 0.70 of 411 N): sway 2.1/0.4 cm; NO free-stand claim — runner states COM sits outside the support polygon, this body cannot statically stand [w2l_mujoco\balance_w2l_20261020.log] | 6 perturbed cases at support 0.70: post-event sway 32.5–43.7 cm, mean foot load 10–11% BW (passes the ≥10% gate); support 0.85 FAILS the gate (9%) [same log] | harness carries ~90% (foot load 10–11%) | N/A — tonic open-loop BPA drive, no afferents exist to ablate (stated, not dodged) |

## Plain-language summary

1. NO model in this program free-stands. Every human-family winner
   carries roughly HALF its weight on the support rig (48–52%); the W2L
   robot and BPA walker are harness/held (~50% and ~90%). Falls outside
   the rig: W2L steps-and-falls at 1.18 s; BPA body cannot statically
   stand at all.
2. Perturbations do NOT topple any rigged config: zero falls in every
   winner battery (60 N pushes, 6 cm translate, 6 deg tilt); sway stays
   2–9 cm for the human family. All three s3p winners reproduce their
   DB scores BIT-EXACT in fresh processes.
3. The searched balance mechanisms mostly did NOT beat the incumbents:
   stock and w2lvar winners are the incumbent pose (trial 0). Only syn6
   adopted wide stance + knee bias + co-contraction + strong vestibular
   — and its apparent rig-wean (rig_scale 0.054) is a measurement
   artifact of pose: the model sags (COM 0.575 m) and hangs on the
   slack spring with feet carrying only 30% of weight. The s3p score
   lacks a rig-share term; that is the honest lesson for any rework.
4. Ib (Golgi tendon) load feedback matters where it exists: ablating it
   crumples the stock model's right knee (50 deg flexion) without a
   fall; w2lvar shrugs; syn6 has no Ib in use at its winner at all; the
   BPA walker has no afferents (N/A).

Gates: census defaults re-run PASS this session (410/376/1186 stock,
888/382/7866 w2lvar, 794/382/3033 syn6; coact adds exactly 4 edges and
restores). Anchor −161.56754173676563 verified bit-exact by the M2
gatekeeper 2026-10-21 pass D after these edits. Fresh verifies:
quiet winners (stock 17.2945 bit-exact, syn6 37.4630 bit-exact, w2lvar
25.3098 vs 25.6134 +1.2% fresh-process deterministic) and s3p winners
(all three bit-exact).
