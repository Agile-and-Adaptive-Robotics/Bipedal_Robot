# Walker + replication campaign status (2026-10-20, easteregg2)

Maintained by the standing supervisor automation (hourly). Ben's
directives this campaign: (1) run the corrected W2L connectome in
MuJoCo+SNS and tune it, /supervisor-gate the results; (2) reflect the
contact pathways in the other MuJoCo models; (3) run the literature
template models under supervision until results repeat.

## Work streams and exit criteria

| stream | exit criterion | status |
|---|---|---|
| W2L corrected connectome (bilateralrg2) | optuna winner survives 30 s in a FRESH process + basin gate (±1%) + M1 supervisor gate | M1 SHIP (re-gate verified all fixes; gatekeeper confirmed every number by re-running gates). Tuned + fresh-verified (obj 91.4). Basin: NO 11/12 pass; ch_contra bracket complete (0.08/0.3/0.5/0.8 → 6/6/9/9 of 12 — plateau, structural marginality); accept-or-restructure = Ben's call (DESIGN.md #4). Ladder degenerate (documented, no weaning claim). |
| Contact pathways in other families | census+anchor gates PASS with knobs default-0 | DONE 2026-10-20 (410/376/1186 + anchor bit-exact; w2lvar 888/382/7866; syn6 794/382/3033; ON-probe +38 exact) |
| Deng 2019 template | DENG_TEMPLATE PASS | DONE 2026-10-20 (period 1937.8 ms vs 1938/1940 record; r −0.859 vs −0.86; self-sustained; re-verified by the M1 gatekeeper) |
| Shevtsova 2026 template | SHEVTSOVA_TEMPLATE verdict recorded (PASS, or documented FAIL) | PASS 2026-10-20 (C1+C2; C3 trot→bound scoped out for the lumbar-only reduction; k_tau calibrated stand-in flagged; re-verified by the M1 gatekeeper) |
| Rybak 2006a template | RYBAK2006_TEMPLATE verdict recorded | PASS 2026-10-20 (13/13: drive-asymmetry regimes + non-resetting/resetting deletion diagnostic; deterministic re-runs; independent session re-run exit 0; flags F-A..F-I recorded as DESIGN.md open issue 5; τ_h,max=1200 ms extracted from the local fulltext) |
| Supervisor gates | M1 (tuning results) + M3 (literature fidelity of template builds) run and acted on | BOTH SHIP 2026-10-20 (each after one BLOCK-fix-re-gate cycle; gatekeepers re-ran gates/verifiers themselves - all numbers reproduced first-hand, Rybak cross-process bit-identical) |

## Machine facts for any supervisor run

- python: `D:\Anaconda\envs\myo\python.exe` (bare `python` = Store stub).
- Campaign folder: `Code\MuJoCo_SNS\spinal\reports_20261020\w2l_contact2\`
  (campaign2.log, round_best_r*.json, README.md = the delta record +
  flagged items F1/F2 for Ben).
- Replication: `Code\MuJoCo_SNS\spinal\replication\` (verify_deng.py,
  shevtsova_standalone.py + verify_shevtsova.py, rybak2006_standalone.py
  + verify_rybak2006.py).
- Supervisor gate procedure: plugins\mujoco-sns-walker\commands\
  supervisor-gate.md — spawn SUPERVISOR_AGENT.md FRESH as a
  general-purpose agent at milestone M1 (tuning claims) / M3 (literature
  fidelity). Act on BLOCK verdicts.
- Basin gate on finalists: reports_20261016\w2l_wean\_w2l_wean_basin.py
  pattern (±1% joint perturbations, 12 trials, survival persists).
- Anchor of record for the human stack: s3k −161.56754173676563 +
  census (410,376,1186); w2lvar (888,382,7866); syn6 (794,382,3033).
- NO commits — Ben commits via GitHub Desktop; leave the tree clean and
  tell him what to commit.

## Log (append, newest last)

- 2026-10-20 ~04:00Z session start: bilateralrg2 template generated
  (96/204), W2L_SMOKE PASS r −0.782, MuJoCo runner --net=bilateralrg2,
  first datum: Ben's untuned defaults SURVIVE (obj 90.0, 19/17 cycles,
  phase 0.498). Family knobs contra_flex/contact_ank_mn gated + verified.
  DENG_TEMPLATE PASS. Shevtsova + Rybak agents building.
- 2026-10-20 04:2xZ: r100 round 1 done (best obj20 91.4); fresh-verify
  had omitted EV.BASE protocol knobs (rig/legcoll/sink/walker_yield) —
  FIXED, campaign2 relaunched (adds 3 trials, re-verifies, runs ladder).
- 2026-10-20 04:4xZ: campaign2 r100 winner (te 5.27, tf 7.55, ch_hc 0.48,
  ch_ank 0.91, ch_contra 1.10, ctn 4.61, cref 129.3, ibnA 0.17,
  comm 0.70) SURVIVES 30 s fresh-process; ladder descending (r085
  running, seed already survives at obj 90.0). Standing supervisor
  automation active (hourly).
- 2026-10-20 ~11:0x-11:3xZ supervisor pass (session): campaign2 done —
  all 7 rungs survive30(fresh)=True BUT the ladder is DEGENERATE
  (dscale bookkeeping only; rung verifies byte-identical obj 87.17257).
  Basin gates (12 fresh trials each, ±1%): r030 winner 7/12; r100
  winner 1/12; Ben defaults 6/12 (cycles robust, phase drifts);
  comm×2 4/12; comm×4 1/12 (over-coupling → near-sync); ch_contra=0.5
  9/12 (best; the crossed edge hardens the phase lock). Robust-objective
  study (34 trials) converged back to Ben's defaults (obj_mean 84.03).
  NO basin pass claimed. SHEVTSOVA_TEMPLATE built + PASS (C1+C2; C3
  scoped out; k_tau calibrated 6.67 mV flagged). M1 gate: gatekeeper
  verified all numbers by re-running gates, BLOCKED on missing negative-
  result reporting → 3 ordered doc fixes + advisories APPLIED (README
  basin/wean/CSV sections, DESIGN.md open issues, status rows); re-gate
  in flight. Rybak agent relaunched with lessons (first attempt
  degenerated, no files).
- 2026-10-20 supervisor pass 2: M1 re-gate VERDICT SHIP (all 3 fixes
  verified against raw rows; 1 nit fixed: CSV note timestamp range
  11:01-11:11). verify_deng + verify_shevtsova re-run: both PASS.
  Rybak files exist (rybak2006_standalone.py + verify_rybak2006.py),
  agent still iterating (quick probe did not complete in 240 s — left
  to the agent; verdict next pass). ch_contra bracket completed:
  0.3 -> 6/12, 0.8 -> 9/12 (with 0.08->6, 0.5->9 from pass 1): plateau
  at 0.5-0.8, residual marginality structural. Automation stays active
  (Rybak verdict + M3 gate + Ben's open decisions pending).
- 2026-10-20 supervisor pass 3: Rybak agent DELIVERED —
  RYBAK2006_TEMPLATE PASS 13/13 (drive-asymmetry regimes; deletion
  diagnostic non-resetting PF-level 0.001 T vs resetting RG-level
  0.416 T; deterministic; τ_h,max=1200 ms extracted from the local
  fulltext). Session independently re-ran verify_rybak2006.py: exit 0.
  M3 gate: BLOCK — Rybak flags absent from DESIGN.md open issues +
  stale "in progress" line + Shevtsova F-B..F-E coverage. Fixes
  applied (DESIGN.md items 3+5, verdict line, status row, scratch
  _shev_grep*.py removed); M3 re-gate in flight. If SHIP: final
  report + automation retirement next.

## FINAL CAMPAIGN REPORT (2026-10-20; M3 re-gate SHIP — automation retiring)

Ben's three directives, all executed and gated:

1. CORRECTED W2L CONNECTOME RUN + TUNED (M1 SHIP). Ben's 2026-10-18
   aproj transcribed GUID-level as the bilateralrg2 template (96/204;
   make_bilateralrg2.py; 7-point delta + flags F1/F2 in the campaign
   README). Runs in MuJoCo via test_w2l_ground.py --net=bilateralrg2.
   Results: (a) Ben's UNTUNED defaults survive the wean-campaign walk
   gate fresh-process — obj 90.0, 19/17 cycles, inter-leg phase 0.498;
   (b) tuned winner obj20 91.4, 30 s fresh-survival 31/28 cycles,
   phase 0.502, slope 0.049; (c) HONEST basin state: no 11/12 config
   (best 9/12 at ch_contra 0.5-0.8; bracket 6/6/9/9; residual
   marginality structural — phase excursions with stepping intact);
   (d) the campaign's wean ladder degenerated (dscale bookkeeping
   only) — documented, no weaning claim. The old bilateralrg template
   and every record artifact re-ran bit-identical.
2. CONTACT PATHWAYS REFLECTED IN ALL FAMILIES. params.G contra_flex
   (crossed heel -> contra RG-F, exc) + contact_ank_mn (heel+toe ->
   ankle_pf MNs) in s3k, syn6, w2lvar, spiking; default-0 =
   bit-identical: 410/376/1186 + s3k anchor -161.56754173676563
   BIT-EXACT, w2lvar 888/382/7866, syn6 794/382/3033, ON-probe +38
   exact, one-step finite everywhere.
3. TEMPLATE MODELS REPEATED (M3 SHIP). DENG_TEMPLATE PASS (1937.8 ms
   vs 1938/1940 cross-platform record, r -0.859, self-sustained).
   SHEVTSOVA_TEMPLATE PASS C1+C2 (4.51 Hz in the 3-6 Hz trot band,
   LR phase 0.500, freq -> 6.92 Hz >= 6.5, extensor fraction 0.42 ->
   0.22; trot->bound scoped out for the lumbar-only reduction).
   RYBAK2006_TEMPLATE PASS 13/13 (drive-asymmetry regimes; the
   deletion diagnostic: PF-level non-resetting 0.001 T vs RG-level
   resetting 0.416 T — the paper's two-level claim; deterministic;
   gatekeeper-verified cross-process; tau_h,max=1200 ms extracted).
   Gatekeepers personally re-ran all verifiers — every number
   reproduced first-hand.

AWAITING BEN (DESIGN.md 2026-10-20 OPEN ISSUES 1-5): F1 RE->RE sign
ruling; F2 weak-class strength doctrine; Shevtsova k_tau + F-B..F-E;
basin accept-9/12-or-restructure; Ryak stand-ins F-A..F-I. No
commits made — commit list in the campaign README + gatekeeper notes.
