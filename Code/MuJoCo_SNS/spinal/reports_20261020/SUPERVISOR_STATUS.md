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
- 2026-10-20 wave-2 supervisor pass A: four streams DELIVERED by agents,
  one still running. (1) Contact sensing: feet DO contact (heel strikes
  multi-kN, 13-16x saturation; crossed-edge weakness = pulse duty, not
  air); RA/SA + shear/total modes + tunable presyn threshold/saturation
  built default-off; agent gates PASS; default smoke re-verified this
  pass. Ben committed files (8a125e75). (2) Shevtsova SNS: PASS, matches
  standalone to 2 decimals; fixed-tau backend LATCHES at high alpha (her
  bell asymmetry required) - program-wide finding. (3) syn6 Ib+antagonist:
  per-muscle Ib->E-layers gap wired (ib_e_central) + Ia PF gate + S5
  lamination; censuses/anchor bit-exact; ON-probe +409 exact. (4)
  AnimatLab: no harness exists (Root free, no shelf); old wiring falls
  then steps on its side 63/63 bursts (Ben's observation reproduced);
  corrected model 0 bursts in AnimatLab vs free-running SNS twin -
  suspect V3 edges retyped 2.749 uS (27x weak class) + F1 RE->RE exc;
  tuned drive LATCHES (ratios port, numbers don't). Balance stream
  (agent 5) still running with perturbations/metabolism/Ib-ablation.
  PENDING for next passes: per-stream M1/M2 supervisor gates (spawn
  SUPERVISOR_AGENT.md verbatim), balance agent completion, then final
  wave report + automation retirement. The fixed-tau latch finding and
  the AnimatLab V3 retype should be added to DESIGN.md open issues.
- 2026-10-21 supervisor pass B: the balance agent DIED on the 5-hour
  usage limit mid-stream (limit has reset). Relaunched fresh with the
  full requirement set (perturbations, metabolic cost, Ib ablation,
  mechanism knobs, all 5 models, supervisor gates). Wave-2 streams 1-4
  remain as logged in pass A (delivered, agent gates PASS; formal M1/M2
  gatekeeper spawns still pending next pass).
- 2026-10-21 supervisor pass C: spawned the formal M2 gatekeeper over
  the four delivered streams (bit-identity re-verification, edge
  citations, originals-untouched check, DESIGN.md open-issue coverage
  for the fixed-tau latch + AnimatLab V3 findings). Verdict pending its
  completion. Balance agent (relaunched) still running.
- 2026-10-21 stage-3 balance + perturbations stream, entry 1 (resume
  session, plain language): surveyed what the dead agent left. GOOD
  NEWS: its two background studies SURVIVED and were still running
  (curr_w2lvar_s3p_perturb, curr_syn6_s3p_perturb; left running, not
  killed). Already on disk and now re-verified first-hand: (a) the
  main-model quiet study curr_s3_balance is COMPLETE, 40/40 trials,
  log curr_s3bal_20261020.log, best = trial 0 (score 17.29) = the
  incumbent with vestibular gains 0 and rig 1.0 - an honest negative:
  the searched vest gains found NOTHING better than the incumbent
  quiet stand (40 stand / 0 fall). (b) Fresh-process re-verification
  of all three quiet winners (verify_20261021.log): stock reproduces
  its DB score BIT-EXACT (17.294510739883755), syn6 BIT-EXACT
  (37.46299873272405), w2lvar 25.3098 vs DB 25.6134 (+1.2%, same
  fresh-process value the dead agent got - deterministic in a fresh
  process; the in-study number differs by the documented in-process
  attractor effect). (c) Ib ablation re-run: stock byte-identical to
  the dead agent's json. (d) W2L Hill robot: re-ran intact stands +
  ibnA=0 ablation + full rig ladder (balance_w2l_ibablate_20261021.log):
  unrigged it FALLS (t=1.18 s step regime; latch regime survives
  rig-0.0 only by hanging - sway 48 cm, tilt 4.8 deg, heel duty 0.00,
  harness ~50% weight); rigged (0.25-1.0) all PARTIAL with harness
  ~49-50% weight and heel duty 0.00; ibnA=0 changes the rigged stand
  barely at all and the unrigged fall time 1.18 -> 1.88 s (falls
  either way). (e) BPA walker battery (dead agent's log re-read):
  support 0.70 passes the foot-load gate (>=10% BW after event) on all
  7 battery cases, post-event sway 32-44 cm; support 0.85 fails it
  (9%); the runner itself states no free-stand claim (COM outside the
  support polygon); NO Ib exists to ablate on the BPA walker (tonic
  open-loop drive, no afferents) - stated as N/A, not dodged.
  (f) NEW rig-share analysis (_s3_rigshare_20261021.py/json): the
  three human-family quiet winners carry ~48-50% of body weight on
  the feet (333-354 N of ~706 N) - i.e. the support rig holds up the
  other HALF at their winning rig_scale 0.95-1.0. Tilt envelopes
  19-29 deg. Honest interim verdict: these are rig-assisted stands
  with large lean, NOT free balance. (g) syn6 Ib ablation is EXACTLY
  identical to intact (every metric) - its stage-3 winner has all
  extensor-Ib gains at default 0, so there is no Ib to remove: the
  syn6 quiet stand never used Ib feedback at all. stock Ib ablation:
  no fall, but the right knee drops into 42-50 deg SUSTAINED flexion
  post-warmup (intact holds ~+9..+10 deg near extension) - the
  Ib-removed model partially crumples into the rig without falling;
  w2lvar barely changes (tilt 25.3 -> 19.7 deg, foot load 47 -> 42%).
  Studies still finishing: syn6 s3p (17/18), w2lvar s3p (~11/18).
  Supervisor M1 gates pending until they land.
- 2026-10-21 supervisor pass D: M2 GATEKEEPER VERDICT SHIP for all four
  delivered streams (it re-ran every gate itself: smoke 2.997/-0.782,
  410/376/1186 + anchor bit-exact, syn6 794/382/3033, w2lvar
  888/382/7866, ON-probe +409 reproduced, both Shevtsova verifiers
  PASS). Originals verified untouched; Ben's commits 8a125e75 + 5f829950
  carry the streams. Gatekeeper added DESIGN.md open issues 6 (AnimatLab
  V3 retype) + 7 (fixed-tau latch at high alpha). Remaining for wave
  close: balance agent (running) + its gate, then final report +
  automation retirement.
- 2026-10-21 supervisor pass E: balance agent mid-work and producing
  (reports_20261020\balance\ already holds fresh-verify + Ib-ablation
  JSONs for w2lvar/syn6 winners and a pull experiment). Streams 1-4
  remain M2 SHIP. No action needed this pass; wave close awaits balance
  delivery + its M1 gate.
- 2026-10-21 supervisor pass F (batched firings): no change - streams
  1-4 M2 SHIP; balance agent still mid-run (artifacts accumulating in
  reports_20261020\balance\). No action required; wave close remains
  blocked on balance delivery + its M1 gate. Heartbeat only.
- 2026-10-21 supervisor pass G: balance stream progressing - newest
  artifacts curriculum_s3p_stock.json + curriculum_stage3.json written
  2026-10-21 00:33 (the main-model perturbation study is running /
  seeded). Streams 1-4 M2 SHIP. Wave close still blocked on balance.
- 2026-10-21 supervisor pass H: heartbeat - balance artifacts still
  accumulating (see balance\ timestamps); streams 1-4 M2 SHIP; no new
  action. Wave close remains blocked on balance delivery + gate.
- 2026-10-21 supervisor pass I: balance stream found STALLED (zero
  python processes, no log growth since 00:33, no completion
  notification - second agent death suspected). Fresh resume agent
  launched: finish/relaunch the main-model perturbation study, collect
  the per-model verdict table, run the M1 gate, append the final
  balance report. Streams 1-4 remain M2 SHIP.

## FINAL BALANCE REPORT (2026-10-21 resume session; M1 gate SHIP — wave-2 CLOSED)

Plain language, numbers first-hand. Good news first: the perturbation
studies did NOT need relaunching — the background processes survived
both agent deaths and finished on their own (s3p_stock EXITCODE 0 at
00:24, s3p_syn6 at 00:54, s3p_w2lvar at 01:41; DBs show 20/20, 18/18,
18/18 COMPLETE). The resume session re-verified everything firsthand.

PER-MODEL VERDICTS (full table + sources:
`reports_20261020\balance\VERDICT_TABLE_20261021.md` +
`verdict_table_20261021.json`):

1. HUMAN main (stock): LEANS ON RIG. Quiet winner = the incumbent
   (vestibular gains 0, rig 1.0), score 17.2945, fresh BIT-EXACT;
   40/40 trials stood. Perturbation battery (60 N pushes, 6 cm
   translate, 6 deg tilt): winner is again the incumbent + vest_prop
   0.375, score −46.08, ZERO falls in 7/7 evals, mean push sway 9.0 cm,
   tilt 28.8 deg — fresh BIT-EXACT. The rig carries 49.8% of body
   weight. Ib ablation: no fall, but the right knee crumples into
   50 deg sustained flexion (intact −9 deg) — Ib load feedback is what
   holds knee extension.
2. w2lvar: LEANS ON RIG. Quiet winner 25.61 (fresh 25.31, +1.2%
   documented fresh-process delta). Perturbations: 18/18 trials
   zero-fall; winner = the FIXED incumbent (vest 0, rig 1.0) at −17.49
   (fresh BIT-EXACT) — the vestibularly-tuned quiet winner did NOT
   transfer to perturbed standing. Rig carries 52.8%. Ib ablation
   barely changes anything (tilt 25.3→19.7 deg).
3. syn6: LEANS ON RIG. Quiet winner 37.4630 fresh BIT-EXACT. s3p
   search DID adopt mechanisms (vest 0.89 + stance width 9.2 deg +
   knee bias 4.7 deg + coact 0.23 + rig_scale 0.054) and improved
   −68.42→−54.05 with zero falls — fresh BIT-EXACT — BUT the rig-share
   probe (NEW this session, `_s3p_syn6_rigshare_20261021.json`) shows
   that pose carries only 29.8% of weight on its feet (COM sagged to
   0.575 m): the low sway is a hang on the slack rig, NOT a free
   stand. Honest lesson recorded: the s3p score lacks a rig-share
   term. Ib ablation exactly identical to intact — the syn6 winner
   never used Ib feedback (its extensor-Ib gains are all 0).
4. W2L robot (Hill): PARTIAL, rigged stand only. Rigged (0.25–1.0)
   stands with sway 0.1–1.0 cm and harness ~50% weight; UNRIGGED step
   regime FALLS at t=1.18 s; latch regime survives rig-0 only by
   hanging (sway 48 cm, heel duty 0). Push battery rigged: 0 falls,
   sway 2.2–4.8 cm. ibnA=0 barely changes the rigged stand; unrigged
   fall 1.18→1.88 s (falls either way).
5. BPA robot: LEANS ON RIG, structurally. Retained harness
   (support 0.70): all 6 perturbed cases pass the ≥10% BW foot-load
   gate (10–11%, sway 32.5–43.7 cm); support 0.85 fails it (9%). The
   runner itself states no free-stand claim (COM outside the support
   polygon). Ib ablation N/A — the BPA walker is tonic open-loop with
   no afferents.

GATES: M1 supervisor gate SHIP (per-item evidence in
`reports_20261020\balance\M1_GATE_20261021.md`; census scripts re-run
PASS this session — 410/376/1186, 888/382/7866, 794/382/3033, coact
adds exactly 4 edges and restores; all three s3p winner batteries
re-run in FRESH processes, each score delta +0.0000 vs DB). One
deviation noted in the gate record: no subagent facility exists in the
resume environment, so the gatekeeper checklist was executed in-session
with every verification command re-run firsthand. Anchor −161.56754173676563
stands bit-exact per the M2 gatekeeper (pass D, post-s3p-edits).

WAVE-2 STATEMENT: with this delivery the wave is FULLY GATED — streams
1-4 were already M2 SHIP (pass D) and stream 5 (balance) is now
delivered + M1 SHIP. The 'Walker wave-2 standing supervisor' automation
can be retired by the main session. Ben: no commits were made; the
balance artifacts to commit live in `spinal\` (study logs, winner
jsons, verify/rigshare/verdict JSONs, `reports_20261020\balance\`) —
note the tree ALSO carries your own uncommitted torque-paper +
BiPulley work (untouched by agents), and the dead agent's scratch
`_supervisor_gate_batch1_20261021.ts` (unfired workflow draft) is safe
to delete.
- 2026-10-21 ~10:30Z BALANCE-OBJECTIVE FIX + REAL-BUDGET RERUN (Ben's
  order, wave 3; this session): the s3p objective under-rewarded
  weight-bearing - the syn6 winner sagged onto the slack rig (feet
  29.8% BW, COM 0.575 m) and the incumbents lean ~50%. FIXED in
  _curriculum_s3p.py (change recorded in its docstring) + additive
  runner.py metrics: (a) RIG-SHARE term -250*max(0, rig_share-0.10)
  graded from feet>=90% BW, hard sentinel -149.5 when quiet rig share
  >50%; (b) METABOLIC term -0.5*(bal_act_sum + bal_cocon_act) = mean
  SUMMED muscle activation with the antagonist-overlap mass weighted
  extra (typ. act_sum ~56 over 92 actuators). Every 2026-10-20 metric
  kept; sentinel ladder now nan -200 < all-fall -150 < RIG-SAG -149.5
  < clip -149. New runner metrics bal_foot_load/bal_rig_share/
  bal_act_sum/bal_cocon_act are additive-only; ANCHOR RE-CHECKED after
  the edit: s3k flavor A BIT-EXACT -161.56754173676563. Recorded
  studies untouched - new names curr_s3q_perturb (optuna_walk.db) and
  curr_syn6_s3q_perturb (optuna_syn6.db), 60 trials each, launched
  background 10:26/10:31Z (logs s3q_stock_20261021.log,
  s3q_syn6_20261021.log; seeds = quiet incumbent defaults AND the
  2026-10-20 s3p winner params). Fresh verify script
  _s3q_verify_20261021.py includes the explicit sag check (foot-load
  fraction + COM height). Ben's fallback addendum logged: if 60 trials
  still cannot beat the hand-set seeds under the fixed objective, that
  IS the finding - the next round becomes a deliberate HAND-SET sweep
  (knobs/values proposed in the final entry below) until a pattern
  emerges to reformulate the search around.
  CORRECTION (same session, ~10:50Z): the first launch's foot-load
  metric was WRONG twice over - runner log_contact double-counts each
  foot contact (added to both the heel and toe accumulators; the
  rigshare scripts' 0.5 factor compensates) and the post-warmup 0.6 s
  window includes the release transient (raw sums 1341 N over [0.6,8]
  vs 708 N over [4,8] on the same stand). First 3 stock + 2 syn6
  trials scored with a negative rig share (zero rig penalty) - KILLED,
  the contaminated curr_s3q_* studies PURGED (they were this session's,
  minutes old; Ben's recorded studies untouched), bal_foot_load fixed
  to the rigshare convention of record (0.5 x logged sum, steady
  t >= 4.0 s window), anchor re-checked BIT-EXACT, both studies
  RELAUNCHED clean 10:47Z with the same seeds.
- 2026-10-21 stage-3 balance stream, reconciliation entry (the FIRST
  resume session, woken late - scopes only the DELIVERED s3p-era
  results; the s3q wave-3 rerun above belongs to the other live
  session): the FINAL BALANCE REPORT and M1 gate SHIP were delivered by
  the pass-I resume session; this session adds an INDEPENDENT
  re-verification of every number, run after Ben's commit 30c4261d
  "torque and balance" (which took the stream's artifacts - the
  "no commits" line in the report was true when written, now
  superseded). Re-run first-hand: (1) anchor -161.56754173676563
  BIT-EXACT (_anchor_20261021.log); (2) both census gates PASS
  (410/376/1186 + 888/382/7866 + 794/382/3033; coact adds exactly 4
  edges per family, restore clean); (3) quiet winners: stock + syn6
  bit-exact vs DB, w2lvar 25.3098 deterministic fresh-process value
  (documented +1.2% vs in-study); (4) stock Ib-ablation byte-identical
  to the original run; (5) ALL THREE s3p winner batteries re-run a
  THIRD time in fresh processes (s3p_verify_20261021.log, START
  10:15:08): stock -46.0837, syn6 -54.0504, w2lvar -17.4920 - every
  delta +0.0000 vs DB, falls 0/7 each. Three independent sessions now
  agree on every number. Housekeeping: duplicate scratch
  _s3p_verify_best_20261021.py (accidental second copy of the
  verifier) deleted; the unfired workflow draft
  _supervisor_gate_batch1_20261021.ts KEPT - it compiles clean and
  lets the MAIN session fire a formal INDEPENDENT M1 gate spawn
  (batch 1: stock/W2L/BPA) beyond the in-session gate record, whose
  deviation (no subagent facility in the resume environments) is
  documented in M1_GATE_20261021.md. The s3q foot-load convention this
  wave-3 session standardized (0.5 x logged sum, t>=4 s window) matches
  the _s3_rigshare_20261021.py convention this session used for the
  verdict table - one convention of record. s3p-era stream CLOSED;
  wave-3 s3q in flight under the live session.
- 2026-10-21 ~17:15Z BALANCE-OBJECTIVE RERUN, RESULTS (Ben's wave-3
  order; both studies EXITCODE 0, 60/60 trials each, fresh-verified
  BIT-EXACT). PLAIN LANGUAGE, numbers first-hand:

  STOCK (main human model, curr_s3q_perturb): the fixed objective is
  now honest and the answer is negative - NO trial in 60 beat even the
  clip floor. Distribution: 44 trials clipped at -149.0, 15 hit the
  new RIG-SAG sentinel -149.5, 1 all-fall, 0 above the clip. Both
  seeds (incumbent + prior s3p winner) scored the sentinel (feet ~46%
  BW). Best feet-load seen anywhere: 64.4% BW (t50); the recorded
  winner is trial 8 (first at the clip): feet 53.9% BW, COM 0.753 m,
  sway 0.118 m, 0 falls/7, fresh-process score -149.0000 = DB
  BIT-EXACT. Clean study: 0 degenerate batteries, 0 bounce artifacts.
  MEANING: with vestibular/stance/knee/coact/pelvis/rig-wean knobs,
  this model cannot stand above ~2/3 of its weight on its feet; every
  zero-fall pose needs the rig for ~35-50%. The objective now says so
  instead of rewarding the lean.

  SYN6 (curr_syn6_s3q_perturb): the fix KILLED the sag exploit (both
  seeds sentinel -149.5, feet 34%/27%) but exposed that syn6's stands
  are mostly DEGENERATE: 43/60 trials break (instability) before the
  4 s measurement window opens - all 6 perturbed evals dead early
  (push_sway -1.0), and ALL 34 clip trials are such early-breaks. The
  recorded "winner" (trial 5, fresh BIT-EXACT -149.0000) is one of
  them: quiet run ends at t=1.15 s, the rig-share term then reads a
  SINGLE pre-break foot-slam sample (foot "load" 278.6% BW -> negative
  rig share -> zero penalty; no fall because COM never had time to
  drop). It is NOT a meaningful pose - the winner slot is a lottery
  among tied degenerates. The best HONEST syn6 result remains a
  sentinel rig-hang (-149.5). Third gaming vector for the record:
  early-break runs garbage every windowed metric.

  OBJECTIVE-TERM VERDICT: rig-share graded term + >50% sentinel work
  as designed against sagging; metabolic term behaves (act_sum
  53-57 across trials, small spread); REMAINING HOLES, both fixable
  runner-side next round: (1) early-break runs need a t_end gate
  (quiet eval must reach ~7.5 s else sentinel - no stand happened);
  (2) mean foot force can exceed BW in oscillating runs - use
  weight-bearing DUTY (fraction of steady samples in 0.8-1.2 BW)
  instead of the mean, or cap the credit at 1.0.

  BEN'S FALLBACK (logged per his addendum; this is the finding, not a
  dead end): 60 real-budget trials under the fixed objective did NOT
  beat the hand-set incumbents - stock cannot exceed ~64% feet-load
  with the searched mechanisms, syn6 cannot complete the battery
  un-rigged. NEXT ROUND = deliberate HAND-SET GRID, then reformulate
  the search around whatever pattern it shows:
  (a) add the t_end >= 7.5 s gate + duty-based foot-load first
      (runner-side, ~2 lines, default-off safe);
  (b) STOCK grid, quiet evals only (~40 s each): rig_scale
      {0.6,0.7,0.8,0.9,1.0} x stance_width {0,3,6,9} x knee_bias
      {0,4,8} x pelvis_ty {0.90,0.92}, vest_prop fixed 0.375
      (120 runs, ~1.5 h) - map the feet-load ceiling surface; look
      for WHERE (if anywhere) feet-load > 70% BW with sway < 0.15 m;
  (c) SYN6: first just measure t_end across its stage-3 winner
      neighborhood at rig 0.9-1.0 (the early-break fragility is the
      real finding); hand-set only poses that survive 8 s;
  (d) if the ~65% feet-load plateau is structural everywhere, the
      honest conclusion for the dissertation is exactly that: these
      SNS mechanisms rig-assist but do not free-stand, and the next
      mechanism to add is ankle strategy / foot geometry / trunk
      bias, not more tuning of the current knobs.
  M1 GATE (this session, in-session execution - no spawn facility,
  same deviation as the morning gate): VERDICT SHIP. Checklist
  re-verified firsthand (record:
  reports_20261020\balance\M1_GATE_S3Q_20261021.md): no figures/no
  dissertation/no wiring touched (runner diff in Ben's 45caa365
  verified metrics-only); anchor BIT-EXACT twice; both s3q winners
  reproduce DB-vs-fresh +0.0000; failures reported as failures.
  Advisories: the syn6 winner json mechanically records the degenerate
  t5 (read the status note, not just the json); t_end gate +
  duty-based foot load queued for the next round. NOTE: Ben committed
  mid-flight (45caa365 carries the code edits); still uncommitted for
  him: the s3q study logs/npzs/dbs, winner+verify jsons, gate record,
  and this status file.
