# CAMPAIGN REPORT — bi-pulley escalation ladder, 2026-10-22

Bounded demonstration run of the improved, configuration-aware optimizer
(`tmp/paper_bridge/opt_improve_20261022/Opt_run_BiPulley2_20261022.m`, mode
`escalate`) over the doctrine configuration space of the pulley sweep
directive (`Code/Matlab/Mesh_Optimization/PULLEY_SWEEP_DIRECTIVE_20261021.md`,
§1 and §5): nBPA ≤ 3, tackle line parts G ≤ 3, actuator sizes re-solved per
configuration. Machine easteregg2, MATLAB R2025a `-batch`, real budgets
(surrogateopt 1000 + patternsearch 5000 evaluations per solve, parallel,
`parpool` 10 — recorded in every result mat and in
`Results/BiPulley_escalation_choices_20261022.csv` columns `budgetS`/`budgetP`).

- Launch: 2026-10-22 04:22:13 (`campaign_start_stamp.txt`), command
  `OPT_BPULLEY2_RUN=1 OPT_BPULLEY2_FULL=1 OPT_BPULLEY2_MODE=escalate matlab -batch "run('Opt_run_BiPulley2_20261022.m')"`
  → `bi2_fullrun_escalate_20261022.log`; finished 04:53:49 with
  `Opt_run_BiPulley2_20261022 done (mode escalate, smoke 0)` and
  `MATLAB_EXIT=0` (≈ 32 min wall: 20-anchor replay gate + parpool start +
  15 solves).
- Per-solve wall times (ladder-CSV stamps): 1.2–2.5 min for soleus, tib_ant,
  rect_fem, bifemsh; 1.7–6.8 min for med_gas — inside the observed 90 s–7 min
  envelope of the 2026-10-21 doctrine round at the same budgets
  (`tmp/paper_bridge/SURFACES_INDEX.md:138-141`).
- No Ben file edited. The only code change this campaign is EDIT-6 on the
  campaign's own dated driver (env overrides `OPT_BPULLEY2_RUN` /
  `OPT_BPULLEY2_FULL` for `RUN_BATCH`/`FULL_RUN`; dry-run default unchanged),
  `Opt_run_BiPulley2_20261022.m:72-82`.

## 1. Preconditions: the gate suite (all green, blocking list EMPTY)

Every gate from the previous phase was run 2026-10-22 (logs in
`tmp/paper_bridge/opt_improve_20261022/gaterun_20261022/`), each ending
`EXIT=0`:

| Gate | Log | Result |
|---|---|---|
| G1 regression replay (5 recorded mats, 0.000 tol) | `gate1_regression_replay.log` | 5/5 PASS, worst dTorque = 0 |
| G2 deterministic replay, all anchors | first 7 lines of `bi2_fullrun_escalate_20261022.log` | 20/20 anchor mats reproduced exactly |
| G2 escalation advisor unit tests (incl. slack-gaming blocker fix) | `gate2_advisor_test.log` | ALL PASS |
| G2 kinematics guards (hip envelope + soft corners) | `gate2_kin_guards.log` | ALL PASS |
| G2 live builds (mount deflection, bowden, materials) | `gate2_live_builds.log` | ALL CHECKS PASS |
| G2 path-artifact audit (opensim 4.4.1 pinned) | `gate2_path_audit.log` | done, 11/92 rows flagged |
| G2 sweep-cell + G3 ladder-walk smokes | `gate2_sweep_cell.log`, `gate3_ladder_walk.log` | EXIT=0 both |
| Flexor sweep33 smoke | `gate2_flexor_sweep33.log` | EXIT=0, sanity checks pass |

The gates' blocking list was **empty** — nothing was skipped or waived. The
two review blockers found in the plan review (slack-gamed torque test;
"run-to-run tolerance" regression phrasing) were fixed before the gates ran
(`Opt_run_BiPulley2_20261022.m` EDIT-2 deterministic 0.000-tolerance replay;
`escalation_advisor_20261022.m` taut-only torque test) and are covered by the
advisor gate above.

## 2. What ran

Per site (the five-site single-builder manifest of
`biPulleySpecsFromOpenSim.m:80`: med_gas_r, soleus_r, tib_ant_r, rect_fem_r,
bifemsh_r): ladder seeded at the 2026-09-26 baseline config (nBPA 1, G 1),
solve → diagnose with `escalation_advisor_20261022` (torque test over TAUT
feasible cells only; slack/infeasible cells map to the RoM side) → step
(torque-short → nBPA+1, RoM/slack-short → G+1, both → both, capped at the
directive grid {3}×{3}) → re-solve sizes at the new rung. Gain doctrine of
record: insertion force gain = nBPA/G, travel gain = G, per-crossing gains
`[G 1]`. 15 solves total, one dated full-workspace mat per rung in
`Code/Matlab/Mesh_Optimization/Results/` (`BiPulley_opensim_<site>_ladder_rung<r>_nB<b>_G<g>_<stamp>.mat`,
34–47 KB each).

Anti-gaming integrity: the objective/constraint pair is byte-identical to
`Opt_run_BiPulley.m:387-412` as carried in the copy
(`Opt_run_BiPulley2_20261022.m:608-658`: `-margin + 10*max(0,worstC)`,
strain-floor and zero-pulley-infeasible constraints; margin prices slack
cells as 0 torque — visible as `marginRaw = -1` rows below). Soft-corner
weights were NOT applied in this run, and none were needed: every
slack/infeasible cell that fired a step lies inside the graded requirement
envelope (the soft-corner machinery only tapers cells beyond task
envelopes); both raw and taut margins are reported.

## 3. Per-site result (chosen configuration + escalation rationale)

Torque convention: `minTaut`/`peak` are the vecnorm magnitudes of
`Torque_ins` over the active crossing, taut feasible cells only
(`escalation_advisor_20261022.m` iExtract). `tautMargin = minTaut/target − 1`.
The driver's objective-of-record margin (`marginRaw`, component-wise min over
ALL cells, slack priced as 0) is listed separately and is never quoted as a
physical capability.

| site | ladder (rungs → chosen) | rest/tendon (m) | target (N·m) = fmax×5 cm | minTaut (N·m) | tautMargin | peak (N·m) | taut cells |
|---|---|---|---|---|---|---|---|
| med_gas_r | (1,1)→(1,2)→(2,3)→**(3,3)** | 0.459 / 0.176 | 77.90 | 4.99 | −93.6% | 332.5 | 20/25 |
| soleus_r | (1,1)→(2,1)→**(3,1)** | 0.212 / 0.021 | 177.45 | 34.57 | −80.5% | 84.4 | 25/25 |
| tib_ant_r | (1,1)→**(2,1)** | 0.340 / 0.020 | 45.25 | 59.00 | **+30.4%** | 114.2 | 25/25 |
| rect_fem_r | (1,1)→(2,2)→**(3,2)** | 0.437 / 0.020 | 58.45 | 63.64 | **+8.9%** | 186.1 | 25/25 |
| bifemsh_r | (1,1)→(2,2)→**(3,3)** | 0.078 / 0.020 | 40.20 | 0.84 | −97.9% | 36.6 | 25/25 |

Rationale chains (from `Results/BiPulley_escalation_ladder_20261022.csv`,
campaign rows timestamp ≥ 20261022_042213, and the LADDER CHOICE lines in
`bi2_fullrun_escalate_20261022.log`):

- **med_gas_r** (knee×ankle, the only true two-joint grid): (1,1) direct
  drive had **0 taut cells** (15 slack + 10 infeasible of 25) — the torque
  test itself was not evaluable, so only the G-side step was justified →
  (1,2) still 10 slack and taut torque 3.29 vs 77.9 → both → (2,3) **G = 3
  closed the RoM** (0 slack, 0 infeasible) but pure torque-short (3.43) →
  (3,3) at the grid cap, still torque-short (minTaut 4.99, and 5 slack cells
  returned at the knee = −120° column; worstC +0.26, see caveat C5).
- **soleus_r**: never slack at any rung (short ankle route) — a pure torque
  ladder: minTaut 11.53 → 23.02 → 34.57 N·m = **×2.00, ×1.50 per added BPA —
  exactly the nBPA/G insertion-force doctrine at G = 1** (ratios 2/1, 3/2).
  Still −80.5% at the cap: the 177.45 N·m target (fmax 3549 N × 5 cm) needs
  ≈ 5.1× what three 20 mm BPAs deliver at the binding pose.
- **tib_ant_r**: pure torque ladder, (1,1) −35.2% → (2,1) **+30.4%, meets
  target, converged**. minTaut ratio 29.33→59.00 = ×2.01 (doctrine again).
- **rect_fem_r**: (1,1) 15 slack (deep-flexion hip poses over-feed travel) +
  torque-short → both → (2,2) RoM closed (G = 2 sufficed), torque −27.2% →
  (3,2) **+8.9%, meets target, converged**. minTaut ratio 42.55→63.64 =
  ×1.497 (doctrine 3/2).
- **bifemsh_r**: (1,1) 10 slack + 10 infeasible and torque-short → both →
  (2,2) still 10 slack, minTaut collapsed to 0.60 → both → (3,3) RoM closed
  (25/25 taut, worstC = 0, patternsearch exit 1) but minTaut 0.84 = −97.9%
  at the cap. Note the trade the ladder surfaced: rung 1 had minTaut 19.8
  with 10 dead cells; rung 3 has full coverage with almost no torque — sizes
  re-solved per rung, rest pinned at its 0.30·Lchar lower bound (see C6).

**Headline, stated plainly: 2 of 5 sites (tib_ant_r, rect_fem_r) meet the
fmax×5 cm heuristic target inside the directive grid; 3 of 5 (med_gas_r,
soleus_r, bifemsh_r) exhaust the nBPA = 3 cap still torque-short.** The
ladder's job was to find that out per site without sweeping blind, and the
verdict chains above are the deliverable. Whether those three sites need
> 3 parallel 20 mm BPAs, larger-diameter actuators, or a relaxed envelope is
Ben's call — the directive grid ends at 3.

## 4. RoM / peak table vs the requirement envelopes

Requirements provenance: torque target = fmax × 5 cm nominal arm (fmax =
max_isometric_force from `gait2392_robotbody.osim`, the builder's documented
heuristic — NOT per-joint human torque targets, see C2); RoM requirement =
the graded rom_table ranges in `actuator_map_27.json` (they match the audit
CSV's rom_lo/rom_hi exactly). Full table:
`Results/BiPulley_rom_peak_20261022.csv`.

| site | joint grid (req range, deg) | taut coverage | failing poses | joint taut span (deg) |
|---|---|---|---|---|
| med_gas_r | knee [−120, 10] × ankle [−50, 20] | 20/25 (80%) | slack at knee −120° (all 5 ankle angles) | knee [−88, 10]; ankle full |
| soleus_r | ankle [−50, 20] | 25/25 | none | full |
| tib_ant_r | ankle [−50, 20] | 25/25 | none | full |
| rect_fem_r | hip flex [−15, 110] | 25/25 | none | full |
| bifemsh_r | knee [−120, 10] | 25/25 | none | full |

## 5. Quotable numbers

- Doctrine gain law verified in-run: min taut torque multiplies by the
  nBPA₂/nBPA₁ force-gain ratio at fixed G — exact to within a percent on the
  clean sites (soleus 11.53→23.02→34.57 = ×2.00/×1.50; tib_ant ×2.01;
  rect_fem ×1.497 vs the predicted 2/1, 3/2), and ×1.455 vs 3/2 for med_gas
  (3.43→4.99 at G = 3, off-prediction because its equilibrium cells shift
  between rungs).
- Direct-drive (1,1) with re-solved sizes goes slack on the long/2-joint
  routes: med_gas 15 slack + 10 infeasible of 25, rect_fem 15 slack,
  bifemsh 10 + 10 (the proven fixed-size slack mechanism, re-confirmed with
  sizes free to move). The short ankle routes (soleus, tib_ant) never slack
  at G = 1.
- G = 3 closes the RoM slack everywhere it was stepped to: med_gas 25/25 at
  (2,3) and (3,3)-minus-the-−120°-column; bifemsh 25/25 at (3,3).
- 15 solves, all at budgets 1000/5000, in ≈ 32 min wall on 10 workers
  (1.2–6.8 min per solve).
- Deterministic replay gate reproduced 20/20 recorded anchor mats exactly
  (dMargin = dWorstC = 0) inside the campaign process before any solve.

## 6. Caveats — every one stated plainly

- **C1 Path-artifact flags (directive §6).** Audit of record:
  `gaterun_20261022/path_artifact_audit_gait2392_simbody_20261022_040554.csv`
  (opensim 4.4.1 pinned, model + map sha1 recorded). Flagged among the five
  sites: **med_gas_r** (knee DOF: conditional via point in range,
  activation angle −45.25°, moment-arm step 1) and **rect_fem_r** (knee DOF:
  activation angle −83.75°, step 1). bifemsh_r, soleus_r, tib_ant_r:
  unflagged. Mitigating fact: this campaign consumes NO osim moment arms —
  the requirement is the fmax×5 cm heuristic and the routes are the
  builder's own 3-point paths — so no flagged artifact enters these numbers
  numerically. The flags bite the moment per-joint human torque targets
  replace the heuristic; those two sites must be re-audited then. Also
  inherited: tib_ant_r's osim muscle row is named "Tibialis Posterior
  (duplicate name)" (Gait2392 naming artifact; the map resolves it to
  tib_ant_r — cosmetic, unflagged).
- **C2 Target provenance (coupling-envelope adjacent).** The "requirement"
  here is fmax × 5 cm — `Opt_run_BiPulley.m`'s documented heuristic, not an
  OpenSim per-joint torque/moment-arm curve. The three torque-short verdicts
  are short against that heuristic. Per-joint human targets (Routing27_v2
  CSVs) are the plan's flagged follow-on; the shortfall sites may look
  different (better or worse) under real targets.
- **C3 Hip coupling envelope is PENDING-BEN and NOT consumed here.** The
  coupling table ships 9/9 PENDING-BEN rows and the untrusted-banner rule
  stands (`gate2_kin_guards.log`). None of the five sites sweeps
  hip_adduction (secondary coordinates pinned at zero,
  `biPulleySpecsFromOpenSim.m:29`; rect_fem's grid is hip FLEXION with the
  knee rigid), so no number above depends on the pending table. Future hip
  sites must not trust their envelopes until Ben approves control points.
- **C4 Margin conventions.** `marginRaw` (the driver objective's margin)
  prices slack cells as zero torque, so any slack-bearing rung reads −1 —
  that is the anti-gaming behavior working, not a physical capability. All
  capability numbers quoted here are taut-cell magnitudes; both columns are
  in the CSVs.
- **C5 med_gas_r's chosen rung is not fully feasible.** At (3,3): 5 slack
  cells (the knee = −120° column), worst constraint +0.26 (strain-floor
  violation), patternsearch exit 0 (hit the 5000-eval limit without a
  feasible point). Its margin numbers come from the penalized region; treat
  med_gas's (3,3) row as "best found inside the grid", not a validated
  design. The other four chosen rungs are fully feasible (worstC = 0, slack
  0, infeasible 0; patternsearch exit 1).
- **C6 bifemsh_r sizes are bound-pinned.** At (3,3): rest = 77.7 mm = its
  0.30·Lchar lower bound, tendon = 20 mm = its lower bound. The solve is
  feasible but edge-of-bounds; the −97.9% minTaut is what remains there.
- **C7 Scope.** Ladder mode exercises `moving_via` routing only (the
  directive's §2 regime); the bowden branch and mount-deflection rescoring
  were gate-tested (`gate2_live_builds.log`) but not run in this campaign.
  The flexor-side sweep33 driver was smoke-gated only, not run at full
  budget — this campaign is the five bi-pulley sites.
- **C8 CSV hygiene.** `BiPulley_escalation_ladder_20261022.csv` also carries
  4 earlier SMOKE rows (timestamps 20261022_035431–040759, 60-eval budgets)
  from the gate phase; every campaign row is timestamp ≥ 20261022_042422.
  The two deliverable CSVs contain campaign rows only.
- **C9 Environment note.** The instructed pre-addpath fragment
  `Robot_Data/Functions/Functions/ModernRobotics` does not exist on
  easteregg2; ModernRobotics lives at `Code/Matlab/Functions/ModernRobotics`
  and is picked up by the driver's own `addpath(genpath(Code/Matlab))`. The
  campaign ran clean on the driver's path setup (the addpath warning in the
  log's first line is that dead fragment, harmless).
- **C10 Records of record untouched.** Anchor mats were only read (gate
  asserts equality); no `Opt_run*.m`/`BiPam_pulley.m`/`MonoPam_pulley.m`/
  `biPulleySpecsFromOpenSim.m` file was modified. The pre-existing modified
  files in git status (`Knee_Extensor_20mm.m`, `Results_table_20mm.xlsx`,
  `readserialnumbers2.m`, `Ext20mm_compare_20261022.png`, new ExtTest_3_*
  test data) predate this campaign and were not touched by it. New artifacts
  for Ben to commit: the 15 ladder mats + 2 deliverable CSVs in
  `Code/Matlab/Mesh_Optimization/Results/`, the appended ladder CSV, and the
  `tmp/paper_bridge/opt_improve_20261022/` campaign files (driver EDIT-6,
  logs, tables script, this report).

## 7. Artifacts

- `Code/Matlab/Mesh_Optimization/Results/BiPulley_escalation_choices_20261022.csv` — per-site chosen config, verdict chain, sizes, targets, margins, budgets, mat filename.
- `Code/Matlab/Mesh_Optimization/Results/BiPulley_rom_peak_20261022.csv` — RoM coverage + min/peak torque vs requirement envelopes.
- `Code/Matlab/Mesh_Optimization/Results/BiPulley_escalation_ladder_20261022.csv` — raw per-rung log (15 campaign rows + 4 smoke rows, see C8).
- `Code/Matlab/Mesh_Optimization/Results/BiPulley_opensim_<site>_ladder_rung<r>_nB<b>_G<g>_20261022_<stamp>.mat` × 15 — full-workspace result mats.
- `tmp/paper_bridge/opt_improve_20261022/bi2_fullrun_escalate_20261022.log` — complete campaign log (MATLAB_EXIT=0).
- `tmp/paper_bridge/opt_improve_20261022/campaign_tables_20261022.{m,log}` — table builder + its log.
