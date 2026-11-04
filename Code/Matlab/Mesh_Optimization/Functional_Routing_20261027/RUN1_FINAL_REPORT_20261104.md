# Functional routing campaign — FINAL RUN REPORT (2026-11-04)

**Machine: easteregg2.** Every path below is on easteregg2 at
`D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\campaign_20261103\`.
From EB475WS4 or the laptop the same folder is reachable at
`\\easteregg2\D$\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\campaign_20261103\`.

Campaign: corrected 27-group BiPulley reoptimization per
`MASTER_DIRECTIVE_20261103.md` + `HARDWARE_AND_ZCODE_HANDOFF.md` + Ben's
2026-11-04 eight-item continuation directive and per-muscle request.
This report supersedes STAGE1_REPORT_20261104.md for results (that report
remains the stage-1/infrastructure record).

## What ran (of record)

- Sizing screen v2 (`sizing_search.log`, `sizing_winners.json`): all 27
  groups, real per-group evaluation, full diameter/count/family/gain grid,
  150-eval screen. 14 groups feasible at screen budget; 13 recorded
  best-infeasible with limiting constraints (all listed below); every
  group carries a retained best candidate into run_1's full budgets.
- run_1 four passes (81 checkpoints): true winners, rigid + compliant
  (Xi0–Xi3 per-pose) branches at 600 evals each, layered deep-outward with
  strictly-earlier-layer envelopes, refinement round, per-constraint
  honest feasibility labels.
- Verification: gate r14 PASSED (all checks; log `output/gate_r14.log`),
  dense 25-cell verification on all 54 (group, branch) pairs
  (`dense_verification.log`), per-muscle named-counterpart comparison
  (`per_muscle_compare.log`), signed human-vs-robot curve figures ×27
  (`figures/`), weight sensitivity PASS (`weight_sensitivity.log`).

## Final per-group status (both branches)

13 groups feasible on BOTH branches: bifemsh, ercspn, ext_dig, flex_dig,
grac, iliacus, intobl, lat_gas, med_gas, psoas, soleus, tib_ant, vas_int.

14 groups with at least one best-infeasible branch (best candidate
retained; limiting constraint of record from the checkpoint rep, in-layer):

| group | branch | violEff | limiting constraint (in-layer) |
|---|---|---|---|
| extobl | both | 1.67e-04 | centerline z=0: 0.27 mm (pelvis frame) |
| bifemlh | both | (in-layer) | retained rigid envelope (deeper fittings) |
| semimem | rigid | (in-layer) | retained rigid envelope |
| semimem | compliant | 5.71e-01 | tx slack cells: 57% of grid |
| sar | rigid | 1.00e+00 | tx infeasible cells: 100% of grid |
| sar | compliant | 7.99e-02 | span knee pass1: 3.1 mm |
| tfl | both | (in-layer) | retained rigid envelope: 4.3 mm (rig) / 8.6 mm (comp) |
| rect_fem | both | (in-layer) | retained rigid envelope: 23.5 mm |
| tib_post | both | 4.61e-03 | foot occupancy ban: 4.7 mm |
| tib_ant | compliant | 2.22e-01 | tx infeasible cells: 22% of grid |
| per_brev | rigid | 2.41e-02 | foot occupancy ban: 24.2 mm |
| per_brev | compliant | 6.91e-01 | tx slack cells: 56% of grid |
| per_long | both | 2.82e-02 | foot occupancy ban: 28.3 mm |
| per_tert | both | 6.6e-03 | foot occupancy ban: 6.7 mm |
| flex_hal | both | 3.7e-03 / 4.0e-04 | foot occupancy ban: 3.8/5.0 mm |
| ext_hal | both | 1.12e-02 | foot occupancy ban: 11.4 mm |
| glut_max3 | rigid | 1.11e-01 | tx infeasible cells: 11% of grid |
| glut_max3 | compliant | (in-layer) | retained rigid envelope: 20.7 mm |
| add_mag3 | both | 3.97e-02 | centerline z=0: 39.8 mm (documented seed crossing at hip adduction >= 11.9 deg) |

"violEff 0.00 / none" rows in `deliverable_status.csv` are context-free
re-evaluations; the IN-LAYER violator of record is quoted above from each
checkpoint's own solve report (the layer's obstacle context is what makes
them infeasible — evaluating in isolation removes it).

## Constraint families blocking groups (decision points for Ben)

1. **foot_occupancy_ban (7 lower-leg groups, 3.8–28.3 mm):** the
   sleeve-fitting tip must stay proximal of the ankle plane; these human
   counterparts attach distally. Resolving needs either a design change
   (longer tendon-only span rule relaxation is NOT allowed — the rule is
   your hard constraint) or smaller standoffs/longer sleeves; currently
   honestly infeasible.
2. **centerline z=0 (add_mag3 39.8 mm, extobl 0.27 mm):** add_mag3's
   femur-side seed row crosses the pelvis midline at hip adduction
   >= 11.9 deg (packaging agent documented). Needs an explicit
   on-centerline seed resolution ruling or a bounded patch that stays
   on-side at all poses.
3. **retained rigid envelopes (bifemlh, tfl, rect_fem, semimem rigid,
   glut_max3 compliant, 4.3–23.5 mm):** the deeper layers' rigid fittings
   occupy the corridor these groups want. Refinement did not clear them
   within budget; more budget or patch adjustments are the levers.
4. **transmission feasibility (sar 100%, semimem compliant 57%,
   per_brev compliant 56%, tib_ant compliant 22%, glut_max3 rigid 11%):**
   no taut configuration exists at these sizing configs; the sizing
   search already rejected alternatives as worse — remaining levers are
   config space extension or seed handling for the 10 mm groups.

## Per-muscle named-counterpart comparison (Ben's rerun request)

`per_muscle_compare.csv`: robot route vs EACH named human muscle
(per-constituent signed curves from the targets_v2 pipeline).
Compliant branch, 42 constituent rows: only 8 meet >= 99% of poses.
Representative shortfalls (meet fraction / max deficit): vas_lat_r 11% /
45.5 N*m; soleus_r 22% / 131.0 N*m; ercspn_r 0% / 95.2 N*m; vas_int_r
44% / 34.1 N*m; vas_med_r 44% / 28.8 N*m; iliacus_r 33% / 34.9 N*m;
bifemsh_r 33% / 18.4 N*m. Group-sum meets hide these: the rerun instinct
is data-justified. Per-muscle meet as an OPTIMIZATION TARGET is not wired
into the objective — awaiting your ruling.

## Dense 25-cell verification

54 evaluations verified, 0 missing. 9 jump warnings (largest: psoas
compliant, 782 N*m cell-to-cell at cell 15, coinciding with an 8%
slack-cells flip — a taut/slack transition discontinuity, implementation
class, logged with cell locations for fixing).

## Objectives, weights, budgets (of record)

- `objective_spec.json`: status TUNED-20261104-AUTO; terms
  torque_capacity_deficit w=100 (norm: sum|tau_target|),
  usable_rom_shortfall w=40 (fraction), torque_shape w=0.5
  (mean/mean), bpa_bpa_overlap_soft w=2.0 (rest-diameter depths),
  retained_bpa_corridor_soft w=2.0 (same norm). Capacity share of
  weighted terms 54–69% across representative groups (PASS >= 50%),
  shape <= 1.3% (PASS). Weight overrides for checks flow through
  opts.weightOverride; the spec file is never rewritten by checks.
- Solver: patternsearch, 600 evals/branch (screen 150), MeshTolerance
  and StepTolerance 1e-3, penaltyScale 1e4 (raised from 1e3 after gate
  r10 proved under-enforcement), refinement round 1 at 150 evals with
  non-degradation acceptance and revert.
- Constraints carry explicit units/normalization/tolerances in
  rep.cinfo; feasibility = every hard residual within its stated
  tolerance; top violator named with dimensional residual.

## Files (clickable)

- Summary tables: [deliverable_status.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/deliverable_status.csv),
  [deliverable_displacement.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/deliverable_displacement.csv)
- Dense verification: [dense_verification.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/dense_verification.log),
  [dense_verification.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/dense_verification.csv)
- Per-muscle: [per_muscle_compare.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/per_muscle_compare.csv),
  [per_muscle_compare.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/per_muscle_compare.log)
- Figures (27 + alt-text): [figures/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/figures)
- Torque curve CSVs: [torque_curves/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/torque_curves)
- Sizing: [sizing_winners.json](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/sizing_winners.json),
  [sizing_search.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/sizing_search.log)
- Weights: [objective_spec.json](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/objective_spec.json),
  [weight_sensitivity.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/weight_sensitivity.log)
- Gate: [gate_r14.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/gate_r14.log)
- run_1 state: [run_1/STATUS.json](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/run_1/STATUS.json),
  [run_1/checkpoints/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/run_1/checkpoints) (81 files),
  [run_1/logs/master.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/run_1/logs/master.log)
- Xi3 checks: [tmp/xi3_check_final.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/tmp/xi3_check_final.log)
- Packaging: [packaging_v2/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/packaging_v2) (SDF, registration report, patches, cam README)
- Targets: [targets_v2/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/targets_v2)
- Stage-1 record: [STAGE1_REPORT_20261104.md](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/STAGE1_REPORT_20261104.md)
- run_0 crash: [run_0/RUN0_CRASH_RECORD.md](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/run_0/RUN0_CRASH_RECORD.md)

## Statuses used (definitions)

- **solved** = optimizer returned AND every hard constraint within its
  stated tolerance (per-constraint table in each checkpoint rep).
- **best-infeasible** = best candidate retained (usable warm start) but
  at least one hard residual exceeds its stated tolerance; the top
  violator is named with units.
- **failed** = the solve threw; no checkpoint; retried on relaunch.
- A manifest entry, gate pass, or launched process does NOT count as a
  completed requirement — statuses above are backed by checkpoints and
  re-evaluated reports only.

## Unresolved (explicit)

- Per-muscle meet as optimization target: awaiting Ben's ruling.
- CAD-level obstacle SDF is registered and quantified but not yet wired
  into the campaign constraint path (clouds carry the margins; call site
  documented in fr_obstacles.m:cadHook and packaging_v2/README.md).
- Cam-assisted vas_int evaluator exists; not wired into the equilibrium.
- Hardware ratings unresolved (catalog capacities with unknown
  convention; flag/reject only, no approval claimed).
- 10 mm calibration file not located (own braid-law endpoints carried);
  6061-T6 bracket compliance stand-ins pending Ben.
- Dense-verify 9 jump warnings (taut/slack transition class) unfixed.
- Legacy drivers Opt_run_BiPulley.m:363 / survey2:439 still carry the
  25% contraction limit; one-argument fix required before any reuse.
