# Functional routing campaign — SIGNED-STANDARD RESULTS (2026-11-04, evening)

**Machine: easteregg2.** All paths on easteregg2 under
`D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\campaign_20261103\`
(UNC from other machines: `\\easteregg2\D$\GitHub\Bipedal_Robot\...`).

This report supersedes the results sections of RUN1_FINAL_REPORT_20261104.md
(written before Ben's figure review caught opposite-sign robot curves).
The infrastructure record remains STAGE1_REPORT_20261104.md.

## The signed standard, and what it revealed

Ben's figure review caught robot curves opposite in sign from the human
target (flex_dig example: robot +2..+9 N*m against a -2 N*m signed
target). Root cause: the objective's capacity deficit scored magnitudes
(|robot| vs |target|), sign-blind — wrong-sign designs earned full
credit. Seed diagnostic: anatomical seeds agree with human signs in
211/243 cells (86.8%); failures concentrated in flex_dig, flex_hal,
rect_fem, vas_int, sar. So the optimizer exploited sign-blindness; no
pipeline convention flip.

Fix (of record, applied everywhere): signed deficit per cell
`max(0, |tau_target| - tau_primary * sign(tau_target))` — wrong-sign
torque costs MORE than zero torque. Signed meet in deliverables, dense
verification, per-muscle compare, four-way compare. Spec
`objective_spec.json` definition updated; checkcode clean.

Consequence, stated plainly: **statuses and attainment quoted before the
fix were inflated. The signed standard is the campaign's real scoreline,
and the campaign is far from meeting it broadly at current budgets.**

## Scoreline of record (signed, full budgets, true sizing winners v2)

"Feasible" = every hard constraint within its stated tolerance (in-layer
context). "Signed meet" = robot torque in the target's direction with at
least target magnitude, per pose cell.

| group | feasible (rigid/compliant) | signed meet (rigid) | signed meet (compliant) |
|---|---|---|---|
| ext_dig | yes / yes | 100% | 100% |
| tfl | no / no | 100% | 100% |
| tib_post | no / no | 100% | 100% |
| per_long | no / no | 100% | 0% |
| lat_gas | yes / yes | 100% | 44% |
| ext_hal | no / no | 100% | 0% |
| per_brev | no / no | 56% | 78% |
| per_tert | no / no | 78% | 0% |
| bifemsh | yes / yes | 78% | 11% |
| tib_ant | yes / no | 78% | 56% |
| semimem | no / no | 78% | 0% |
| rect_fem | no / no | 67% | 56% |
| soleus | yes / yes | 44% | 0% |
| med_gas | yes / yes | 44% | 22% |
| psoas | yes / yes | 33% | 44% |
| bifemlh | no / no | 44% | 22% |
| iliacus | yes / yes | 11% | 22% |
| flex_hal | no / no | 44% | 0% |
| sar | no / no | 11% | 11% |
| grac | no / yes | 0% | 44% |
| ercspn | yes / yes | 0% | 0% |
| intobl | yes / yes | 0% | 0% |
| vas_int | yes / yes | 0% | 0% |
| flex_dig | yes / yes | 0% | 0% |
| extobl | no / no | 0% | 0% |
| glut_max3 | no / no | 0% | 11% |
| add_mag3 | no / no | 0% | 0% |

(Meet fractions from `deliverable_status.csv` re-evaluation; torque is
context-independent so these are valid. Constraint status is in-layer.)

## 2019 robotbody-routing BASELINE (same signed standard, as-is geometry)

`baseline_compare.csv` — gait2392_robotbody's existing path points,
robotbody kinematics, no optimization, both branches:

- 100% signed group-sum meet: tib_post (both), per_brev (rigid),
  ext_dig (both), glut_max3 (rigid).
- Everything else 0–89%. Notables: vas_int 0% (159 N*m max deficit),
  ercspn 0% (148), add_mag3 0% (114), soleus compliant 0% (611 N*m —
  the compliant settle at as-is geometry collapses), semimem compliant
  0% (186 N*m).

**Honest headline: neither the 2019 routings nor the current optimized
designs meet the signed standard broadly.** The optimized campaign beats
the 2019 baseline on most groups but reaches full signed coverage on
only a handful. Closing the gap is a solver-budget/strategy question,
not a bookkeeping one: the constraint-feasible "solved" status of 24
pairs does NOT imply target attainment, and under the signed objective
the optimizer is stuck near seeds for the large-torque groups (vas_int
181 N*m, ercspn 116, soleus 148, glut_max3 113 sums).

## Limiting constraints on feasible-but-not-meeting groups

- **Wrong-sign seeds (flex_dig, flex_hal):** seeds produce opposite-sign
  torque; the optimizer must re-route around bone to flip. flex_dig
  still lands at 0% signed meet despite constraint-feasible status.
- **foot_occupancy_ban (7 lower-leg groups, 3.8–28.3 mm):** several meet
  their torque (per_long rigid 100%, ext_hal rigid 100%, tib_post 100%)
  but violate the fitting-tip rule.
- **retained rigid envelopes (tfl 4.3–8.6 mm, rect_fem 23.5 mm,
  glut_max3 compliant 20.7 mm):** tfl and tib_post MEET torque 100% but
  are blocked by deeper layers' fittings.
- **centerline z=0 (add_mag3 39.8 mm, extobl 0.27 mm):** documented seed
  crossing.
- **transmission (glut_max3 rigid 11% cells, semimem compliant 57%
  slack, tib_ant compliant 22%, per_brev compliant 56% slack,
  sar/grac 100% infeasible):** no taut configuration at these configs.

## Dense 25-cell verification

54 evaluations, 0 missing, 9 jump warnings (no declared conditional
nearby) — largest psoas compliant 782 N*m at cell 15 coinciding with an
8% slack-cells flip. Taut/slack transition discontinuity class; unfixed.

## Per-muscle named-counterpart comparison (signed)

42 compliant constituent rows: 8 meet >= 99%. Group-sum meets hide
per-muscle shortfalls (vas_lat_r 11%, ercspn_r 0% / 95 N*m, soleus_r
22% / 131 N*m). Per-muscle meet as an optimization target remains
Ben-ruling-gated, not wired into the objective.

## Files

- [RUN1_FINAL_REPORT_20261104.md](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/RUN1_FINAL_REPORT_20261104.md) (infrastructure-stage report; results superseded by this file)
- [deliverable_status.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/deliverable_status.csv)
- [baseline_compare.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/baseline_compare.csv)
- [four_way_compare.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/four_way_compare.csv)
- [per_muscle_compare.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/per_muscle_compare.csv)
- [dense_verification.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/dense_verification.log)
- [deliverable_displacement.csv](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/deliverable_displacement.csv)
- [figures_signed/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/figures_signed) (27 signed figures + alt-text)
- [torque_curves/](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/torque_curves)
- [sizing_winners.json](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/sizing_winners.json)
- [objective_spec.json](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/objective_spec.json)
- [gate_r14.log](D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027/campaign_20261103/output/gate_r14.log)

## Decision points for Ben

1. Solver strategy for the deficit: current 600-eval patternsearch from
   anatomical seeds stalls near seeds under the signed objective.
   Options: multi-start/higher budgets, warm starts from high-tension
   designs, cam wiring for vas_int. No rerun started without a ruling.
2. foot-rule vs lower-leg groups: keep the hard rule (7 groups stay
   infeasible) or change the sleeve/standoff design geometry.
3. add_mag3 centerline seed: explicit on-centerline resolution needed.
4. Per-muscle meet as optimization target: yes/no.
5. The 9 dense-verify jumps: fix before trusting any psoas-class curves.
