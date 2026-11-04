# targets_v2 — SIGNED per-pose human torque targets (directive item 2)

Built 2026-11-04 on easteregg2. Supersedes the constant-Fmax target curves of
`demo/fr_targets.m` (Σ Fmax·(−dL/dθ), no force-length modulation) with fully
POSE-DEPENDENT signed targets from the STOCK human model's own mechanics.

- Human reference of record:
  `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim`
  (stock human ONLY — never robot routes; opensim 4.4.1,
  `D:/Anaconda/envs/myo/python.exe`).
- Grid of record: `tmp/functional_routing_20261103/grid_export.json`,
  exported by `../fr_export_grid.m`, which calls `fr_cases27(fr_setup(),
  fr_kinematics(...))` itself and dumps each case's pass-1 qStructs
  VERBATIM — 27 groups × 9 cells, bit-identical to the consumer grids by
  construction (the targets are indexed by cell, in fr_cases cell order).

## Files

| file | content |
|---|---|
| `<group>_signed_targets.json` (27) | per-cell signed `tau_target_Nm` + full per-cell constituent breakdown (path length, rigid-tendon nce, fl_active, F_m, moment arm, per-muscle torque contribution) + convention/force-length documentation strings + peak fields |
| `signed_targets.mat` | one struct array `rows` (243 rows = group × cell), fields per cell incl. per-constituent vectors; `meta` carries provenance, convention text, and the verification result. NEW file; no stage-1 file touched; `-v7` |
| `../../fr_loadSignedTargets.m` | MATLAB loader (see below) |
| `../../fr_export_grid.m` | grid exporter (writes the tmp grid_export.json) |

Generator scripts (scratch, reusable):
`D:/GitHub/Bipedal_Robot/tmp/functional_routing_20261103/targets_v2/` —
`build_signed_targets.py`, `verify_signed_targets.py`, `test_loader.m`,
`build_signed_targets_mat.m`, `probe_targets2.py`, plus the
`sign_check_20261104.txt` / `verify_signed_targets_result.json` /
`targets_summary_20261104.json` records.

## Torque definition and sign convention (EXACT)

```
tau_target(cell) = SUM over constituent stock muscles m of  F_m(q) * r_m,dof(q)
r_m,dof(q)       = OpenSim Muscle::computeMomentArm(state, primaryDof)
F_m(q)           = Fmax_m * fl_active(nce(q))
```

`computeMomentArm` returns **−dL_m/dq**: verified this ask against a central
finite difference of `Muscle::getLength` (h = 1e-5 rad) for psoas_r,
semimem_r, tib_ant_r, vas_lat_r, add_long_r — ratio `computeMomentArm / FD
dL/dq` = −1.0000 for all five. It is therefore ALREADY the torque-convention
moment arm, `tau = F·r` equals the consumer convention
`fr_targets.m: tau = −Fmax·dL/dq`, and the meet-or-beat comparison in
`fr_objective.m` is sign-consistent. **Positive tau drives the primary
coordinate in its positive direction per the model** (see sign map).

## Force-length treatment (EXACT)

- `Fmax_m` = the model's own `max_isometric_force` (all 92 muscles are
  Thelen2003Muscle).
- `fl_active` = the model's own Thelen2003Muscle **ActiveForceLengthCurve**
  evaluated at the CE length derived by the **rigid-tendon approximation**
  `nce = (pathLength(q) − tendonSlackLength_m) / optimalFiberLength_m`,
  where `pathLength(q)` is `Muscle::getLength(state)` at the exact cell pose.
- NO muscle overrides `active_force_length_curve` in
  gait2392_simbody.osim (grep count 0), so the class-default curve object
  IS the model's own active force-length behavior. The curve clamps flat to
  `min_value = 0.1` outside its native domain (verified: fl(0.0)=0.1,
  fl(0.35)=0.1, fl(1.0)=1.0, fl(1.4)=0.5437, fl(1.8)=0.10001, fl(2.5)=0.1);
  the curve function reports its domain as ~[0.4441, 1.8123].
- Activation = 1 (maximum isometric), zero fiber velocity. Passive
  force-length is EXCLUDED. **Pennation is NOT applied** — this keeps the
  target on the same path-moment-arm convention the fr_objective/fr_targets
  consumer uses (`Fmax·|dL/dθ|` has no pennation factor either); each
  muscle's `pennation_angle_at_optimal` is recorded in the JSONs for
  auditability.

## Pose protocol (identical to geometry_20261103 and fr_fillQ27)

- Every model coordinate at its default except the group's single primary
  coordinate, set per cell. No coordinate couplers exist in the model
  (verified geometry_20261103); all joints are CustomJoints.
- Clamping: the model's `clamped` flags are all false, so OpenSim itself
  clamps nothing — the geometry-README protocol clamp is applied by the
  evaluator: requested radians clamped into the model's stored Coordinate
  range. The ONLY cell that clamps is the knee −120° request, which lands on
  the model's stored range endpoint (−2.0943951 rad = −119.99999986°),
  recorded per cell as `requested_*` vs `realized_*` with a `clamped` flag.
  Both the loader comparison (against `requested_rad`, exact fr_cases
  values) and the OpenSim evaluation (at `realized_rad`) are documented.
- Grids: 9 evenly spaced pass-1 cells from actuator_map_27 `rom_table_deg`
  (lumbar_extension −20..15°, hip_flexion_r −15..110°, hip_adduction_r
  −45..20°, ankle_angle_r −50..20°, subtalar_angle_r −25..35°, mtp_angle_r
  −40..40°), except knee primary groups, which use the campaign knee range
  **0 to −120°, no hyperextension** (fr_cases27's `C.kneeRangeDeg` override
  of the map's −120..+10). Coverage passes (pass 2, hip/ankle) exist in the
  cases but are NOT part of the target (only pass 1 is scored).

## Group mappings (verbatim from target_contract.json / handoff items 3–7)

- `vas_int` = vas_int_r + vas_med_r + vas_lat_r, elementwise signed sum at
  each same pose.
- `semimem` = semimem_r + semiten_r, elementwise signed sum at each same
  configuration.
- `add_mag3` (legacy alias preserved as the group name of record) =
  add_mag1_r + add_mag2_r + add_mag3_r + add_long_r + add_brev_r — the
  COMPLETE constituent list. Adductor census of the stock model re-verified
  this ask: the adductor-named muscles are exactly
  {add_mag1, add_mag2, add_mag3, add_long, add_brev} × {r, l}; the only
  adductor-named muscles outside the group are the LEFT-side mirrors
  (add_mag1_l, add_mag2_l, add_mag3_l, add_long_l, add_brev_l), which are
  not pulled in (campaign optimizes the right side; left is its mirror).
  No other muscle is literally named as an adductor.
- `glut_max3` = ALL nine compartments glut_max1..3_r + glut_med1..3_r +
  glut_min1..3_r as ONE aggregate baseline. **The 3-actuator compartment
  split (one actuator per max/med/min compartment) is exploratory and NOT
  approved** (handoff item 6 / contract glutealConsequence: adopting it
  would replace one group by three, changing 27 groups to 29, and requires
  a new explicit approval plus downstream re-derivation).
- `psoas` and `iliacus` = each group carries its OWN muscle's torque curve
  (psoas from the psoas_r path, iliacus from the iliacus_r path). They are
  NOT summed and iliacus force is not double-counted. The README-of-record
  statement: the ROBOT routes will share the psoas→iliacus junction
  (shared attachment/load-transfer topology per handoff item 7); the two
  target curves represent each actuator's own force contribution at that
  junction, and the orchestrator integrates the shared-junction topology on
  the robot side.

## Sign map (mandatory verification, task step 4)

Coordinate positive directions, established from the model's own conventions
and cross-checked with the produced torque signs (neutral = all coordinates
at defaults):

| dof | positive direction per model | evidence at neutral (tau, N·m) | strongly-flexed check |
|---|---|---|---|
| hip_flexion_r | flexion | psoas_r +42.88, iliacus_r +41.72 (hip flexors drive flexion positive) | at +110°: −3.38 / −2.78 (active insufficiency past the conditional range) |
| knee_angle_r | **extension** (flexion NEGATIVE, repo convention) | vas_int_r +13.45 (knee extensor, positive); semimem_r −30.72 (knee flexor, negative) | at −120°: vas_int_r +15.54, semimem_r −1.41 |
| hip_adduction_r | adduction | add_long_r +35.41 (adductor drives adduction positive) | at +20°: +18.70 |
| ankle_angle_r | dorsiflexion | tib_ant_r +33.89 (dorsiflexor, positive); med_gas_r −72.11 (plantarflexor, negative) | tib_ant_r at +20°: +27.30; med_gas_r at −50°: −45.77; med_gas_r at knee −120° (ankle 0): −1.40 (knee flexion slackens the gastroc — biarticular behavior carried by the model's own path) |
| subtalar_angle_r | inversion | tib_post_r +15.87..+29.89 (invertor, positive); per_brev/per_long/per_tert all negative (everters) | consistent at every cell |
| mtp_angle_r | toe extension (dorsi) | flex_dig/flex_hal negative (toe plantarflexors); ext_dig/ext_hal positive (extensors) | consistent at every cell |
| lumbar_extension | extension | ercspn_r +77.53..+115.85 (extensor, positive); intobl/extobl negative (flexion side) | monotone within range |

Full record: `tmp/functional_routing_20261103/targets_v2/sign_check_20261104.txt`.

Magnitude sanity: for every cell the build asserts
`|tau_target| <= SUM_m |Fmax_m · r_m|` (triangle bound with the SAME per-cell
arms), and peaks are reported per group (e.g. vas_int peak 180.82 N·m at
knee −45° vs bound 195.5 N·m; add_mag3 peak 130.48 N·m vs 137.47 N·m).

## Independent verification (task step 6)

`verify_signed_targets.py` recomputes 5 sampled (group, cell) values —
(vas_int,6), (semimem,9), (psoas,3), (tib_post,5), (med_gas,6) — through a
different API route, rebuilding each pose independently from the model's own
defaults:

- **PRIMARY route (the gate):** moment arm via CENTRAL FINITE DIFFERENCE of
  `Muscle::getLength(state)` (h = 1e-5 rad, r = −dL/dq) instead of
  `computeMomentArm`; fl from the same model curve object at the same nce.
  **Max abs difference = 6.645e-09 N·m — PASS (< 1e-6, pure numerical
  noise).**
- SECONDARY route (reported + explained): fl recomputed through
  `ActiveForceLengthCurve.printMuscleCurveToCSVFile` sampling + 5-point
  quartic interpolation; max abs difference 7.116e-05 N·m. EXPLAINED: the
  CSV printer disagrees with `calcValue` at the CSV nodes themselves (max
  6.868e-06 in fl; it evaluates the raw SmoothSegmentedFunction without the
  curve's domain remapping), which at Fmax·r ≈ 50 N·m per muscle reproduces
  exactly the observed residual. It is an approximation WITH systematic
  bias, not a faithful independent curve evaluation; the primary route
  carries the gate. Record:
  `verify_signed_targets_result.json`.

## Loader (`fr_loadSignedTargets.m`)

```
[tau, T] = fr_loadSignedTargets(groupKey, grid)  % grid = fr_cases27 case
                                                 % struct, or a numeric
                                                 % radian vector
```

Returns the nCells×1 signed `tau_target` (drop-in for `cs.target.tau`)
AFTER validating: JSON exists, swept coordinate matches `T.primaryDof`,
cell count matches, and every cell matches the fr_cases27 pass-1 `requested_rad`
within 1e-9 rad — errors loudly (`fr_loadSignedTargets:Missing/Dof/Count/
GridMismatch`) otherwise, so a misaligned grid can never silently score
against a wrong curve. Validated 2026-11-04: 27/27 groups load and validate
against the REAL `fr_cases27` case structs; perturbed-grid, missing-group,
and wrong-dof negative tests all error as designed (`test_loader.m`).

## Regeneration

```
"D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027'); addpath('campaign_20261103'); fr_export_grid('D:/GitHub/Bipedal_Robot/tmp/functional_routing_20261103/grid_export.json')"
D:/Anaconda/envs/myo/python.exe -I tmp/functional_routing_20261103/targets_v2/build_signed_targets.py
D:/Anaconda/envs/myo/python.exe -I tmp/functional_routing_20261103/targets_v2/verify_signed_targets.py
matlab -batch: test_loader; build_signed_targets_mat
```

## Not done / handoff notes

- The psoas→iliacus SHARED JUNCTION itself lives on the robot-route side;
  this pass produces the two separate own-muscle target curves (per the
  contract's targetDef) — the orchestrator integrates the topology.
- The exploratory gluteal 3-actuator split is NOT built (not approved).
- Integration of the loader into `fr_objective`/case builds is the
  orchestrator's follow-up (this task deliberately did not modify any
  `fr_*.m` beyond the two new files).
- Pass-2 coverage poses (hip/ankle) have no signed targets yet — only the
  scored pass-1 grids carry targets, matching fr_objective's scoring.
