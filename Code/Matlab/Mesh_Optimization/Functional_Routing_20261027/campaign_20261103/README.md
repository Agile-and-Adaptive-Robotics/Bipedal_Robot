# campaign_20261103 — stage-1 machinery for the Functional Routing campaign

Built 2026-11-04 on easteregg2. Implements the stage-1 campaign machinery
(MASTER_DIRECTIVE_20261103.md items 3/5/6/7/8/9 as corrected by
HARDWARE_AND_ZCODE_HANDOFF.md items 8/9/10) and the deciding gate
`sensitivity_gate(round)`.

## Runnable gate (the deciding check)

```
"D:/Program Files/MATLAB/R2025a/bin/matlab.exe" -batch "cd('D:/GitHub/Bipedal_Robot/Code/Matlab/Mesh_Optimization/Functional_Routing_20261027'); sensitivity_gate(1)"
```

`Functional_Routing_20261027/sensitivity_gate.m` is a thin shim (addpath +
call); the full gate lives at `campaign_20261103/fr_sensitivity_gate.m`.
The gate writes `campaign_20261103/output/gate_r<round>.log`, calls
`error()` (nonzero exit) on ANY failed check, and runs ~2 min (budget
ceiling 25 min). Dev-run of record: the log above, all checks (a)–(e) PASS.

## Layout

| path | role |
|---|---|
| `fr_setup.m` | paths, physics constants (620 kPa, contraction 16.5/24 %, diameter endpoints 20→38 / 13.6→17.837 mm, adopted Xi0–Xi3 from the contract), hardware families S/M/L, solver budgets |
| `objective_spec.json` | PLACEHOLDER-PENDING-BEN objective weights + hard-constraint declarations; re-read by `fr_objective` at EVERY solve (never hardcoded) |
| `kin/` | `fr_kinematics` (parses the .osim: per-joint scoped offset frames, SpatialTransform axes incl. SimmSpline knee translations), `fr_bodyT`/`fr_chainT` (validated composition: translate-then-rotate about offset-frame origins), `fr_splinemm` (natural cubic spline = SimmSpline), `fr_validateKin` (exact-transform check vs `geometry_20261103` bodyTransformsInGround, 6 poses, max 0.037 mm) |
| `routes/` | `fr_muscleroute` (declared path points incl. ConditionalPathPoint ranges and MovingPathPoint SimmSpline functions, straight from the XML), `fr_evalRoute` (implements conditional/moving behavior — nothing frozen at default pose), `fr_patch` (bounded surface patch + tangent-slide projection + mount standoff in mm), `fr_muscleProp` (max_isometric_force) |
| `mech/` | `fr_mechBuild` (maps candidate design + discrete tx config onto the mech4-verified `BiPam_pulley2_20261022`; rigid branch = Xi1/Xi2 = Inf per the class's own convention; compliant branch = adopted Xi0–Xi3 + mount/Bowden cantilever stand-ins + Xi3 bend-loss fixed point), `fr_forceLaw` (maxBPAforce + extrapolation flag), `fr_diameter` (calibrated endpoints, 10 mm own law), `fr_spans` (resolves the contract's tendon-span placeholders: max chord over the motion range + 2 ferrules (12 mm PLACEHOLDER) + stop allowance) |
| `obj/` | `fr_objective` — spec-driven terms (torque capacity deficit, usable-RoM shortfall, small shape penalty, soft BPA–BPA overlap) + hard constraints (z=0 centerline, foot-occupancy ban, required spans, transmission feasibility, rated hardware loads, bone exclusion, retained deeper-layer envelopes). NO off-axis torque penalty anywhere; secondary torques recorded in `rep.tauSecondary` |
| `pack/` | `fr_envelope` (pose-dependent sleeve/fitting/tackle primitives, WORLD frame, physically simultaneous poses), `fr_obstacles` (bone-cloud exclusion `d ≥ r + clearance` + retained-envelope capsules/spheres; documented CAD-registration hook, NOT faked), `fr_layers` (contract layerOrder), `fr_runLayered` (deep-outward pass 1 + refinement rounds with revert rule) |
| `demo/` | `fr_cases` (7 demonstrator cases, seeds = STOCK gait2392_simbody routes parsed from the XML), `fr_solveCase` (patternsearch, budgeted), `fr_targets` (human target torque curves = Σ Fmax_m · (−dL/dθ) on the STOCK routes), `fr_locForPass`/`fr_unpackX`/`fr_packX` |
| `fr_sensitivity_gate.m` + helpers | checks (a)–(e) below |
| `output/gate_r1.log` | dev-run of record |

## Demonstrator cases (ask item 5)

`hip_flex_ext` (psoas, hip −30..+90°, G=2 tackle family L), `hip_flex_extrot`
(psoas, flexion+external rotation paired profile), `hip_rot` (adductor group,
hip rotation ±30°), `knee_bowden` (med_gas, knee 0..−120°, tackle G=2 family
L, minimum central span enforced), `foot_span` (ext_hal 10 mm, ankle and MTP
passes, minimum foot tendon span, no BPA hardware on the foot), `trunk_ercspn`
(erector spinae per the geometry report's seed decision: STOCK two-point
route), `sar_layering` (sartorius, superficial demonstrator for check (d)).

## Gate checks (ask item 6)

- (a) finite-difference per variable: attachment variables must move the
  evaluated geometry AND an output; length/split variables must move an
  output (the series split's mechanical equivalence is documented in
  `fr_mechBuild`, so its check is the recorded section-length change).
- (b) bounds extent + 150-eval patternsearch on `hip_flex_ext` must not
  restore the seed.
- (c) six representative cases × {rigid, compliant}; rigid-limit recovery
  (probe pushes Xi0→0, Xi1/Xi2→1e6×, Xi3→0, wraps→1e4 — the class's rigid
  convention is INFINITE stiffness); both bracket contributions; parallel
  load sharing `Ftendon == nBPA·F/G` (0 residual); transmission work
  `nBPA·F·travel == Ft·(ΔL+gama)` on the rigid build (2.8e-16) and the
  compliant closure with the class's own cAdd channel (3.2e-16); diameter
  endpoints (20→38 handoff pair; 10 mm own law 13.6→17.837 via bpaR);
  wrapping transitions (conditional-point activation flips at the declared
  ranges; gait2392 has NO wrap objects — the conditional census is the
  functional equivalent, geometry_20261103).
- (d) LAYERING proof: deep psoas solve → retained envelope (WORLD-frame
  primitives at the superficial case's poses) → sar guide recentered ON the
  corridor (bounded 12-mm patch, construction stated in the log) → seed
  genuinely violates (−6.8 mm) → constrained solve (envelope = HARD
  obstacle, warm-started, 1e4 penalty) deflects 61 mm to feasibility →
  refinement round: superficial improves, deep not degraded (revert rule).
- (e) per-evaluation wall seconds + extrapolated full-campaign budget
  (~2.2 h for 27 groups × 2 branches at stage-1 budgets).

## Output discipline

A passing run contains ZERO occurrences of the token FAIL (all verdicts
are PASS; informational prints never embed the token — `fr_infeasTf` says
"constraint NOT binding"), zero Warnings and zero Errors on stdout, and
exits 0. The gate calls `error()` (nonzero exit) only on a real failed
check. Verified with three consecutive fresh-process runs.

## Known limitations / blockers (also in GateReport.blockers)

1. Xi3 coupling is exact-per-pose in the MonoPam knee pipeline; the bi-route
   wrapper couples the max-envelope scalar through the class's Xi0 channel
   (per-pose profiles stored on `mech.bendMeasure`); exact per-pose coupling
   needs a per-cell Xi3 channel in the class = later campaign item.
2. Mechanics grid rows drop conditional points that are not active over the
   whole pass (med_gas keeps its condyle-wrap row FIXED with the conditional
   behavior tracked in `fr_evalRoute` and the transition quantified); true
   variable-row-count grids = later campaign item.
3. Human target torque curves use Σ Fmax·(−dL/dθ) on the stock routes
   (magnitude verified against the contract's mifN as an informational
   print); the meet-or-beat optimization itself is stage-1-seeded, NOT the
   final 27-group campaign (weights are PLACEHOLDER-PENDING-BEN).
4. Hardware: families S/M/L are catalog-dimensional candidates; the rated-
   load constraint uses catalog capacity with the handoff's unresolved-
   convention caveat; the trunk case (3×20 mm at G=2) deliberately flags the
   handoff's documented overload case instead of silently resizing.
5. The 10 mm inflated-diameter endpoint (17.837 mm) comes from the repo's
   own bpaR 10 mm branch (mech2 log) — Ben's calibration was not located
   (contract unresolved[3]).
6. Objective weights are PLACEHOLDER-PENDING-BEN (objective_spec.json).
