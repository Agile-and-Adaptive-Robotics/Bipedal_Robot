# fr_camVasInt — vas_int cam-assisted transmission candidate

Task 3 of the physical-packaging pass (directive item 5; handoff item 3).
Function of record: `campaign_20261103/fr_camVasInt.m`. Case of record:
`fr_cases27` group 11 (`vas_int`, pbKneeMono: femur sleeve rows +
fixed-tracked moving-point tibial insertion; knee rides T1 =
`fr_chainT(kin,q,'tibia_r','femur_r')`).

## Topology

- nBPA parallel 20 mm BPAs run STRAIGHT from the femur origin bracket (the
  case's movable-1 patch, seed of record) to a cam bracket. Both are
  femur-fixed; the sleeve length is pose constant.
- A cable (artificial tendon) leaves the cam bracket exit ON the cam
  profile, wraps the taut path around the cam, exits at a tangent vertex,
  and runs straight (3D) to the tibial terminal attachment.
- Terminal attachment: 15 mm surface patch centered at the seed
  patellar-tendon insertion (the vas_med_r MovingPathPoint's neutral
  evaluated value in `tibia_r` — geometry_20261103 attention item 1
  confirms this is the model's quadriceps mechanism insertion, NOT the
  renderer's too-anterior point). Guard: the attachment may never sit
  proximal of the seed + 2 mm in the tibia frame (closure check 4), so the
  candidate is "at or inferior to the patellar-tendon attachment" per the
  handoff.
- Cam profile (POLAR, on a femoral pivot): `r(th) = r0 + a1*cos(th) +
  b1*sin(th) + a2*cos(2th) + b2*sin(2th)`, sampled at 2 deg into the
  contact polygon the cable wraps. Cam plane: the femur sagittal plane
  `z0` = the seed insertion's femur-frame z at neutral (constant — the
  OpenSim knee flexion axis is the femur-frame z, so the insertion's z does
  not drift; check 5 quantifies it, bound 5 mm).
- Cable contact: exact taut-string-around-polygon computation. From the
  bracket exit vertex the chain may run CCW or CW around the boundary to
  any tangent vertex, plus the direct no-wrap chord; the minimum feasible
  total is taken. The final straight segment must clear the polygon
  (sampled legality test, tolerance ~ half the 2 deg polygon spacing).
  Nonconvex profiles are handled (the taut path lifts off concavity).

## Design vector (14 x 1) — variables, bounds, units

| idx | variable | unit | bounds | meaning |
|---|---|---|---|---|
| 1–2 | pivot px, py | m | px in [-0.03, 0.12]; py in [-0.50, -0.33] (femur frame; knee ≈ y −0.40) | femoral pivot, cam plane z0 |
| 3 | r0 base radius | m | [0.008, 0.030] | polar profile base radius |
| 4–7 | a1, b1, a2, b2 | m | ±0.008 / ±0.008 / ±0.006 / ±0.006 | profile harmonics; r evaluated exactly per pose |
| 8 | thB bracket exit angle | rad | [0, 2π] | where the cable/sleeve bracket sits on the profile |
| 9 | sB bracket standoff | mm | [1, 12] = C.standoff_mm | radial mount standoff off the cam surface |
| 10–11 | u, v | – | [−1, 1] | tibial attachment patch coords (fr_patch convention) |
| 12 | s_tib standoff | mm | [1, 12] | attachment mount standoff |
| 13 | L0 | m | [0.4, 2.5]×cs.L0Seed | BPA rest length (bookkeeping, as the family) |
| 14 | tendon | m | [2·fit+2 mm, 2·fit+50 mm] | cable bookkeeping length |

Bounds via `[lb, ub] = fr_camVasInt('bounds', C, cs)`. The pivot box is the
distal-anterior femoral bracket zone — a mechanical-mount window, NOT a
muscle-attachment patch; the terminal attachment patch IS anatomically
tight (15 mm) per the registration RULE (clearance use never justifies
large excursions).

## Closure checks (c ≤ 0 feasible, fr_objective style)

1 `cam_min_radius` — r_min ≥ 6 mm (conservative bound on the sampled
  minimum; the exact sampled minimum feeds the taut path anyway).
2 `cam_max_radius` — r_max ≤ 45 mm (off the inflated-sleeve/skin envelope;
  the 20 mm sleeve reaches 38 mm at full contraction).
3 `cam_wrap_window` — swept contact arc ≤ 240 deg (no full wraps).
4 `attachment_proximal_guard` — attachment y ≤ seed y + 2 mm (tibia frame,
  y superior; inferior-only).
5 `cam_plane_drift` — |z_attach − z0| ≤ 5 mm over the sweep.
6 `travel_band_G1` — span range ≤ kmaxF·L0 (G = 1; tackle variants are the
  later driver's discrete choice; kept distinct from nBPA per master item 2).
7 `cable_path_clearance` — path samples ≥ 3 mm clear of the registered CAD
  obstacles (fr_obstacles_v2, clearance-only use).
8 `rest_vs_travel` — L0 ≥ spanRange/kmaxF + 30 mm.

## Plug contract (into fr_cases / fr_evalRoute / fr_objective)

`out = fr_camVasInt(C, cs, xCam, kin, qList [, opts])` consumes the case
structure and pass-1 pose list and returns:
- `out.spanCam(i)` — knee-crossing cable span per pose (the crossing-span
  analog of fr_mechBuild's exit→insertion chord);
- `out.dSpan` — dSpan/dθ on the case grid, central differences (the
  fr_objective `fr_dDeltaL` convention, fr_objective.m:297-306);
- `out.tau_perT = -dSpan` — per-unit-tension knee torque under the
  fr_objective line-86 convention `tau = -T*dSpan/dtheta`; a driver wires
  `T = mech.Ftendon` (or the cam's own transmission tension) to this.
- `out.wrapAngle(i)` — swept arc per pose, feed to `fr_hardwareLoads`
  (R = 2·T·sin(β/2)) for the cam-pivot sheave-style reaction check;
  `out.pivotReaction_N_perT` = |u_in + u_out| per unit tension.
- `out.c / out.clabels` — closure vector (nonlcon-ready).
- `out.geometry(i)` — pivot, polygon, bracket, path, exit, attachment per
  pose (visualization-ready export elements, master item 9).
- `out.mlTotal` — sleeve + cable path length (class bookkeeping analog).
- `out.femurViaNote` — the cam topology substitutes for the femur via rows
  (rows 2..cross−1 are not tension-path points here); the substitution is
  explicit, not silent.

`opts.useObstacles` (default true) enables check 7; `opts.uvOrigin`
(default [0 0 5]) sets the origin patch seed point.

## Lineage notes (Opt_run_ExtCam conventions carried over)

- Flat design vector + lb/ub + closure vector, like xCam =
  [p6x p6y p7x p7y p8x p8y tendon] in `objective_KneeExtCam20mm.m` /
  `nonlconExtCam20mm.m`; feasibility failure surfaces as a large c, and a
  driver can wrap this in patternsearch/surrogateopt exactly like
  `Opt_run_ExtCam.m` (surrogate first, patternsearch polish).
- The ExtCam clip-on points p6/p7/p8 are REPLACED here by the polar
  profile + pivot + bracket — this is a different topology (a true cam on
  a femoral pivot), per handoff item 3's "Cam profile (polar r(theta)
  parameterization), pivot location on the femur".
- The family hard constraints that transfer unchanged: no-stretch
  bookkeeping (rest/tendon bounds), standoff bounds (C.standoff_mm),
  condyle-side clearance (here: obstacle SDF check 7 instead of the knee
  context's polygon). The VTUG floor constraint (nonlconExtCam20mm) is a
  pressure-feasibility item for the later driver, not this geometry file.

## Honest status

- Implemented and verified per the check log (packaging_v2/logs/
  check_camVasInt.log): per-pose geometry, arm curve, closure, per-variable
  finite-difference sensitivity, bounds extent.
- NOT done here: wiring into the BiPam_pulley2 equilibrium (the class call
  belongs to the later driver once Ben approves the topology), the VTUG
  pressure floor, and the optimized split/harmonic search (this file is the
  evaluator + bounds of record).
