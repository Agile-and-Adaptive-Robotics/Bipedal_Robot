# REPORT — Geometry and configurations behind dissertation Table 4.12

Generated 2026-10-31 on easteregg2. Read-only provenance extraction; no
optimization, refit, simulation, manuscript edit, or commit. The four
aggregate percentages of Table 4.12 (37.28 / 65.19 % torque coverage,
92.44 / 98.81 % active-force coverage) are quoted from
`verification.json` and were not recalculated.

## What changed vs what was inherited (27 sites, baseline vs selected)

- **Attachment/via routing: inherited at 27/27.** `spec.points`, `spec.frames`,
  `spec.cross`, grading ranges, hinge pivots, and hinge axes are bit-identical
  between the baseline and selected MAT of every site
  (`route_configurations.csv`, column `points_identical_to_baseline`).
- **Transmission configuration: changed at 27/27** (`pulley_identical_to_baseline`
  false everywhere). Selected combos (nBPA x G): nB3_G1 x10, nB3_G2 x6,
  nB2_G1 x4, nB1_G1 x3, nB2_G2 x2, nB3_G3 x2. Tackle line parts are applied on
  the primary crossing only (`tackleLineParts` per crossing: 17 sites 1/1,
  8 sites 2/1, 2 sites 3/1; the secondary crossing is always 1). Routing mode is
  `moving_via` at all 27 sites; the bowden branch was never exercised in these
  saved designs.
- **Diameters: 19 sites unchanged; 8 sites are diameter substitutions.** At
  ercspn_r, semimem_r, vas_int_r, soleus_r, med_gas_r, tib_post_r, glut_max_r,
  add_mag_r the baseline saved design is four parallel 20 mm actuators
  (`survey3_orig_*_nB4_G1` mats); the selected side retains the inherited
  40 mm specification (`diameter_identical_to_baseline` false exactly there).
- **Sizes re-solved at 25/27 sites** (`rest_delta_m`/`tendon_delta_m` nonzero);
  two sites carry identical saved sizes (see CSV).
- **Physical parallel-bundle offsets are recorded** for nBPA > 1 sites
  (`spec.bundleOffsets`, nBPA x 3); they are listed numerically in the
  configuration CSV. Physical sheave counts are NOT specified by the abstract
  transmission record (G is the tackle line-part/travel multiplier, not a
  pulley-wheel count); no packaging layout is implied.

## Saved approximations (disclosed, not corrected)

- rect_fem_r and vas_int_r: the master's builder substitutes fixed tibia point
  `t3` for the OpenSim patella function points (rigid-patella approximation);
  verbatim in each MAT's `spec.notes`. The figures show the saved survey path.
- All route points are fixed masters waypoints; the OpenSim conditional/moving
  via-point functions are not part of the saved robot survey geometry.
- Four hip sites (glut_max_r, add_mag_r, iliacus_r, psoas_r) were graded on the
  2026-10-23 coupled hip adduction/abduction envelope: 10 of 25 grid cells
  excluded per site (`excluded_cells` in the configuration CSV).

## Figures and their composition (disclosed)

- `Figures/wholeleg_saved_routes.pdf` (A sagittal x-y, B frontal z-y):
  26 routes composed into the pelvis frame. Each route's points were folded to
  its own proximal frame with its own recorded T grid at the stored pose
  nearest neutral; femur/tibia/foot frames were placed with the recorded
  rect_fem_r transforms (hip -15.00 deg, knee +10.00 deg) and soleus_r
  transforms (ankle +2.50 deg, subtalar +5.00 deg). Skeleton links connect
  recorded pivots only (rolling-knee center is inside the recorded transforms;
  the drawn knee vertex is the recorded neutral fixed-pivot fallback).
- **Route 27 (psoas_r) is not overlaid in the overview**: its chain root
  (`torso`) has no saved transform registering it to the pelvis frame. It is
  drawn in closeup panel F in its own recorded torso frame.
- `Figures/route_group_closeups.pdf`: per-group closeups with origin/via/
  insertion markers; all 27 sites appear across the figure set.
- No bone/CAD geometry is drawn (no verified registration for the masters-era
  frames); recorded pivots and points only.

## Unavailable fields

- Physical sheave count and packaging layout: not specified by the saved
  abstract transmission records.
- Per-actuator physical bundle mount geometry beyond `bundleOffsets`: not
  recorded; offsets are reported numerically, not drawn as mounts.
- torso->pelvis registration for psoas_r: not present in any saved MAT.

## Location and deliverables

Package of record (Ben, 2026-10-31: keep MATLAB work under Mesh_Optimization/Results):
`Code/Matlab/Mesh_Optimization/Results/KompareXXVII/`.
`Sources/` holds the 20261028 roster/replay/verification files fetched from eb475ws4 so the package is self-contained.

## Deliverables in this folder

- `route_configurations.csv` (54 rows: 27 sites x baseline/selected),
  `route_points.csv` (256 rows: ordered route points, body-local frames,
  roles, sources).
- `table412_configurations.tex` + `table412_configurations_review.pdf`
  (standalone compiled), label `tab:wholeleg_saved_configurations`.
- `Figures/wholeleg_saved_routes.(pdf|png)`, `Figures/route_group_closeups.(pdf|png)`,
  editable `.fig` sources + both render scripts in `Components/`.
- `provenance_manifest.json`: sha256 of every input (roster, canonical replay,
  verification, scripts, builder, masters ingredients, all 54 MATs) and every
  output; pose definitions; completeness check 27/27 sites, 54/54 MATs.
- `captions_draft.tex`: ready-to-review caption + callout (NOT inserted into
  any manuscript).

## Acceptance checks

1. All 27 Table 4.12 sites represented (26 overview + psoas closeup; 27/27 in
   both CSVs). PASS
2. Exact selected-MAT provenance per row (paths + fields in CSVs, hashes in
   manifest). PASS
3. No invented configuration or attachment update: routing inherited 27/27,
   flags recorded per site. PASS
4. Point/frame/function integrity: class rowFrame rule applied verbatim;
   frames recorded per point. PASS
5. Diameter substitution explicit (8 sites flagged, baseline 4x20 shown). PASS
6. G = tackle line parts distinguished from sheave count (stated in table
   caption, report, and manifest). PASS
7. Existing four aggregate percentages unchanged and quoted, not recomputed. PASS
8. No optimization/refit/simulation; no manuscript/bibliography/Chapter 3/
   figure-result overwrites; no commits. PASS
