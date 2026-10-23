# opensim_limits_checklist_20261023

Limits of the OpenSim Gait2392-lineage human reference (stock
`gait2392_simbody.osim`, OpenSim pinned at 4.4.1) that the capability-table
human surfaces must carry as per-cell flags, and the exact flag rules the
human-surface generator must apply. Binding context: `FIX_SPEC_20261023.md`
section 2 — produce a per-cell flag layer on the human surfaces; **flagged
cells still count** in the headline percentage (conservative direction);
with-flag and without-flag numbers are **both** reported; **curves are never
corrected**. Section 1 fixes the surface definition used below: OpenSim
moment arm × Fmax on a 2-D grid over the two DOFs the route crosses, axis
ranges from `rom_table_deg`, grading grid 5×5 per actuator, viewing export
9×9, non-surface DOFs pinned at 0 deg.

## 0. Inputs of record and environment (all verified live this session)

| Item | Value | How verified (command → output) |
|---|---|---|
| OpenSim version | 4.4.1 (pinned) | `"D:/Anaconda/envs/myo/python.exe" -c "import opensim; print(opensim.__version__)"` → `4.4.1` |
| Model | `Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_simbody.osim`, sha1 `b89062c763e6d8fa3d7651c0a4fd2820055d807c` | sha1sum; byte-identical (diff + sha1) to the bundled 4.6 distribution copy `C:/Users/Ben Bolen/Documents/OpenSim/4.6/Models/Gait2392_Simbody/gait2392_simbody.osim` |
| Actuator map | `Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json`, sha1 `a97afffeeae85687810986273db6ba0151eba1a2` | sha1sum |
| 1-D per-muscle audit of record | `tmp/paper_bridge/opt_improve_20261022/path_artifact_audit_gait2392_simbody_20261022_034030.csv` (92 muscle-DOF rows, 11 flagged) | re-ran its generator this session (command in its docstring); field-by-field diff of the 7 flag-relevant columns → **0 diffs**, both files 11 flagged |
| New 2-D probe (this checklist) | `probe_osim_limits_2d_20261023.csv` + `probe_osim_limits_2d_20261023.log` (40 muscle-surface rows) | `"_probe_osim_limits_2d_20261023.py"` run → `PROBE DONE`, read-only check passed |
| Step-location + slack contours | `_probe_step_locations_20261023.py` stdout (this session) | run → `DONE` |
| Determinism / sequence check | `_probe_grac_check_20261023.py` | run → three call sequences and a repeat run bit-identical (`identical=True`) |

The 27 actuators' target muscle sets are the TARGET_RULE of
`gen_human_torque_curves_singlemuscle_20261022.py` (Ben's 2026-10-22 ruling,
lines 77-122 of that file): same-named single muscles except actuator 5 =
`semimem_r` alone, 11 = `vas_int_r` alone, 24 = `glut_max1-3_r` only, 25 =
`add_mag1-3_r` only, 15 = `tib_post_r`, 16 = `tib_ant_r`.

## 1. L1 — Conditional (via) path-point activation jumps

**What it is.** A `ConditionalPathPoint` is "a point that is active only for
a specified range of a coordinate"; its range property is "the minimum and
maximum values that the coordinate can range between, for which the
PathPoint is active" (official source: `ConditionalPathPoint.h`,
opensim-core). When the point switches on or off, the path polyline gains or
loses a vertex: the muscle length, and therefore the moment arm
(arm = −dL/dq), changes discontinuously at that boundary. OpenSim's own Moco
guidance tells users to "check if there are any large discontinuities in the
model's moment arms ... e.g. from wrapping geometry or conditional path
points engaging/disengaging" (forum, source F3), and the opensim-core team
shipped topology-based moment-arm detection for `PolynomialPathFitter` in
2026 precisely because this artifact class is real (issue #4348, official
tracker).

**Live behavior on our model (verified this session).** The boundary is a
straight line in the 2-D graded box: `act_boundary_span_deg = 0.0` on every
probed surface (`probe_osim_limits_2d_20261023.csv`), because
`isActive` depends only on the point's own coordinate. This is exactly what
the ConditionalPathPoint semantics imply, now confirmed numerically. Two
consequences the 1-D audit could not show:

1. A conditional point keyed on the **secondary** DOF draws its boundary
   across the surface at constant secondary angle. The gastrocnemii are the
   case: `med_gas_r` / `lat_gas_r` carry `knee_angle_r:-45.00..10.00`
   conditional points, so on the ankle × knee surfaces the via-point line is
   knee = −45 deg (all knee angles below −45 deg, i.e. flexed past 45 deg,
   activate the via point) at every ankle angle. The 1-D audit flagged these
   only on the knee-secondary sweep (`activation_angles_deg = -45.25`);
   on the 2-D surface the line spans the full ankle axis.
2. The **jump magnitude depends on the other angle**. Gracilis: the knee
   −25 deg boundary jump grows from 1.60% of column-peak arm at hip = −15 deg
   to 7.04% at hip = 80 deg (measured at 1 deg sampling; check log of
   `_probe_grac_check_20261023.py`). The 1-D audit at hip = 0 deg reported
   only 0.81%, so the 1-D number understates the jump the 2-D surface
   actually contains at flexed hip.

**Which of the 27 it touches** (boundary lines inside the graded boxes; all
constant in the other angle):

| Act | Muscle | Boundary line (XML range; audit sample midpoint) |
|---|---|---|
| 5 | semimem_r | knee = −32.00 deg (audit −32.25) |
| 8 | rect_fem_r | knee = −83.65 deg (audit −83.75) |
| 9 | grac_r | knee = −25.00 deg (audit −25.25); jump grows with hip flexion |
| 11 | vas_int_r | knee = −81.36 deg (audit −81.25) |
| 13 | med_gas_r | knee = −45.00 deg on the ankle × knee surface (audit −45.25) |
| 14 | lat_gas_r | knee = −45.00 deg (audit −45.25) |
| 26 | iliacus_r | hip flexion = 45.00 deg (audit 45.25) |
| 27 | psoas_r | hip flexion = 45.00 deg (audit 45.25) |

All other 19 actuators' target muscles have zero conditional points (audit
`n_cond_pp = 0`).

**Per-cell flag rule (FLAG_VIA).** For each surface (primary P, secondary
Sd), precompute a fine sweep at 0.5 deg (audit-consistent) along P at each
graded secondary value and along Sd at each graded primary value, recording
the path-point `isActive` pattern. Flag a cell `FLAG_VIA` if any
isActive-set change falls inside the cell's angular bin (primary bin at the
cell's secondary value, or secondary bin at the cell's primary value). The
boundary is a line of constant angle (verified), so this is equivalent to
flagging the cell row/column crossed by the line.

## 2. L2 — MovingPathPoint nonsmoothness (patella-less quadriceps, jitter regimes)

**What it is.** A `MovingPathPoint` relocates every instant as a function of
its coordinate. The Gait2392/2354 documentation states: "Seth removed the
patella to avoid kinematic constraints; insertions of the quadriceps are
handled with moving points in the tibia frame" (official model page, source
D1). Moving points make muscle length a piecewise-smooth function with
kinks; the moment-arm formula for moving points was wrong before OpenSim
3.3 and was fixed to the correct equation (forum t=18774, A. van den Bogert;
forum-labeled but from an authoritative author; fixed long before our 4.4.1,
so 4.4.1 arms are correct but nonsmooth).

**Where it touches the 27.** Only the target muscles with moving points:
`rect_fem_r` and `vas_int_r` (audit `n_moving_pp = 1` each). Their step
clusters sit near knee extension: vas_lat/vas_med (dropped from our targets
but documented in the audit) step at 3.75-7.25 deg; `vas_int_r` shows 6
events at 1 deg sampling with max 3.55% of column peak, and 1.84% at the
0.5 deg audit pass (audit `step_angles_deg` for vas_int_r is empty at 0.5
deg, `max_step_pct = 1.8372`, just under the 2% threshold). Separately, the
iliacus/psoas "jitter" regimes are **not** via-point chatter: their event
magnitudes scale with the sampling step (1 deg max 8.55% vs 0.5 deg max
4.41% for iliacus, ratio 1.94; psoas 8.21% vs 4.11%, ratio 2.00 — a true
discontinuity would not scale), so they are smooth-but-steep gradients of
the arm as the path sweeps near its via point at high hip flexion. Iliacus:
25 events at 1 deg spanning hip 45.5 deg and 86.5-98.5+ deg; psoas: 27
events spanning 83.5-96.5+ deg (fine-pass lists in
`probe_osim_limits_2d_20261023.csv` columns `fine_steps_sec0_deg`; the 0.5
deg maxima 4.4058% / 4.1096% match the audit's numbers exactly).

**Per-cell flag rule (FLAG_STEP).** Same fine 0.5 deg sweeps as L1. Let
A* = max |arm| over the whole surface (surface-wide normalization, matching
the audit's detector; do **not** normalize per column — per-column peaks
inflate ratios near arm sign crossings, seen this session in exploratory
runs). Flag a cell `FLAG_STEP` if any adjacent-sample |Δarm| > 2% of A*
inside the cell's bin on either axis. Expected to fire at 0.5 deg for
actuators 26, 27 (jitter regimes), 18, 17 (subtalar −24.5 to −11.5 deg at
plantarflexed ankle), 15 (subtalar 12.5-27.5 deg at plantarflexed ankle),
and 9 at flexed hip; borderline (may or may not cross 2% at 0.5 deg): 8, 11.
Marginal 1 deg-only regimes (6 sartorius near knee extension, 10, 19, 21,
7, 24, 25 at ~2.1-2.6%) are expected **not** to fire at 0.5 deg with
surface-wide normalization; whatever the generator's detector reports is
what gets flagged — no manual overrides.

## 3. L3 — Wrapping-solver artifacts: structurally absent in this model file

**What it is (in general).** Path wrapping (`PathWrap` over `WrapObject`
surfaces) is the mechanism GeometryPath uses to shape paths around
obstacles (official source: `GeometryPath.h`, opensim-core). Wrapping
solvers have a documented artifact history: Release 1.5 fixed "moment arms
for some muscles exhibit glitches at the min/max of a coordinate's range of
motion, due to limitations of the difference-based calculation method"
(bug 547) and "when a muscle wraps less than 2mm over a cylinder or sphere,
the wrapping algorithm ignores the wrapping" (bug 548) (official Release
History, source D6); the legacy SimTK tracker item #2014 "Unexpected Moment
Arm Values" exists on the same theme (title seen in search results; page
returned HTTP 403 this session, content **not verified**). Post-4.1, the
wrapping calculation itself changed: a user documented that "something
changed in OpenSim 4.2 involving wrapping" that did not appear in the 4.2,
4.3, or 4.4 release notes, produced "wonky" fiber lengths and moment arms
and paths "through the shoulder" in wrapped models, and that after removing
the offending wrap surfaces "I got identical force length/moment arm curves
in 4.2/4.3/4.4 as the original models had in 4.1" (forum t=16953; forum
evidence, but the mechanism is corroborated by the official 4.4 release
notes listing multiple wrapping fixes — performance, multithreading, and
visualization items only, none claiming to change moment arms; source D7).

**Why it does not apply to our surfaces.** The stock `gait2392_simbody.osim`
contains **zero `PathWrap` entries and 13 empty `WrapObjectSet`s** (grep of
the model XML this session: 0 `<PathWrap name` matches; 184 muscle tags, 52
`ConditionalPathPoint`, 16 `MovingPathPoint`), and the live API confirms it:
`live_wrap_count = 0` on all 40 probed muscle-surfaces
(`probe_osim_limits_2d_20261023.csv`). All path shaping in this model is
static + conditional + moving points. A wrapping-solver artifact therefore
cannot occur in the human surfaces as generated from this file; `wrap_flips
= 0` in the audit of record by the same construction.

**Per-cell flag rule.** None (nothing to flag). Instead the generator keeps
two **halt gates**: (a) live `getWrapSet().getSize() > 0` on any target
muscle, or (b) any wrap-flip detected, → abort generation and escalate.
The version pin (below) is the second half of this limit's protection.

## 4. L4 — Muscles too short for Fmax: tendon slack and sub-physiological fiber lengths

**What it is.** The human surface is τ_h = arm × Fmax, which assumes the
muscle can deliver full isometric force at every configuration. The
Thelen2003Muscle implementation says otherwise at short lengths: tendon
force is exactly zero at or below slack length (`fse = 0` for strain ≤ 0;
official source `Thelen2003Muscle.cpp`), the fiber-length clamp "does not
prevent the muscle fiber from becoming shorter than is physiologically
possible (that is shorter than approximately half a normalized fiber
length)" (official doxygen, source D5), and the active force-length curve
"is a Gaussian, which is always greater than 0" with shape
`exp(−(lceN−1)²/0.45)` (doxygen D5 + source D4), so deliverable active force
decays smoothly but never reaches zero — the surface overstates worst where
the tendon cannot even carry tension.

**Live measurements (this session, full graded boxes).** Normalized fiber
length was estimated as `lceN_est = sqrt((l_mt − tsl)² + (lo·sinφ0)²)/lo`
(constant-width pennation model, φ0 = pennation at optimum from the model).
This is an **upper bound** on true fiber length when the tendon carries load
(tendon strain only shortens the fiber further), so `lceN_est < 0.5`
guarantees the documented sub-physiological regime. Tendon-slack cells
(`l_mt ≤ tsl`, 1 deg × 2 deg grid):

| Act | Muscle | Slack region in the graded box | Cells |
|---|---|---|---|
| 5 | semimem_r | knee −120 to −119 deg with hip −15 to +15 deg (also slack at hip = 0 on the primary-only slice: min l_mt 0.3445 < tsl 0.359) | 537 |
| 7 | tfl_r | knee = +10 deg (full extension edge) with hip ≥ ~94 deg | 669 |
| 13 | med_gas_r | ankle = −50 deg (plantarflexion edge) with knee ≤ −60 deg | 384 |
| 14 | lat_gas_r | ankle = −50 deg with knee ≤ −78 deg | 150 |

Minimum `lceN_est` over each graded box (values below 0.5 = sub-physiological
fiber; from `probe_osim_limits_2d_20261023.csv`):

- Below 0.5: tfl 0.052 (knee −106, hip 104), lat_gas 0.139 (ankle −44, knee
  −96), rect_fem 0.121 (knee 10, hip 110), semimem 0.259, bifemlh 0.263,
  glut_max1 0.288, med_gas 0.292, add_mag1 0.339, psoas 0.363, iliacus
  0.376, glut_max2 0.448-0.448 (adduction/rotation boxes 0.448/0.409),
  soleus 0.423, vas_int 0.446.
- Borderline 0.5-0.55: intobl 0.532, flex_dig 0.510, flex_hal 0.515.
- All other target muscles stay above 0.55 (see the CSV for every value,
  with the configuration of each minimum in `min_lcen_est_at_deg`).

No target muscle's minimum dips below 0.01·lo (the implementation clamp
floor), but a dozen sit below the physiological half-length floor somewhere
in their graded boxes.

**Per-cell flag rules.** Flag a cell `FLAG_SLACK` if `l_mt(θp, θs) ≤ tsl`
for any target muscle in the actuator's set (for the three-compartment sets
24 and 25, flag if any member triggers; record which). Flag a cell
`FLAG_SHORT` if `lceN_est(θp, θs) < 0.5` for any target muscle. Both flags
are caveats on cells that still count; the reported with/without-flag
numbers in FIX_SPEC §2 make the conservative direction explicit.

## 5. L5 — Muscle paths routing through bones

**What it is.** With zero wrap surfaces in the model (L3), every muscle path
is a straight-segment polyline between path points; nothing in the file or
the solver prevents a segment from intersecting bone geometry. The official
model page documents only the patella workaround (D1) and is silent on bone
clearance for the other paths. The forum record shows what wrap degradation
does to wrapped models (paths "through the shoulder", F1), which is
structurally inapplicable here, but the class "path crosses a bone" can
neither be affirmed nor excluded for this model without a geometry check.
Related caution from the official tracker: "the moment arm solver will
occasionally return [non-zero results for] coordinates or joints that the
path does not physically cross" (opensim-core #4347) — one more reason our
surfaces take their DOF sets from Ben's actuator map (the routing of
record), never from computed arms, and one more reason the STANDING TRAP
rule of FIX_SPEC §1 (never route human muscles on robot/BPA paths) stands.

**Status of the check.** The bone-penetration test against the bone
geometry (point clouds in `Code/Matlab/Bone_Mesh_Plots/Open_Sim_Bone_
Geometry/` or the MyoConverter STLs) is **TODO and not run** — this is
stated, not papered over, exactly as the audit of record states it
(`audit_osim_path_artifacts_20261022.py` docstring, item 4).

**Per-cell flag rule.** None until the penetration check exists and is run;
no flags are invented for this limit. When it is built, the rule will be:
flag `FLAG_BONE` for cells whose configuration's path segments intersect a
bone surface; until then the capability table carries this checklist as the
disclosure.

## 6. L6 — GeometryPath behavior across OpenSim versions (pinned at 4.4.1)

**What it is.** GeometryPath calculations have changed across releases more
than once: the pre-3.3 moving-point moment-arm formula was wrong (F2);
GeometryPath/wrapping behavior changed in 4.2 without a release-note entry
(F1); the official 4.4 notes list wrapping fixes that are
performance/multithreading/visualization-scoped only (D7). Moment-arm values
are therefore version-dependent for models with wrapping, and the pin is
load-bearing even for a wrap-free model.

**Standing protections (all verified live this session).** The generator's
`PINNED_VERSION` gate (abort on `opensim.__version__ != '4.4.1'`), the sha1
pins on model and map (values in section 0), and this session's replay
check: the audit of record re-ran bit-consistently (0 field diffs across
the 7 flag-relevant columns of all 92 rows). Stock-model identity to the
distribution copy is sha1-verified. Determinism of the API path itself was
checked directly: three call sequences (arm-only; arm interleaved with
`getLength`/`isActive`; fresh-state-per-column) and a repeat run produced
bit-identical gracilis columns (`identical=True`).

**Per-cell flag rule.** None; this limit is enforced by halt gates
(version, both sha1s) already in the generator, which must stay.

## 7. Consolidated per-actuator flag expectations (gate table)

The generator must produce at least these flags; if any listed flag is
missing after its 0.5 deg surface-wide run, its detector is broken — halt,
do not publish. "—" = no flag expected from L1/L2/L4.

| Act | Muscle target | FLAG_VIA | FLAG_STEP regime | FLAG_SLACK | FLAG_SHORT |
|---|---|---|---|---|---|
| 1 | ercspn_r | — | — | — | — |
| 2 | intobl_r | — | — | — | borderline (0.53) |
| 3 | extobl_r | — | — | — | — |
| 4 | bifemlh_r | — | — | — | yes (0.26; deep knee + extended hip) |
| 5 | semimem_r | knee −32 deg | — | knee −120/−119 deg corner | yes (0.26) |
| 6 | sar_r | — | marginal only (near knee extension) | — | — |
| 7 | tfl_r | — | marginal only | knee +10 deg, hip ≥ ~94 deg | yes (0.05; extreme) |
| 8 | rect_fem_r | knee −83.65 deg | borderline at boundary | — | yes (0.12; knee 10, hip 110) |
| 9 | grac_r | knee −25 deg (jump grows with hip flexion to ~7% at hip 80) | expected at flexed hip | — | — |
| 10 | bifemsh_r | — | marginal only | — | — |
| 11 | vas_int_r | knee −81.36 deg | borderline near extension | — | borderline (0.45) |
| 12 | soleus_r | — | — | — | borderline (0.42) |
| 13 | med_gas_r | knee −45 deg line (secondary axis) | — | ankle −50 deg, knee ≤ −60 deg | yes (0.29) |
| 14 | lat_gas_r | knee −45 deg line | — | ankle −50 deg, knee ≤ −78 deg | yes (0.14) |
| 15 | tib_post_r | — | subtalar 12.5-27.5 deg at plantarflexed ankle | — | — |
| 16 | tib_ant_r | — | — | — | — |
| 17 | per_brev_r | — | subtalar −24.5 to −11.5 deg at plantarflexed ankle | — | — |
| 18 | per_long_r | — | same regime, strongest (10.2% of column peak at 1 deg) | — | — |
| 19 | per_tert_r | — | marginal only | — | — |
| 20 | flex_dig_r | — | — | — | borderline (0.51) |
| 21 | flex_hal_r | — | marginal only | — | borderline (0.51) |
| 22 | ext_dig_r | — | — | — | — |
| 23 | ext_hal_r | — | — | — | — |
| 24 | glut_max1-3_r | — | marginal only (1/2) | — | yes (0.29-0.45; hip extension + adducted corner) |
| 25 | add_mag1-3_r | — | marginal only (3) | — | yes (0.34-0.55; adducted + flexed corner) |
| 26 | iliacus_r | hip 45 deg | 86.5-98.5+ deg jitter (steep-smooth) | — | yes (0.38) |
| 27 | psoas_r | hip 45 deg | 83.5-96.5+ deg jitter (steep-smooth) | — | yes (0.36) |

"Borderline" entries are within sampling noise of the threshold; the
generator's own detector output decides them, and either outcome is
acceptable — only the "yes" and named-boundary cells are hard gates.

## 8. Standing rules the generator must keep (from FIX_SPEC sections 1-2 and §6)

1. Flagged cells still count in the headline meet percentage; report
   with-flag and without-flag counts side by side; never correct a curve.
2. Flags are reason-coded (`FLAG_VIA`, `FLAG_STEP`, `FLAG_SLACK`,
   `FLAG_SHORT`; multiple allowed per cell) and carried in CSV columns
   alongside the existing audit caveat strings.
3. Halt gates: version ≠ 4.4.1; model or map sha1 mismatch; any live wrap
   set non-empty; any wrap flip; gate-table flags missing.
4. Detector spec: 0.5 deg sampling, surface-wide peak normalization, both
   axes swept at each graded value of the other (per L1/L2 above).
5. STANDING TRAP: human surfaces are computed on OpenSim model paths only;
   robot/BPA paths are never substituted.
6. Knee finite-difference validation stays informational-only (existing
   generator behavior; wrap/moving-point nonsmoothness makes FD tolerances
   inapplicable as a hard gate at the knee).

## Sources

Official documentation (docs/doxygen/source code/issue tracker):

- D1 OpenSim documentation, "Gait 2392 and 2354 Models" page
  (opensimconfluence.atlassian.net, page id 53086215), fetched 2026-10-23:
  "Seth removed the patella to avoid kinematic constraints; insertions of
  the quadriceps are handled with moving points in the tibia frame."
- D2 `ConditionalPathPoint.h`, opensim-core (GitHub), fetched 2026-10-23.
- D3 `GeometryPath.h`, opensim-core (GitHub), fetched 2026-10-23.
- D4 `Thelen2003Muscle.cpp`, opensim-core (GitHub raw), fetched 2026-10-23:
  `fal = exp(-x/KshapeActive)`, `KshapeActive = 0.45`, `FmaxTendonStrain =
  0.04`, `fse = 0` below slack, `clampFiberLength`, `MuscleCannotEquilibrate`.
- D5 OpenSim/Moco doxygen, class Thelen2003Muscle, fetched 2026-10-23:
  Gaussian active force-length "always greater than 0"; clamp "does not
  prevent the muscle fiber from becoming shorter than is physiologically
  possible (that is shorter than approximately half a normalized fiber
  length)".
- D6 OpenSim documentation, Release History page (page id 53091041),
  fetched 2026-10-23: Release 1.5 bug fixes 547 and 548 (moment-arm
  glitches at range extremes; sub-2 mm wrap ignored). Historical context
  only; not a current known issue at 4.4.1.
- D7 OpenSim documentation, "What's New in OpenSim 4.4" (page id 53088140),
  fetched earlier this session: wrapping items are performance,
  multithreading, and visualization fixes; none claim moment-arm changes.
- D8 opensim-core issue #4347 "Add AbstractGeometryPath utility for
  detecting dependent coordinates" (official tracker), fetched 2026-10-23:
  moment-arm solver "will occasionally return" non-zero results for
  coordinates the path does not physically cross.
- D9 opensim-core issue #4348 "PolynomialPathFitter: use topology-based
  moment arm detection" (official tracker, closed 2026-08-18), fetched
  2026-10-23.

Forum (labeled as such; not official documentation):

- F1 simtk.org forum thread t=16953, "Change in Muscle Geometry Path
  Calculations post-4.1", fetched 2026-10-23.
- F2 simtk.org forum thread t=18774 (moving-path-point moment-arm formula
  fix history), fetched earlier this session.
- F3 simtk.org forum post p=34178 (Moco guidance to check moment-arm
  discontinuities from wrapping/conditional points), fetched earlier this
  session.
- F4 Legacy SimTK tracker item #2014 "Unexpected Moment Arm Values": title
  from search results only; page returned HTTP 403 this session; content
  not verified.

Local artifacts (this repository, this session):

- L1 `tmp/paper_bridge/opt_improve_20261022/path_artifact_audit_gait2392_
  simbody_20261022_034030.csv` and its generator `audit_osim_path_artifacts_
  20261022.py` (re-run; 0 diffs).
- L2 `tmp/paper_bridge/capability_table_20261022/probe_osim_limits_2d_
  20261023.csv` (+ log; generator `_probe_osim_limits_2d_20261023.py`).
- L3 `tmp/paper_bridge/capability_table_20261022/_probe_step_locations_
  20261023.py` output (event locations; its "worst" lines carry a percent-
  vs-fraction comparison bug and are superseded by L2/L4 magnitudes).
- L4 `tmp/paper_bridge/capability_table_20261022/_probe_grac_check_
  20261023.py` output (sequence-independence, determinism, hip-dependence
  of the gracilis boundary jump).
- L5 `Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json`
  (RoM table, DOF pairs) and `gen_human_torque_curves_singlemuscle_
  20261022.py` (TARGET_RULE, PINNED_VERSION gate, caveat machinery).
