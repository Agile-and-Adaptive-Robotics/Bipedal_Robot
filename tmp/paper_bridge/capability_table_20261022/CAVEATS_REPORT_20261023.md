# CAVEATS_REPORT_20261023

Plain-language caveats for the v2 27-actuator capability table
(capability_table_20261023.md / .tex / .xlsx) and the per-muscle rows
(per_muscle_capability_20261023.csv), assembled 2026-10-23 on easteregg2
under `FIX_SPEC_20261023.md`. Every number re-derives from the CSVs of
record in this folder: `human_surfaces_5x5_grading_20261023.csv` (the single
grading source for BOTH rows), `hand_tuned_metrics_20261023.csv`
(hand-tuned row), `BiPulley_survey2_perMuscle_20261023.csv` (optimized row),
cross-checked against `human_surfaces_manifest_20261023.json`. The assembler
(`assemble_deliverables_20261023.py`) asserts every count below before
writing; all asserts passed in this run.

## 1. What changed from the v1 table (2026-10-22)

- Metric: v2 grades each actuator on the AREA SHARE of its coupled-DOF 5x5
  surface (the two joint coordinates the route crosses), replacing v1's 1-D
  primary-DOF slice (FIX_SPEC sections 1/1b).
- Optimized row: fresh solves against the single-muscle surfaces (all 27
  built and graded; v1 had 23 built, 16 gradeable).
- Hip adduction/abduction coupling (FIX_SPEC 1b, Ben 2026-10-23) is applied;
  v1's quarter-ellipse placeholder is retired.
- v1 files (`*_20261022*`) are untouched.

## 2. Covered counts

- **Hand-tuned row: 27 of 27 covered, 0 excluded.** Coverage source is
  `connor-corrected` on every row (old-code routing, old-code BPA sizes,
  620 kPa).
- **Optimized row: 27 of 27 built and graded.** `built = true` and
  `converged = 1` on every row of the survey2 CSV (builder coverage 27/27).
- No muscle is excluded from either row.

## 3. Exclusions that ARE in the metric: hip-envelope cells (not muscles)

On the four surfaces whose adduction axis is graded (Gluteus Maximus,
Adductor Magnus, Iliacus, Psoas X, acts 24-27), 10 of 25 grid cells each lie
outside the coupled hip adduction/abduction envelope and are REMOVED from the
area denominator for BOTH rows (never counted as failures): 40 cells total.
Rule (FIX_SPEC 1b): flexion axis -15 deg (max extension) to +110 deg (max
flexion); adduction maximum tapers linearly from +20 deg at full extension to
+10 deg at full flexion, abduction maximum from -10 deg to -60 deg; a cell is
graded when its center is inside the closed envelope. Physical origin: BONE
INTERFERENCE (femur-pelvis contact, per Ben), a real anatomical/hardware
constraint. Graded cells per hip actuator: 15. `hip_rotation` stays the
static rom-table box (pinned at 0 on act 24); act 7 (TFL) has hip_adduction
pinned at 0, which is inside the envelope at every graded flexion, so no
envelope exclusion applies to its surface.

## 4. Artifact-flagged cells, and both percentage versions

Flag layer on the human surfaces (per-cell; curves are NEVER corrected):
FLAG_VIA (conditional path point active-set change in the cell bin),
FLAG_STEP (adjacent-sample moment-arm jump over 2% of the surface-wide peak),
FLAG_SLACK (musculotendon at or under tendon slack length), FLAG_SHORT
(estimated normalized fiber length under 0.5). FLAG_BONE (bone penetration)
is NOT APPLIED: the check is still TODO per the checklist, no flags were
invented.

Totals (re-derived from the grading CSV; equal to the manifest): 635 graded
cells, 40 envelope-excluded cells, **126 flagged graded cells**. Per-actuator
flagged counts are in `surfaces_index_20261023.md` and the per-muscle CSV
column `cells_flagged_in_domain`.

**Flagged cells count in the headline percentages (conservative). Both
versions, as row averages of the per-muscle shares:**

- Hand-tuned row: WITH flags meet **72.49%** / active **96.64%**;
  WITHOUT flags meet **73.21%** / active **95.23%**.
- Optimized row: WITH flags meet **65.19%** / active **74.27%**;
  WITHOUT flags meet **67.00%** / active **72.75%**.

## 5. Transcription suspects (hand-tuned row; carried over VERBATIM from v1)

Nine actuators carry a transcription-suspect flag because the corrected
Connor script differs from Ben's reference values at one or more routing
points (the hand CSV header states "suspect+note carried over VERBATIM from
hand_tuned_metrics_20261022.csv"; this pass verified 9 `suspect = yes` rows).
**All were graded as-is**; none was silently corrected. A suspect flag means
the hand-tuned percentage could shift once Ben rules on the value; it does
not mean the grade is wrong. Hand percentages below are the v2 (2026-10-23)
re-grade values.

- **Sartorius** (act 6, hand meet 36.0000%, active 96.0000%): VALUE MISMATCH: col2 vs Ben sar_wr1: Connor [[0.091 -0.105 0.037]] Ben [[0.091 -0.105 0.038]]; col7 vs Ben t6: Connor [[0.038 -0.382 -0.004]] Ben [[0.038 -0.382 0.004]]
- **Gracilis** (act 9, hand meet 20.0000%, active 100.0000%): origin differs from Ben p8; Connor y z equal rect_fem origin p10's (likely transcription slip); VALUE MISMATCH: col1 vs Ben p8: Connor [[-0.074 0.125 0.079]] Ben [[-0.074 -0.119 0.028]]
- **Bicep Femoris, Short Head** (act 10, hand meet 100.0000%, active 100.0000%): VALUE MISMATCH: col3 vs Ben twr3: Connor [[-0.023 -0.056 0.029]] Ben [[-0.023 -0.056 0.034]]
- **Tibialis Posterior** (act 15, hand meet 88.0000%, active 96.0000%): MIF 3549 = map master_mif_n but Ben PAM_Prop has tib_post_mif 1588 (map documents 3549 as soleus copy-paste artifact; same 40 mm diameter either way); routing values match Ben
- **Peroneus Longus** (act 18, hand meet 36.0000%, active 96.0000%): VALUE MISMATCH: col1 vs Ben t2: Connor [[0.005 -0.062 0.047]] Ben [[0.006 -0.062 0.047]]; col5 vs Ben per_long_wr4: Connor [[0.048 0.011 0.028]] Ben [[0.068 0.011 0.028]]
- **Peroneus Tertius** (act 19, hand meet 80.0000%, active 96.0000%): three components differ from Ben (values appear shuffled between t2/per_tert_wr1/a2); VALUE MISMATCH: col1 vs Ben t2: Connor [[0.006 -0.062 0.079]] Ben [[0.006 -0.062 0.047]]; col2 vs Ben per_tert_wr1: Connor [[0.023 -0.407 0.022]] Ben [[0.023 -0.407 0.016]]; col3 vs Ben a2: Connor [[0.079 0.016 0.034]] Ben [[0.079 0.022 0.034]]
- **Extensor Digitorum Longus** (act 22, hand meet 40.0000%, active 96.0000%): VALUE MISMATCH: col5 vs Ben ext_dig_wr4: Connor [[0.005 0.005 0.015]] Ben [[0 0.005 0.015]]
- **Adductor Magnus** (act 25, hand meet 93.3333%, active 93.3333%): Location columns 2 and 3 both equal add_mag_wr2; Ben routing uses add_mag_wr1 then add_mag_wr2 (wr1 absent; likely transcription slip); VALUE MISMATCH: col2 vs Ben add_mag_wr1: Connor [[-0.059 -0.108 -0.03]] Ben [[-0.136 -0.093 0.059]]
- **Psoas X** (act 27, hand meet 66.6667%, active 100.0000%): PLocation undefined-variable slip (values = Location4; map: presumed Location4); b1/wr1 frame roles swapped vs Ben's inline formula; point values match Ben

## 6. Unconverged or capped rows, and budget disclosure

- **Convergence: none failed.** `converged = 1` and `stopReason = converged`
  on all 27 optimized rows.
- **No row hit the evaluation cap.** Evaluations per muscle: min 68
  (Vastus Intermedius), median 163, max 252 (External
  Oblique), all under the 300-evaluation cap. Recorded wall time 0.12-0.37
  min per muscle, far under the roughly 10-min cap. Parpool was capped at 8
  per the FIX_SPEC; the survey CSV of record does not record the pool size
  actually used, so that is stated as the cap, not as a measured value.
- Ladder: nBPA {1,2,3} x tackle G {1,2,3} with sizes re-solved per
  configuration, resting lengths as upper limits, envelope guards; the
  recorded winner is the best-meeting design (an active tie-break never
  records a lower-meeting design; FIX_SPEC section 4).
- **Configuration ceiling:** only 4 of 27 winners reach 100.00%
  meet (Erector Spinae (act 1), Bicep Femoris, Short Head (act 10), Soleus (act 12), Tibialis Posterior (duplicate name) (act 16)). The rest fall short of full coverage within the
  allowed configuration set and budget; the ladder chains in the survey2 CSV
  are the record. This is a capability survey, not a design campaign: the
  budget answers "can this routing class reach the human surface at all",
  and the percentages are not flight-ready designs.

## 7. Active-force zeros in the optimized row

Six optimized winners record active 0.0000% while still meeting part of the
surface: Erector Spinae (act 1, opt meet 100.0000%), Semimembranosus (act 5, opt meet 32.0000%), Medial Gastroocnemius (act 13, opt meet 20.0000%), Lateral Gastrocenemius (act 14, opt meet 20.0000%), Gluteus Maximus (act 24, opt meet 53.3333%), Adductor Magnus (act 25, opt meet 60.0000%). The winners are taut where they meet but have zero
positive-strain cells on the grading grid. The values are recorded as-is
from the CSV of record; treat those rows' active cells with care until the
strain-window question is resolved (the same class of issue v1 flagged on
Bicep Femoris Long Head). The hand-tuned row has no active zeros (row minimum
80.00%).

## 8. Two-hinge scope for 3-plus-joint muscles, and constant-secondary rows

Each actuator is graded on the TWO coordinates the route crosses; remaining
joints the muscle crosses are pinned at 0 deg: act 1 lumbar_rotation; act 7 hip_adduction_r; act 20 subtalar_angle_r; act 21 subtalar_angle_r; act 22 subtalar_angle_r; act 23 subtalar_angle_r; act 24 hip_rotation_r. Moment arms
the human model develops at other joint configurations are not represented;
this is a per-muscle capability slice, not a coupled-posture prediction.

Constant-secondary rows: Bicep Femoris, Short Head (act 10: knee_angle_r graded, hip_flexion_r robot surface constant), Vastus Intermedius (act 11: knee_angle_r graded, hip_flexion_r robot surface constant), Soleus (act 12: ankle_angle_r graded, knee_angle_r robot surface constant), Psoas X (act 27: hip_flexion_r graded, hip_adduction_r robot surface constant). For these the optimized robot
surface is constant along the secondary axis because that coordinate is not
a BiPam class hinge (the survey2 CSV notes say "remaining joint pinned at 0,
robot surface constant along that axis (same as the hand-tuned row)"); the
hand-tuned robot surface has the same property on those rows. Their meet
shares are effectively 1-D over the primary coordinate, evaluated on the
same 25-cell grid against the full 2-D human surface.

## 9. Grid, percentage, and sign conventions

- Grading grid: 5x5 INCLUSIVE uniform nodes per axis over the rom-table
  ranges (25 cells per actuator); cell values are DIRECT OpenSim evaluations
  at the nodes; cell bins are midpoint-Voronoi intervals. Percentages are
  computed on the grading grid, never on the 9x9 exports.
- 9x9 exports (`human_surfaces_9x9_20261023/`) are bilinear interpolations
  from the 1-deg OpenSim samples, for VIEWING only.
- Human torque: tau_h = sum over chosen muscles of computeMomentArm (about
  the primary coordinate) x Fmax, OpenSim native signs (moment arm = -dL/dq;
  positive torque acts in the positive coordinate direction); knee flexion
  NEGATIVE; hip adduction POSITIVE.
- Per cell: meet = robot torque has the human torque's sign AND
  |tau_robot| >= |tau_human|; active = positive strain AND taut torque > 0.
- Row averages are UNWEIGHTED means of the per-muscle area shares (the
  printed 4-decimal shares of the two row CSVs), displayed rounded to two
  decimals.
- Human surfaces are STOCK gait2392_simbody routing only (human data on
  human paths; robot/BPA paths never substituted).

## 10. Provenance and environment

- OpenSim 4.4.1 (pinned) against `gait2392_simbody.osim` (sha1
  b89062c763e6d8fa3d7651c0a4fd2820055d807c) via the actuator map of record
  (`Code/Matlab/Mesh_Optimization/Routing27_v2/actuator_map_27.json`, sha1
  a97afffeeae85687810986273db6ba0151eba1a2) and the 1-D audit of record
  (sha1 47912eb6c1e858e9e62ceefe36b8a510e052c417); provenance headers are in
  the grading CSV.
- Machine easteregg2; Python `D:/Anaconda/envs/myo/python.exe` (the same
  interpreter the surface generator used); xlsx written with openpyxl 3.1.5,
  already present in that env (verified by import in this pass; no package
  was installed).
- Assembler: `assemble_deliverables_20261023.py` (this folder) re-derives
  every cell from the three CSVs of record, asserts the counts (27/27 both
  rows; 635/40/126 cell totals; 9 suspects; manifest flag equality) before
  writing, verifies the xlsx by read-back, and writes ONLY the six dated
  _20261023 deliverables. Inputs and all v1 (20261022) files are unmodified.
- Git state at assembly time (checked with `git status --porcelain`, branch
  `KneeTestSetup_BenBo_stw`): the working tree was NOT clean before this
  assembly (tracked edits to AGENTS.md, the torque-paper build files under
  `Documentation/Reports and Papers/Knee_Torque_Test/`, and
  `tmp/paper_bridge/build_flx20_20261023/` + `figscripts/`; untracked
  BiPulley_survey2 result mats under
  `Code/Matlab/Mesh_Optimization/Results/` and other files). This pass added
  only the six dated deliverables plus the assembler in this folder. No
  commits or pushes were made; Ben commits via GitHub Desktop.


## Exclusion ruling (Ben, 2026-10-23) — the table averages 19 of 27

Eight actuators are EXCLUDED from both rows of the final table (values kept
in the per-muscle CSV with reasons; nothing deleted):

1. **Knee extensors (Rectus Femoris, Vastus Intermedius):** Ben ruled that
   the legacy hand-code patella location moment arms were calculated
   incorrectly — the literature's ds/dl (dL/d-theta) derivative form is the
   correct moment-arm calculation — and the optimized row uses the rigid
   fixed-tibia bypass, so neither row is trustworthy for the quads.
2. **Six solver-record-inconsistent rows (Erector Spinae, Semimembranosus,
   Medial Gastrocnemius, Lateral Gastrocnemius, Gluteus Maximus, Adductor
   Magnus):** optimized winners record torque in cells with zero
   positive-strain cells anywhere, with resting lengths that contradict the
   recorded strain (e.g., medial gastrocnemius rest 0.139 m on a ~0.4 m+
   route). Root cause unresolved by ruling of no-time; the re-grade pass was
   stopped and the rows excluded rather than fixed (Ben: "We don't have time
   to fix everything").

Resulting averages over the 19 graded actuators: hand 66.32% meet /
97.05% active; optimized 75.30% / 95.23%. Strictly
better 11 / ties 3 / worse 5 on meet share.
