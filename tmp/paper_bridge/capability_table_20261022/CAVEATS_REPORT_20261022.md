# CAVEATS_REPORT_20261022

Plain-language caveats for the 27-actuator capability table (capability_table_20261022.md / .tex / .xlsx) and the per-muscle rows (per_muscle_capability_20261022.csv). Every number in the table re-derives from the CSVs of record in this folder: hand_tuned_metrics_20261022.csv (hand-tuned row), BiPulley_survey_perMuscle_regraded_20261022_223448.csv (optimized row), and human_torque_curves_20261022.csv (human targets, both rows).

## 1. Covered counts for the two rows

- **Hand-tuned row: 27 of 27 covered, 0 excluded.** Its averages (72.89% meet, 98.81% active) run over all 27 actuators.
- **Optimized row: 23 of 27 built; its averages (52.00% meet, 82.25% active) run over the 16 built muscles that could be graded.** Ben asked for all 27; the missing are stated, never faked:
  - **Could not build (4):** Erector Spinae, Internal Oblique, External Oblique, Vastus Intermedius. The master muscle path crosses no joint (all routing points sit on one body), so a joint actuator cannot be formed. This is a property of the OpenSim reference paths, not a solver failure.
  - **Built but grade-unable (7):** Rectus Femoris, Tibialis Posterior, Peroneus Brevis, Peroneus Longus, Peroneus Tertius, Adductor Magnus, Psoas X. Each was built from its seed route, but the regrade spec's hinge chain does not model its primary coordinate (subtalar_angle_r for Tibialis Posterior, Peroneus Brevis, Peroneus Longus, and Peroneus Tertius; knee_angle_r for Rectus Femoris; hip_adduction_r for Adductor Magnus; hip_flexion_r for Psoas X), so no regraded percentage exists. Their opt columns are blank in the per-muscle CSV.

## 2. Transcription suspects (hand-tuned row)

Nine actuators carry a transcription-suspect flag because the corrected Connor script differs from Ben's reference values at one or more routing points. **All were graded as-is**; none was silently corrected. A suspect flag means the hand-tuned percentage could shift once Ben rules on the value; it does not mean the grade is wrong.

- **Sartorius** (act 6, hand meet 28.00%): VALUE MISMATCH: col2 vs Ben sar_wr1: Connor [[0.091 -0.105 0.037]] Ben [[0.091 -0.105 0.038]]; col7 vs Ben t6: Connor [[0.038 -0.382 -0.004]] Ben [[0.038 -0.382 0.004]]
- **Gracilis** (act 9, hand meet 12.00%): origin differs from Ben p8; Connor y z equal rect_fem origin p10's (likely transcription slip); VALUE MISMATCH: col1 vs Ben p8: Connor [[-0.074 0.125 0.079]] Ben [[-0.074 -0.119 0.028]]
- **Bicep Femoris, Short Head** (act 10, hand meet 100.00%): VALUE MISMATCH: col3 vs Ben twr3: Connor [[-0.023 -0.056 0.029]] Ben [[-0.023 -0.056 0.034]]
- **Tibialis Posterior** (act 15, hand meet 100.00%): MIF 3549 = map master_mif_n but Ben PAM_Prop has tib_post_mif 1588 (map documents 3549 as soleus copy-paste artifact; same 40 mm diameter either way); routing values match Ben
- **Peroneus Longus** (act 18, hand meet 40.00%): VALUE MISMATCH: col1 vs Ben t2: Connor [[0.005 -0.062 0.047]] Ben [[0.006 -0.062 0.047]]; col5 vs Ben per_long_wr4: Connor [[0.048 0.011 0.028]] Ben [[0.068 0.011 0.028]]
- **Peroneus Tertius** (act 19, hand meet 96.00%): three components differ from Ben (values appear shuffled between t2/per_tert_wr1/a2); VALUE MISMATCH: col1 vs Ben t2: Connor [[0.006 -0.062 0.079]] Ben [[0.006 -0.062 0.047]]; col2 vs Ben per_tert_wr1: Connor [[0.023 -0.407 0.022]] Ben [[0.023 -0.407 0.016]]; col3 vs Ben a2: Connor [[0.079 0.016 0.034]] Ben [[0.079 0.022 0.034]]
- **Extensor Digitorum Longus** (act 22, hand meet 40.00%): VALUE MISMATCH: col5 vs Ben ext_dig_wr4: Connor [[0.005 0.005 0.015]] Ben [[0 0.005 0.015]]
- **Adductor Magnus** (act 25, hand meet 76.00%): Location columns 2 and 3 both equal add_mag_wr2; Ben routing uses add_mag_wr1 then add_mag_wr2 (wr1 absent; likely transcription slip); VALUE MISMATCH: col2 vs Ben add_mag_wr1: Connor [[-0.059 -0.108 -0.03]] Ben [[-0.136 -0.093 0.059]]
- **Psoas X** (act 27, hand meet 56.00%): PLocation undefined-variable slip (values = Location4; map: presumed Location4); b1/wr1 frame roles swapped vs Ben's inline formula; point values match Ben

**Code-level note (not a per-muscle value):** the repo copy `Code/Matlab/Previous Optimization Code/PamData.m` line 271 reads `size = delta/(1 - 1/def);9hi`. The stray `9hi` token is syntax corruption: re-verified in this assembly on easteregg2 with MATLAB's Code Analyzer (`checkcode` on the repo file), which reports **"Parse error at hi: usage might be invalid MATLAB syntax." at line 271, severity error**; the hand-tuned grading stage additionally verified at runtime on eb475ws4 (MATLAB R2025a) that the class fails to load with `Unexpected 'hi'` (documented in `ht03_grade_hand_tuned_20261022.m`'s header). Ben's original `Size.m` in the old-code copy has the clean line. The run used a staged dated copy with only that token removed (diff-verified single-line change); the repo file was NOT modified. Fixing it in place is Ben's call.

## 3. Artifact-flagged human curves (caveat only)

Six human target curves carry via-point artifact flags from the path-artifact audit: **Semimembranosus (semimem_r), Rectus Femoris (rect_fem_r), Gracilis (grac_r), Vastus Intermedius (vas_int_r), Iliacus (iliacus_r), Psoas X (psoas_r)**. The flags are printed in the curves CSV's `artifact_flag_caveat` column and in the per-muscle CSV's notes. They are **caveats only: the curves were never corrected**, and grading used them as generated.

## 4. Unconverged small-budget runs

**None.** Every surveyed site's checkpoint records stopReason `converged` (the regrade CSV's converged column is 1 for all 16 graded rows). The cannot-build and grade-unable rows are structural exclusions, not convergence failures.

## 5. Budget disclosure: this is a capability survey, not a design campaign

Each optimized design was chosen with a few hundred function evaluations per muscle (graded rows used 34-195 evaluations, median 130). That budget answers "can this routing class reach the human curve at all"; it does not produce flight-ready designs. Several rows are additionally **capability-capped**: the ladder's final rung still reported the verdict `torque_short` at the directive configuration cap (nBPA, G) of (3,1), (3,2), or (3,3), meaning more actuators or tackle might still help. The capped muscles are flagged in the per-muscle CSV notes (11 of 16 graded rows).

## 6. Still pending: hip adduction/abduction coupling

The hip adduction/abduction coupling (Ben's 9 control points) is **still pending** and is not represented anywhere in this table. Any row whose actuator crosses the hip is graded with that coupling absent.

## 7. Target-change disclosure

The optimized designs were **chosen against the earlier, harder group targets** (full survey run 14:06-14:49 against the group-target curves) and then **re-graded against the single-muscle targets** that both rows now use (Ben's 2026-10-22 ruling: each actuator is graded against the same-named muscle only; no cross-muscle group sums). Optimizing against harder targets and grading against easier ones is a **conservative over-design**: the optimized percentages are, if anything, understated relative to a like-for-like optimization. The per-muscle CSV states the target version for every row in two dedicated columns: `hand_target_version` (`single-muscle curves` on all 27 rows) and `opt_target_version` (`group-targets chosen; regraded vs single-muscle curves` on the 16 graded rows, `n/a (not graded)` on the other 11); the optimized rows' notes also carry it verbatim (`optimized vs GROUP targets; regraded vs SINGLE-MUSCLE csv`). Two optimized rows also record Ben's target-change effects explicitly (Semimembranosus meet 44.0% -> 48.0%; Gluteus Maximus unchanged). Provenance of the hand-tuned row's target version: its grader loads `human_torque_curves_20261022.csv` (`tmp/hand_tuned_eval_20261022/ht03_grade_hand_tuned_20261022.m:72`) "copied fresh from easteregg2 immediately before this run" (script header), the run is time-stamped 23:05:59 on EB475WS4 (`final_run2_20261022.log:4`), and easteregg2's file of that name has carried the single-muscle content since its 21:59 regeneration.

## 8. Grading grid and the neutral-pose pinning convention

Each actuator is graded on **25 inclusive uniform points over its primary coordinate's range of motion** (the campaign `rom_table_deg` table; verified against the actuator map by the grading scripts). **Every non-graded joint the muscle crosses is pinned at 0 deg**: a one-DOF neutral-pose slice. Moment arms the human model develops at other joint configurations are therefore not represented; this is a per-muscle capability slice, not a coupled-posture prediction.

## 9. Sign conventions

Human torques come from OpenSim `computeMomentArm` times Fmax with **OpenSim native signs** (moment arm = -dL/dq, so positive torque acts in the positive coordinate direction); **knee flexion is negative** (vasti positive, hamstrings negative). A grid point counts as **meet** when the actuator is taut and its torque magnitude is at least the human magnitude with matching sign; **active force** counts taut points with positive BPA strain and nonzero torque.

## 10. One data-quality flag inside the optimized row

**Bicep Femoris, Long Head** shows meet 100.00% with active 0.00%: the winner is taut everywhere it meets but has zero positive-strain cells (flagged `SUSPECT` in the regrade CSV). The grade is included as-is from the CSV of record; treat that row's 100% with care until the strain-window question is resolved.

## 11. Environment notes

- Human curves generated with OpenSim 4.4.1 (pinned) against `gait2392_simbody.osim` (sha1 b89062c7...) via the actuator map of record (`Routing27_v2/actuator_map_27.json`, sha1 a97afffe...); provenance headers are in the curves CSV.
- The xlsx was written with openpyxl **3.1.5, already present in `D:/Anaconda/envs/myo` when this assembly ran** (verified by import immediately before writing; the earlier 23:29 assembly pass reports having pip-installed it, 3.1.5 plus et-xmlfile 2.0.0). No package was installed during this pass.
- All five deliverables were rewritten in place by this pass's assembler `assemble_deliverables_20261022.py`, which re-derives every cell from the three CSVs of record and asserts the counts (23 built, 16 graded, 4 not buildable, 7 grade-unable, 9 suspects, 6 flagged curves) before writing; the assembly modified no file outside this folder. Relative to the earlier 23:29 pass, the table row labels now state the count actually averaged ("Hand-tuned (27 of 27 graded)", "Optimized (16 of 27 graded; 23 built)") and the per-muscle CSV gained explicit `hand_target_version` / `opt_target_version` columns. All numbers are unchanged.
- **Repo-state disclosure (re-verified in this assembly on easteregg2, git CLI on PATH, `git status --porcelain`, branch `KneeTestSetup_BenBo_stw`):** the working tree is NOT clean, and this campaign folder holds two of the modified tracked files: `human_torque_curves_20261022.csv` and `human_torque_curves_manifest_20261022.json` were regenerated in place by the earlier stages of this campaign (13:12 group-target version -> 21:59 single-muscle version, per Ben's 2026-10-22 target ruling; these are inputs of record for this table, not products of the assembly). The other two tracked modifications are `Code/Matlab/Mesh_Optimization/Results/Knee_Flexor_Data_20mm.m` (`tendon = xBest(8)-.02`, `KMAX` 0.255 -> 0.235) and `Testing_Data/2026_06_Festo/readserialnumbers2.m` (`prf` [0.085, -0.296, 0] -> [0.095, -0.325, 0]); both are hardware-side edits outside this campaign's scope and their attribution could not be verified from the working tree alone. Untracked non-tmp files are also present: `Flx_20mm_pressure_for_human_torque.csv` and `Knee_Flexor_20mm.gif` at the repo root, and `Testing_Data/2026_06_Festo/Flx_20mm/`. Ben commits via GitHub Desktop; nothing was committed by the agents.
