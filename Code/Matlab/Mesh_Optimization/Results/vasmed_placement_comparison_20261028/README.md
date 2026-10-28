# Evaluator checkpoint: placement comparison stopped

## Outcome

The actual saved two-BPA baseline reproduces exactly. The constructed G=2 class control passes its scalar closure test. G=1 equivalence fails. Therefore the question of a shorter BPA or a smaller package under matched placement optimization remains unanswered. This is an implementation discrepancy, not evidence against the reverse-pulley transmission. No optimization was run, and no candidate-comparison plots or design claims are supplied.

All new diagnostic files are in this folder. The previous result folder, canonical classes, identification coefficients, measured data, and manuscripts were not edited. No commits, pushes, or messages to other chats. Pre-existing working-tree changes were left alone.

## Compact checkpoint table

| Check | Actual configuration / scope | Result |
|---|---|---|
| Saved baseline | 2 parallel 20 mm BPAs, 620 kPa; rest 0.8277378599107651 m; tendon 0.09249 m; saved cam context | PASS: max absolute difference from stored predBest.TorqueZ = 0 N m on all 100 saved poses |
| Saved torque range | Canonical predictKneeExt20mm -> MonoPamDataExplicit_balanceX3 | 40.2035989333 to 105.904702045 N m |
| G=1 old alternative vs canonical | One BPA on both sides, same saved geometry, pressure, Xi, and tendon length; no target enters force calculation | FAIL: max relative force error 0.08560681497436577; max absolute error 73.83267008147448 N; max effective-length error 0.01170209562813351 m |
| G=1 aligned projection diagnostic | Canonical Pbr, force-scaled bracket transform, and 18-strand-multiplier stiffness | FAIL: relative error 0.09439804745849593; absolute 81.41477867650497 N; effective-length error 0.01288798772109201 m |
| G=2 actual BiPam class control | Fixed five-row route, two crossings, 2 BPAs, active distal gain 2; 4 x 4 x 2 x 2 transform array | PASS for tested closure: force and contraction residuals exactly zero; all four states taut and feasible; slack/infeasible flags match |

Gate tolerances declared in checkpoint.m: relative force 1e-3 with denominator max(1 N, abs(reference force)), effective length 1e-6 m, and torque 1e-3 N m. G1 is not close to passing. At its maximum relative error, -119.999999863 degrees, alternative single-BPA force is 788.629858206525 N versus canonical 862.4625282879995 N. Force alone produces a 4.743765381732729 N m maximum torque discrepancy when both use the SAME canonical moment arm; this is an isolation calculation, not a validated alternative torque curve.

The previous saved_as_is_1_1 designation is incorrect: it changes BPAcount from 2 to 1. Here that variant is used only for the one-BPA regression diagnostic, never as the saved baseline. The actual baseline is always the saved count 2.

## Exact failing implementation and minimum required correction

The old rp6 psiSolve (lines 204-238) treats solved contraction s as the force-effective contraction, using one rigid-geometry Xi3 loss. Its bracketCb (lines 239-272, particularly 269) point-transforms a unit direction without the canonical pre-force scaling/normalization; its centroid and spring convention also differ.

The canonical MonoPamDataExplicit_balanceX3 pipeline instead:

1. Calculates initial Xi3 loss in updateStiffnessGeometry (lines 321-326).
2. Solves deformation in Lok/fortz (lines 529-719), applies a three-dimensional bracket displacement to route rows, and recomputes the actual segment lengths (lines 338-348 and 594-613).
3. Recalculates force-effective strain from that deformed route and tendon stretch (lines 352-365; Contraction_k at 442-525). The loss used here is not necessarily the stored obj.delta_L, which retains the initial pass.

The force-effective loss reconstructed from rest*(strain_f-strain_p) differs from the initial loss by as much as 0.01275111516543816 m. The apparently zero initial-vs-final loss check in the first diagnostic was invalid because it compared against the stored initial field; the final checkpoint corrects that check. The class has a LOCAL normalize helper at lines 866-873 performing row-unit normalization. MATLAB's global normalize is irrelevant inside that class. Neither initial suspicion should be treated as a confirmed normalization bug or a zero wrapping-loss difference.

The smallest safe correction for the G=1 adapter is to call predictKneeExt20mm / MonoPamDataExplicit_balanceX3 directly and consume its final force, strain, effective length, and torque, without reconstructing them from psiSolve. For G>1, an adapter must preserve the canonical staged deformation, route-length recomputation, final Xi3 contraction calculation, and the chosen shared/parallel tendon convention while introducing tackle travel and reaction compliance. Merely changing Pbr or kSpr does not pass. A complete minimal scalar-equation correction has NOT been established within this checkpoint; no coefficients were tuned to conceal the discrepancy. Do not replace the canonical direct-drive evaluator or existing biomimetic results with this alternative.

## Target discrepancy

- The saved ctxCamBest target is vas_med_r, and exactly matches the present OpenSim_Vasti_Results.txt column 5 at the comparison poses. Range 20.247835026..50.32136527 N m.
- The newer surface manifest identifies gait2392_simbody.osim, OpenSim 4.4.1, muscle vas_med_r, Fmax 1294 N, native knee coordinate, and torque = computeMomentArm * Fmax with all other DOFs at defaults. Range on the saved grid 30.1797690408..59.912901 N m.
- Both use negative flexion over -120..+10 degrees. The legacy text header says inDegrees=no, but its actual angle column is numerically degrees; no radian conversion is applied. The saved builder likewise consumes those values directly as degrees.
- target_comparison.csv compares every saved pose, plain and 1.05 enlarged. The maximum pointwise plain-target difference is 25.1424429611 N m. A 5% enlargement does not explain it.
- target_force_arm_audit.csv compares the corresponding legacy moment-arm file and new OpenSim moment arms. The maximum arm difference is 0.000277237641891 m. Applying the SAME constant 1294 N to the legacy arms differs from the newer torque by at most 0.358747748606 N m. Thus route-arm differences cannot explain the 25 N m discrepancy.
- Legacy torque divided by its supplied moment arm implies a varying 643.153479548..1292.69890204 N force rather than constant Fmax. This is a numerical inference from the supplied files, not recovered original OpenSim generation provenance. The exact source model, activation, and force-length/passive-force settings that produced the legacy text are not established here. Do not silently substitute a fixed-Fmax requirement for it.
- Comparison CSVs use linear interpolation explicitly; the saved cam result script uses pchip for its margin calculation. The saved target and evaluation grid coincide to numerical precision, but this checkpoint does not assert interpolation equivalence between arbitrary dense grids. Both interpolation conventions must be retained explicitly in any resumed campaign.

## G=2 control scope and limits

checkpoint.m calls the actual BiPam_pulley2_20261022 constructor. Both dimensions of its transform array are handled explicitly: third dimension has two orientations, and fourth dimension contains the two joint transforms while also defining the second grid extent in gridSize. This yields four grid cells. Outputs are FsingleBPA [2 2], Torque_p [2 3 2 2], and TransformationMat [4 4 2 2]. gate2_actual_class.mat preserves the object, inputs, configuration, and independently rooted closure outputs.

The root uses independently reconstructed bracket coefficients, but consumes canonical MuscleLength and deltaL outputs. This verifies scalar equilibrium, force scaling, effective contraction, and taut-state flags on this fixed control, not an independent route/span adapter, torque reconstruction, Xi3 handling, or compliant mount model. Mount compliance is off only in this algebraic control. It is not an admissible optimization model for the requested compliant arrangement. gate2.csv includes delivered force, which also agrees exactly. All control states are feasible and non-slack; clamp-boundary coverage is not claimed.

## Poster inspection

Components/poster_inspection.png is a rendered inspection of the supplied poster. It depicts two parallel BPAs and two pulley blocks with four wheels each; its text describes increased tendon travel and decreased transmitted force. The image does not establish the precise reeving, number of load-bearing parts, termination locations, pulley dimensions, or full package envelope. No generic replacement schematic was drawn or asserted as the actual arrangement. The requested final schematic and packaging comparison remain conditional on a valid evaluator and resolved reeving.

## Files to use and repeat

- checkpoint.m: editable, self-contained checkpoint driver plus explicitly copied old psiSolve/bracketCb helpers; all repository paths derived from this file.
- checkpoint.txt and console.log: final completed fresh-process MATLAB R2025a U1 outputs.
- saved_baseline.csv: all 100 stored and replayed two-BPA torques, force, strain, and effective length.
- target_comparison.csv and target_force_arm_audit.csv: pointwise target, force, and arm audit.
- gate1.csv and gate1_aligned.csv: per-pose failures and aligned-parameter residuals.
- gate2.csv and gate2_actual_class.mat: actual-class fixed-route control.
- checkpoint.mat: canonical predictions, old closure result, saved inputs, and targets.
- INPUT_SHA256.csv and OUTPUT_SHA256.csv: input and code hashes, and diagnostic export hashes.

Run from the repository root:

    matlab -batch "addpath('Code/Matlab/Mesh_Optimization/Results/vasmed_placement_comparison_20261028'); checkpoint"

No solver settings, seeds, incumbents, placement winners, dense finalist replay, package advantage, or publication figures exist from this checkpoint because the required equivalence gate failed. The saved baseline replay is on its original 100 poses, not a fresh dense-grid finalist verification. The exact inherited-routing extension remains unverified. No final LaTeX design-result fragments are supplied because they would imply completion of a comparison that did not run.
