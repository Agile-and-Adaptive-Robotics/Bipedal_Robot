> Superseded by [the current committee revision](../committee_revision_20261027/REVIEW_REPORT.md). Working source is now ProofFinal; do not follow the old source-directory or outstanding-figure instructions below.

# Dissertation checks for easteregg2

These are bounded tasks for the remaining committee-copy checks. Deadline: October 28, 2026, 9 a.m. Pacific. The source candidate is in `Documentation/Reports and Papers/Dissertation/Notes/submission_review_20261027/`. Transfer that folder or its source ZIP to easteregg2 before editing it there. Do not assume it has already been synced. Do not overwrite the live Overleaf project, original downloaded archive, or ProofFinal. Ben explicitly excluded **Table 4.3 and its Hunt comparison** from corrections.

## Task 1 Verify flexor force transformation and adopted run provenance

Start by checking hostname. On easteregg2, use `D:\GitHub\Bipedal_Robot`; follow that machine's AGENTS.md bullet. This task does not need parallel optimization.

The current `Testing_Data/2022_02_Festo/minimizeFlxPin.m` creates forces in newtons at lines 257–258, then uses:

```matlab
Fbrk(ii,:) = RowVecTrans(Tkbr\eye(4),Fk(ii,:));
Fbrh2(ii,:) = RowVecTrans(Thbr\eye(4),Fh(ii,:));
```

at lines 314–315. `Code/Matlab/Functions/ModernRobotics/RowVecTrans.m` appends homogeneous coordinate 1, so it applies translation as well as rotation. That is appropriate for an attachment point, but a three-component force vector should receive only the rotation. A full wrench needs the appropriate adjoint transformation and separate moment components; simply translating a force vector is not that operation.

Work in copied evaluator files under a new scratch directory. First determine which version actually generated `minimizeFlxPin10_results.mat`, adopted row 77, and which code path is used by `buildKneeFlexorContext20mm.m`. The header/standing context says pitch-only insertion and two-rotation origin; executable lines 289 and 306 currently use pitch-only at both. Do not infer saved-run provenance from the current header.

Then run the existing adopted parameter set through the unmodified evaluator and a copy that changes **only the force-vector coordinate conversion to rotation-only**. Keep bracket geometry, frame choices, data, fold membership, and scalar parameters identical. Do not refit first. Report maximum changes in force, effective compliance, torque, and per-test RMSE/FVU, separated by training, held-out validation, and unused test. If the deployed 20 mm prediction uses the same operation, compare its torque envelope at the fixed adopted design as well. Confirm whether any reported pass/fail margin changes. Search active equivalent implementations for the same operation, without modifying them globally.

Required output: a short result table with exact source hashes and adopted MAT identifiers, the two evaluator copies, and a recommendation stating whether the dissertation numbers remain defensible, require a documented correction, or require a fresh identification. Preserve one shared stiffness pair across flexor/extensor configurations. Do not resurrect retired result MAT files or launch a broad optimizer campaign. Do not change Table 4.3 or silently change dissertation numbers.

## Task 2 Correct figure legends without changing the evidence

For the extensor physical-validation graphic, the prepared generator is `tmp/paper_bridge/figscripts/gen_Ext20mm_cam_20261026.m`; its audited series are in the torque-physical-audit inbox's `csv/Ext20mm_series_audit_20261026.csv`. Change the embedded legend labels to **Vastus medialis reference** and **Prediction at measured pressure**. Keep the nine measurements, pressure-specific predictions, optimized curve, and human reference identical. Save only into the review candidate. Check the regenerated values against the audited CSV and the known maximum full-pressure curve difference before replacing the asset. Retain explicit measured-versus-predicted marker/style distinction and readable labels.

For Appendix A, use `Code/MuJoCo_SNS/spinal/_muscle_force_compare.py` and the archived force-comparison evidence. Replace the OpenSim legend's **compliant** claim with **OpenSim static optimization**. The existing script computes poses and matched activations without assigning trajectory velocities; it cannot support a dynamic tendon-compliance claim. Do not reframe it as a new dynamic validation. Prefer replotting archived arrays. If a rerun is required to recover arrays, compare every saved statistic with `muscle_force_compare.csv` before replacing the graphic. If they do not reproduce, report the mismatch rather than overwriting the dissertation's saved result. Increase plot text where practical without changing normalization. The current candidate caption already explains the limitation; remove the editorial warning about the inherited legend only after the legend is corrected.

For either figure, render the resulting dissertation page and inspect legend, symbols, and caption. A scientific figure correction must preserve the plotted evidence, not approximate it with generated artwork.

## Task 3 Quantify uncertainty if the records support it

Locate the calibration and geometric measurement records for the nine extensor and fifteen flexor points. Use the audited reduction and its load-cell-angle convention, not a mixture of old and new conventions. Establish uncertainty in force, angle, reaction-point position, and pressure, and whether repeated trials exist. Propagate only supported uncertainties and report their effect on the minimum nominal margins, +0.24 N m extensor and +0.42 N m flexor. If the needed uncertainty data were not recorded, return that conclusion explicitly; the candidate already limits the claim to nominal measured capacity and does not require invented confidence bounds.

## Return package

Provide a concise verdict, changed-file list, exact before/after numerical values, updated figures if verified, and any source patches in a separate inbox folder. End any dissertation revision with a fresh full build and per-directive line evidence. No commits or pushes. No messages to other chats unless Ben separately authorizes them. New locomotion searches, Hyfydy migration, and defense-deck work are not required for these checks.
