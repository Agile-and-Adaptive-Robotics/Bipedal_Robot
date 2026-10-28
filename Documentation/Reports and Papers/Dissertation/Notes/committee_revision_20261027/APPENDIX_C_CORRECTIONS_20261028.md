## Measured test dimensions and target comparison (2026-10-28)

Latest measured-dimensions/target clarification: MEASURED_DIMENSIONS_AND_TARGET_20261028.md. Table C7 now records individual resting/tendon lengths: flexors 561/28 and 560/30 mm; extensors 828/92 and 826/85 mm. Shared prediction inputs remain separate. Both manuscripts explicitly state that the -86.5-degree measurement exceeds the human target; the 2.06 N m shortfall is against the pressure-specific prediction only. All nine extensor measurements retain their meet/exceed conclusion. No data, coefficients, fit values, or curves changed. Local dissertation 321 pages/live 322; paper 33 + supplement 10. No commits/pushes.

# Appendix C drawings, dimensions, and transmission notation

Completed the drawing, table-formatting, terminology, and reverse-pulley schematic corrections in the active sources. Actual per-actuator physical dimensions remain unconfirmed; a question is pending with Ben. No commits or pushes, and no ZIP restored over working edits.

## Drawings

- Figures C1, C5, and C6 now show registered femur and proximal-tibia references without a foot. The tibial drawing is limited to the proximal shaft. The source meshes contain femur and tibia only; no foot geometry is loaded.
- Corrected the reference tibia's transform: the old helper used the human joint transform, while the route and clearance circles use the robot t1-to-ICR-to-femur transform. The new illustration uses that same robot transform.
- Presentation registration places the native femoral condyle reference [6.52, -408.00] mm at the saved femoral ellipse center, and the native proximal tibial groove reference [6.00, -31.255] mm at the saved lower tibial circle center. These are illustration registrations, not changes to the robot route, collision geometry, or measured anatomy. The proximal tibia is cropped at local y=-250 mm, excluding the distal region.
- The matching A/H route panels in the combined main-chapter extensor figure and the paper's supplementary cam panels were updated. The combined figure and hardware photograph remain intact. Saved route points, release states, clearances, cam radius, optimizer results, and torque curves were not changed.

The native MATLAB renderers load existing MAT results and only redraw. Generators and registration notes are in Figures/94-AppendixC/Components/registration_20261028; before figures and rendered checks are in Archive/appendix_c_corrections_20261028.

## Transmission notation and reverse-pulley schematic

In Figure C4, G is the ideal tendon-travel multiplier and the number of moving-block cable parts. It is not force gain. With n parallel BPAs, ideal insertion force is n*F_b/G and force gain relative to one BPA at force F_b is n/G. Thus all three displayed count/transmission pairs (1,1), (2,2), and (4,4) have ideal force gain one, while travel multipliers differ. Each configuration has separately re-solved actuator sizes, so equal transmission force ratios do not imply identical actuator forces or joint torques.

The former generic Transmission axis label now reads nBPA, G, and the caption defines the two quantities. Surface data are unchanged.

The easteregg2 schematic existed but had not been included in the manuscripts. Its captured version was crowded and unsuitable as a final figure. A new generic vector schematic is now included in both the dissertation torque methods and the paper (paper Figure 7). It shows parallel BPAs pulling the moving sheave, the two supporting cable parts, cable travel at the insertion, the force/travel equations, the fixed tackle's vector cable load, the opposite support reaction, and the separate BPA origin mount. No individual vastus optimization result was reintroduced. The schematic explicitly assumes ideal cable and rigid blocks; compliance and pulley losses are omitted.

Active shared asset: reverse_pulley_principle.pdf in Dissertation/Figures/30-torque and Knee_Torque_Test/Figures. Editable generator and alternative text are in their Components/reverse_pulley_principle folders. The earlier easteregg2 schematic and its MATLAB/FIG/EPS/PDF assets are preserved in the pass Archive.

## Tables and terminology

- Added C.1 tables for the actual saved actuator objects loaded by minimizeFlx and minimizeExt: new Tables C5 and C6. The obsolete FlxBioBPASet snapshot was inspected and rejected as the source of current input values because it differs from the files the active evaluator loads. No evaluation or optimization was rerun to derive these tables.
- Added Table C7 with four rows (two flexors, two extensors), explicitly reporting shared model inputs per actuator. Flexors: resting length 561.49 mm and effective tendon 5.00 mm (25 mm route parameter minus 20 mm fitting adjustment). Extensors: resting length 827.74 mm and approximate tendon 92.49 mm. These are prediction inputs, not four independently measured installed-actuator dimensions. The flexor test plot's nominal 562 mm label is retained.
- Original Table C7 is now C10. It separately lists optimizer variables, optimizer bounds, and physical mappings: x0 in [0,2] corresponds to chi0 in [0,20] mm; log stiffness coordinates map through 10^x; the extensor's [-2,0.5] coordinate corresponds to [-20,5] mm. The physical bounds and stiffness lock are unchanged. The caption now agrees with the columns.
- Original Tables C8/C9 are now C11/C12. Fit and Held out are replaced with Training and Validation. Existing values, fold membership, unused-test classifications, and all fit statistics are unchanged.

### Physical-dimension limitations and metadata differences

The final-test worksheet has empty per-actuator resting-length cells; separate actual-length cells are also empty. Template tendon entries elsewhere lack adequate measurement/unit provenance and were not promoted to actual final dimensions. The extensor analysis itself calls its tendon length approximate. Ben was asked for separate measured resting/tendon lengths for all four installed BPAs. That request remains open.

Two existing discrepancies were identified and disclosed rather than silently reconciled: the active 10 mm flexor object has 11 mm tendon versus 21 mm in the physical description; the extensor object has 520 mm resting length versus 518 mm in the description. New tables label these as saved evaluator inputs. This pass does not change either the existing prose dimensions or the model objects, predictions, coefficients, measurements, or fit results. The older FlxBioBPASet snapshot also differs from the actively loaded objects and must not be substituted for them.

## Verification and delivery

- Local dissertation: fresh five-pass build succeeds, 321 pages, no errors or undefined references. One inherited 0.669-point overfull line remains. The new tables and figure fit within their pages; the new final-dimension table is kept before its explanatory paragraph.
- Local paper and supplement: all compilation/BibTeX steps succeed, 33 main pages and 10 supplement pages; no undefined references. Existing template/font and small header/footer overflow messages remain.
- Live dissertation: 322 pages, Errors 0 / Warnings 0. Live paper: 33 pages, Errors 0, with the existing 19 font warnings. No undefined-reference warning was found. The supplement has a verified local build; it was not selected as a separate live main document in this pass.
- Both dissertation sources and paper Main.tex were read back from Overleaf and match local sources after normalizing line endings. Auxiliary files were refreshed for supplementary cross-references.
- Delivery PDF/ZIP files updated. Source bytes and ZIP integrity checked. Dissertation includes 90 active graphic files. Conservation checks show all existing main-chapter and paper table blocks unchanged, all prior Appendix C table numbers unchanged except the explicitly re-expressed bounds, and sixteen other existing Chapter 30 graphic files byte-identical. Only the saved route illustrations were redrawn. Protected Chapter 20 prose and figures match their pre-pass hashes.
- No tracked changes were accepted/rejected and no review-comment dispositions were changed. Older ZIPs, before sources, original figures, and the three locked prior paper delivery files are preserved in Archive.

## Directive evidence

| Request | Working evidence |
|---|---|
| No foot near p9; align femur and tibial groove | Figures C1/C5/C6; appendix_c_bone_background.m and registration.txt in Components; native rendered before/after assets in Archive. |
| Clarify G in Figure C4 | fig:PulleyMargins_resolved caption and nBPA, G axis label; Chapter 4 and paper reverse-pulley schematic. |
| Biomimetic actuator tables in C.1 | tab:app_flxbio_tests and tab:app_extbio_tests; source_constants.json and audit_appendix_c_inputs.m. |
| Final two-flexor/two-extensor dimensions | tab:final_test_dimensions; explicit model-input caption and still-pending actual-measurement request. |
| Optimizer versus physical bounds | tab:app_opt_bounds, now Table C10. |
| Training/validation wording | tab:app_gof and tab:app_gof_ext_all, now Tables C11/C12. |
| Reverse-pulley schematic | fig:reverse_pulley_principle in both manuscripts; active PDF and reproducible vector generator. |

Builds and detailed verification are archived outside ProofFinal after delivery. The single remaining data-dependent item is confirmation of actual per-actuator final-test dimensions, including resolution of the recorded/evaluator metadata differences above.
