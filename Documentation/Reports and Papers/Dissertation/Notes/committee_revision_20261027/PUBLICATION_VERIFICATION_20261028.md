# Publication pass verification — October 28, 2026

Local editing and verification are complete except for the externally produced, verified vastus medialis attachment-placement example and actual-system reverse-pulley schematic. That conditional result has not been received. The paper supplies the motivation and scope, without inventing a result. The dissertation explicitly identifies the saved five-node study as fixed-route sizing, reports its saved outcomes as provisional, and requires dense force/slack/feasibility grading before a design verdict.

## Source authority and preservation

The 05:14 paper and 05:16 dissertation downloads remain preserved under Archive/publication_pass_20261028. They were not extracted over subsequent edits. Working sources remain in their established ProofFinal folders; active figures remain in the sibling Figures folders. All 19 dissertation source files and both paper source files are packaged byte-for-byte from the working sources. Bibliographies were not rewritten. Protected Chapter 20 and all 15 corresponding supplied graphic/cache files match the authoritative dissertation ZIP exactly (protected_chapter_verification.json). No commits, pushes, or Overleaf uploads occurred.

## Comments and tracked changes

All 28 new comments were read, including expanded text. Twenty-seven are addressed locally; P26 remains conditional on ZCode's verified placement result and actual-system schematic. publication_comment_compliance.json contains the full comment, disposition, and final line evidence for every item. All 25 tracked-change entries were reviewed (16 Main, seven supplement, one bibliography, one latexmkrc). None was accepted or rejected. The 28 live comments remain open because the revised local sources have not been uploaded. Earlier zero-open-comment statements describe the preceding comment-only pass, not these 28 new comments.

## Numerical changes and coefficient audit

No measured data or adopted correction coefficients were altered, and no new parameter fit or optimization was run.

The old dissertation biomimetic-extensor paragraph reported RMSE 0.8920 N m, FVU 0.3396, and maximum absolute residual 1.4188 N m. The paragraph now reports 1.1702 N m, 0.6681, and 2.2088 N m, matching the existing current tables and a fresh minimizeExt.m evaluation at adopted row 32. Its coefficients are approximately [-0.006378022 m, 39978.73 N/m, 14733.75 N/m, 0.1582887]. The old Figure 10 generator instead called [-0.0109 m, 102000 N/m, 10600 N/m, 0.123]. A separate check of those older generator coefficients under the current evaluator gives [1.5393626659 N m, 1.1560479552, 2.8599739318 N m], so it does NOT establish the provenance of the three old prose values. Their exact historical source remains untraced. The figure was regenerated at the adopted coefficients; its points are unchanged. The report does not conflate historical evaluator versions with current reproducibility.

The requested baseline GoF additions are disclosed: 10 mm biomimetic flexor [7.4512924977 N m, 5.5611939073, 13.3163103416 N m], and 20 mm full-pressure flexor [14.8064358082 N m, 3.3672028448, 23.7218255272 N m]. They were computed from the canonical stored baseline M at measured angles using the evaluator's linear interpolation. Partial-pressure NaN baseline was excluded. These six values were added to the paper, torque chapter, and Appendix C tables; optimized values were retained.

Fifteen gold flexor points were newly computed at the recorded angles and pressures using the saved redesigned route, without refitting. Their provenance is flexor_pressure_publication.json and the editable MATLAB source beside the figure. They are model predictions at measured pressure, not a hybrid calculation using measured BPA length and moment arm. The freshly calculated 620 kPa curve matches the existing plotted curve exactly (maximum difference 0 N m). The measured blue points and human reference were copied unchanged. The physical-figure caption's nominal 56 cm was clarified to the already reported actual 56.2 cm. The extensor plot already contained nine gold pressure-specific predictions; these were retained and its existing dissertation export supplies the correctly named vastus-medialis reference in the paper.

The shortened paper segue removes the old 27-route implementation-status narrative, including its 19/27 first-pass and later 12 failed/15 finalized/eight feasible/six passing counts and associated margins. This is an editorial removal under the shorter-segue request, not a changed optimization result or a new success claim. Original text is preserved in the source archive. Redundant test lengths, repeated rounded GoF values, and equation algebra tokens also change in the textual conservation diff. content_conservation.json records every added/removed numeric token. No existing labels or figure floats were removed: paper 12 figures, supplement four, dissertation torque chapter 14, and Appendix C six.

## Figure and prose audit

Torque axes in the active torque-paper/chapter/Appendix-C figure set use exactly “Torque, N m.” Torque table captions use N m. Mathematical multiplication remains in equations. Figure 9's y label clears panel D; Figure 10's panel-B legend is southwest; the physical-flexor legend moved up/right and includes gold pressure predictions. The published muscle chapter was excluded from all label changes. Replaced assets were moved into Deprecated/publication_20261028 after checking all current source references. Native FIGs, generators, and numeric JSON evidence for the new plots are in Components/publication_20261028.

Prose changes distinguish Hunt/improved-BPA comparisons, approximate hybrid calculations, fitted attachment compliance, the 10 mm extensor's 25 mm wrapping radius, force/moment-arm sensitivity, cam-assisted extension, and the inherited-route survey. Orientation/transform, Hooke relations, and the three GoF metrics have separate numbered equations. The permitted bracket-series summary remains grouped. The dynamic-test discussion explains controlled-speed and controlled-load testing and the less-slack motorized rig proposal. Resting-length test labels remain in cm, tendons in mm, and geometry/parameter calculations in SI units with explicit units. First-person dissertation wording and its fuller placement methods were preserved. Two mojibake dashes in paper Author Contributions were repaired without changing contributor assignments.

## Builds and visual verification

Fresh MiKTeX/BibTeX builds completed with zero errors, zero unresolved references/citations, and no missing figures: dissertation 302 pages, paper 32, supplement 10. Paper layout polish removes a largely empty wrapping-figure page, excessive paragraph stretching, and an orphaned supplementary geometry heading. Paper still has 12 figures and four main tables; no result was removed to enforce a 15-item limit. Logs retain the Frontiers title/footer box and font-substitution warnings, and the dissertation's preexisting 0.669 pt overfull line. This is not a zero-warning claim.

All paper and supplement pages were rendered and visually reviewed; selected dissertation torque figure/table, equation, and Appendix-C pages were inspected. This is not a complete reread of all 302 dissertation pages. Final hashes, build results, and visual-review page scope are in delivery_manifest.json and final_build_verification.json. Build products are archived outside ProofFinal. Root delivery PDFs/ZIPs now contain this revised local edition; the original downloads are preserved separately. Live Overleaf still awaits Ben's drag-and-drop.

## Delivery locations

- Paper sources: Knee_Torque_Test/ProofFinal/Main.tex and frontiers_SupplementaryMaterial.tex.
- Dissertation sources: Dissertation/ProofFinal/chapters/30-torque.tex and 94-AppendixC.tex.
- Paper graphics: Knee_Torque_Test/Figures, including the active *_publication.pdf assets referenced directly by the sources.
- Dissertation graphics: Dissertation/Figures/30-torque and Figures/94-AppendixC.
- Root PDFs/ZIPs: Bolen_Dissertation.pdf/.zip and Design_Factors_for_Isometric_Robotic_Knee_Torque.pdf/.zip, plus the separate Supplement.pdf.

No verified ZCode result is included. ZCODE_VASMED_REQUEST_20261028.txt remains the governing request; integrate its verified result only after it is returned.
