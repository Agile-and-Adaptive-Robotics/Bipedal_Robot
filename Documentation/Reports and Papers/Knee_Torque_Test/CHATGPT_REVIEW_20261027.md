# Torque paper and dissertation review, 27 October 2026

## Result

The paper and dissertation torque material have been revised locally and rebuilt. The paper is 41 pages, the supplement 8 pages, and the dissertation 320 pages. All three final logs have zero LaTeX errors, zero undefined references or citations, no duplicate labels, and no oversized floats. Existing font and overfull-box warnings remain in inherited material; this is not a claim that every page of the dissertation has received a fresh scientific or typographic audit.

No optimization was run, no MATLAB process was started, and no study script, result mat, or experimental dataset was changed. Nothing was committed or pushed. The subsequent Overleaf upload is documented at the end of this report. The old build directories and original figure assets remain available.

## Main changes

- Restored the paper's step-by-step local force, Hooke-law displacement, component displacement, projection, effective stiffness, and series-compliance development. The derivation now states the assumptions used by the scalar reduction.
- Redrew both schematics as vector artwork, with editable SVGs, PDF assets, PNG previews, and alt text. The projection explains off-axis displacement and series shortening. The wrapping figure depicts a braided BPA crossing an actual flexed-knee schematic, with femur, tibia, fittings, tendon, R, and the full contact angle alpha. In the dissertation these are separate Figures 4.5 and 4.6; in the paper they are Figures 5 and 6. The diagram is explicitly schematic and does not present its drawn dimensions as measured geometry.
- Renamed the dissertation subsection to **Identification of the Torque-Model Correction Terms**. Removed the obsolete 48.0 cm exception and 47.9 cm exclusion narrative from the chapter's identification/results exposition. The appendix still preserves the historical test inventory and exclusion record. The adopted folds and shared stiffness pair are unchanged.
- Corrected the claim that every stiffness array is y-symmetric. The actual axis assignments differ by configuration. Both pinned-flexor frames execute the one-rotation construction in the current evaluator, despite an older file-header description; the appendix already records the executed assignments. No stiffness parameter or evaluator was changed.
- Corrected the original placement equation against **costFunctionKnee.m**: minimize negative torque surplus plus alignment mismatch, with the torque contribution divided by **3N**. The earlier equation's maximized mismatch term would reward the wrong direction. Added the original method's term/weight/meaning table, retained the later cost-function evolution, and removed the Morrow placement-paper citation from Background, placing it with the dissertation work in the torque chapter.
- Removed the incorrect interpretation of RMSE, FVU, and maximum residual as separate measures of level, shape, and phase. The definitions remain explicit. Removed the unsupported assertion that the bracket solve necessarily reaches adjacent floating-point values; the text now describes its actual default displacement tolerance. See [MathWorks fzero documentation](https://www.mathworks.com/help/matlab/ref/fzero.html).
- Reframed the bending correction as an empirical usable-length loss rather than a derived braided-wall stress model. Retained the correction equation and fitted values.
- Reconciled the capability comparison to the saved 27-site results, added the axis-average comparison, and added the site/DoF/diameter inventory to Appendix C and the supplement. The diameter column reports the inherited robot sizing rule used in the saved survey, rather than claiming a new diameter extraction from the OpenSim XML. The coordinate columns are the two graded coordinates, not a complete anatomical crossing inventory.
- Restored the contact/release derivation to the paper supplement, linked it from the paper, and corrected its flexor-route caption from mirrored to staggered. Removed the stale causal 4.51 N m redirect claim still present in Appendix C, which Ben had already rejected in the paper.
- Repaired the paper's broken concluding sentence, corrected the unsupported “two thirds or more” claim in its abstract and introduction, and aligned the conclusion with the measured 15 flexor and nine extensor points. Removed another stale statement that the results figures use the initial fit; their captions and curves use the adopted two-stage terms.
- Preserved all 17 main-paper display items, including both original frame figures, per Ben's preference to consolidate later. There are 12 figure floats and five table floats. No figure was removed to enforce the journal ceiling.
- Made the detailed torque-weight table paginate, wrapped Appendix C tables within the margins, reduced an oversized background figure slightly, and repaired a database table whose natural-width columns previously produced a float over 1,500 pt taller than the page. These last two edits change layout only.

## Numerical changes and provenance

No recorded measurement was altered. The material replacements in exposition are explicit here; the complete session diffs accompany this report.

| Item | Previous exposition | Revised exposition | Evidence |
|---|---|---|---|
| Paper capability summary | 19 sites; 66.32/75.30 percent coverage, 97.05/95.23 percent active force | 27 sites; 37.28/65.19 percent cell coverage; 38.39/66.17 percent axis mean; -1.93/-3.63 mean torque-ratio score | `capability_v3_three_metric_20261027.csv` |
| Per-site coverage comparison | 11 wins, three ties, five losses on the old 19-site assembly | 23 wins and four ties on the current 27-site assembly | Same saved CSV |
| Actuator-count conclusion | Claimed a minimum three-actuator scale for the full envelope | No tested combination meets all five angles; the best saved combination meets four. The study does not establish a minimum count that meets the full envelope. | `VasMed_minBPA_study1_summary_20261027.csv` and `study1_vasmed_canon_20261027.m` |
| Routing comparison | Read as an attachment-geometry improvement | Same underlying master's attachment routing, with transmission/sizing changes and different baseline diameter substitutions. The updated side still uses nominal 40 mm actuators at eight sites. | CSV headers and `study2_capability_v3_20261027.m` |
| Meaning of negative mean score | Described as deeper torque shortfall | Negative clipped torque ratios represent wrong-sign torque; inactive cells contribute zero. | Saved metric definition in CSV header |
| Axis averaging | Not stated in the paper | Equals cell averaging on complete rectangular grids; differs only on four clipped hip domains | Saved CSV's M1, M2, and Fubini checks |
| Initial fit narrative | Exceptions referring to 48.0 and 47.9 cm tests in the adopted identification exposition | Actual adopted four-test flexor and six-test extensor folds; history remains in Appendix C | Evaluators, fit tables, and appendix |
| Unity-gain causal claim | Attributed an extra 4.51 N m to the redirect | Removed this causal interpretation; retained the independently optimized run outcomes | Ben's 26 October audit ruling and saved REVIEW.md |

The zero-correction comparison, fitted parameters, physical-validation points, 5 percent design requirement, and cam design-of-record values remain unchanged. The former 40.2 N m at 804 N clause was **not** restored. The stale pA/pB-swap material in Revision_C was **not** applied. The current convention remains **pA = origin, pB = insertion**.

## Directive compliance

| Directive | Source evidence |
|---|---|
| Restore bracket-deflection development in the paper | `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:249 |
| Restore an independent projection figure | `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:315 |
| Separate the knee-wrapping figure, with radius and full alpha angle | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:293 |
| Rename and correct the identification subsection | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:296 |
| Define the Xi objective vector and its components | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:317; `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:361 |
| Keep pA = origin and pB = insertion | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:101; `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:205 |
| Explain the original Morrow/Bolen cost and its weights | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:328; `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:335 |
| Remove the placement-paper citation from Background | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:325 |
| Preserve evolution of stiffness-aware redesign costs | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:365; `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:375 |
| Include contact/release calculation in both documents | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\94-AppendixC.tex`:256; `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\frontiers_SupplementaryMaterial.tex`:281 |
| Use saved actuator-count results without rerunning optimization | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:685 |
| Compare cell coverage, axis means, and torque-ratio score | `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:606; `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\30-torque.tex`:699 |
| List all 27 sites, graded DoFs, and diameters | `Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\94-AppendixC.tex`:392; `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\frontiers_SupplementaryMaterial.tex`:337 |
| Retain all frame figures in the main paper | `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:224; `Documentation\Reports and Papers\Knee_Torque_Test\ProofFinal\Main.tex`:231 |


## Verification and review scope

The scientific/editorial review covers the torque paper, supplement, dissertation torque chapter, and relevant Appendix C material. The full dissertation was rebuilt, and the two identified layout defects in Background and the database chapter were repaired; this does not certify all other chapters' claims or references. Every pre-session label in the edited sources is still present. Bibliography generation succeeded with absolute BIBINPUTS/BSTINPUTS, without rewriting thesis.bib. New figures and affected equations/tables were inspected in rendered PDF pages. No optimizer was used for verification.

Build hashes and warning counts are recorded in `verification.json`. Source snapshots are preserved under `tmp/torque_review_20261027/*_before.tex`; the `.diff` files in this review folder compare against those exact session-start files, so they separate this review from ZCode's existing uncommitted edits.

## Files for review and Overleaf

The sources in ProofFinal are the edited sources. The paper review build lives separately under `Knee_Torque_Test/build_review_20261027/` because the original in-place supplement log was locked. The dissertation review build is `ProofFinal/build_r12_torque_review_20261027/`.

`Overleaf_changed_files_20261027.zip` contains only changed sources and the two new figure sets, grouped by project. Map them to the existing project folders before uploading; the archive does not replace either complete project. The frame figures remain in the main paper, so a later submission pass can decide how to consolidate the 17 display items.

The reproducible schematic generator is `Knee_Torque_Test/Figures/Figure components/draw_torque_schematics_20261027.py`. No existing asset was overwritten. Ben can commit the six changed source files, the new figure sets in both figure trees, the generator, this report/review bundle, and the new review builds as desired. Existing unrelated changes in the working tree belong to earlier work and were left in place.


## Overleaf synchronization, 27 October 2026

Uploaded the revised sources and schematic PDFs to both user-selected projects and recompiled them:

- [Torque paper](https://www.overleaf.com/project/622e9536b8ea1e0a7f167cf5): Main.tex, frontiers_SupplementaryMaterial.tex, and both new schematic PDFs in Figures/. The manuscript compiles to 41 pages with **0 errors** and **0 undefined-reference warnings**. The 19 remaining warnings are font substitutions; typesetting messages also remain. Refreshed Main.aux and frontiers_SupplementaryMaterial.aux because older manually uploaded copies masked the new supplementary labels. The project's main document is left as Main.tex.
- [Dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb): chapters/15-background.tex, 30-torque.tex, 50-database.tex, 94-AppendixC.tex, and both new schematic PDFs in figs/30-torque/. Added the missing figs/30-torque/ fallback to Appendix C's graphicspath, both locally and on Overleaf, so its existing pulley and cam PDFs resolve there. The dissertation compiles with **0 errors** and no missing-figure errors. Its live project also contains Appendix E, which duplicates three database labels from Chapter 5: fig:sadb_explorer, tab:sadb_afferents, and tab:sadb_pathways. Overleaf reports these three conflicts plus their summary warning (4 warnings). Those other live appendix contents were left intact.

No optimization was run. No comments were resolved, and no existing tracked changes were accepted or rejected. The browser file helper lacked file-URL access, so file uploads used the standard Windows picker. Native mouse/keyboard control was stopped when Ben requested his screen back; subsequent checks used the Chrome extension.

Build screenshots are saved under tmp/torque_review_20261027/overleaf/. The changed-files ZIP and Appendix C session diff include the Overleaf path correction.
