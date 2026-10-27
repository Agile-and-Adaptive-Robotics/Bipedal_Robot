# ZCode: torque replot and reconciliation requests

Ben requested these changes in the live torque-paper/dissertation reviews. Read `../ChatGPT_reportTorque.md` and the saved live sources first. Do not alter the published isometric-force chapter. Do not upload whole chapters, resolve comments, accept/reject tracked changes, or launch new optimization campaigns. Replots should use existing data and the adopted correction terms. Ben commits.

## P1 — Correct the physical torque data/series before styling

Paper `fig:Ext20mm` uses `Figures/Ext20mm_cam_20261023.eps`; Ben says blue measurements occupy positions expected for the yellow back-calculated series. The dissertation `fig:extensor_physical_final` uses `extensor_20mm_physical.png` and has nine measurements, versus the paper's six. Do not fix this by swapping colors. Export a table containing angle, load, load-cell angle convention, pressure, reaction-point geometry, measured torque, predicted torque at measured pressure, 620 kPa prediction, and human target. Classify pressure-specific force-model predictions separately from a hybrid calculation using measured muscle length/moment arm.

Generator lead: `tmp/paper_bridge/figscripts/gen_Ext20mm_cam_20261023.m`, especially lines 41–137 (six-point arrays/reduction) and 140–159 (blue scatter). Current source lead: `Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.m`, especially lines 582–653: `K_deg` contains nine values. Data: `Testing_Data/2026_06_Festo/Results_table_20mm.xlsx`, `ExtTest20mm_2`. Follow AGENTS.md's load-cell convention and the current measured rig constants. Audit old flexor-template sign-flip machinery rather than assuming it is authoritative. Recalculate both documents' claims from the same declared dataset. Replace “holds” with data measurements in descriptive prose/captions.

## P1 — Hunt equation, fit stages, and GoF tables

Ben questions the Hunt equation in paper `fig:FlxPin`; the dissertation `fig:FlxPin_preopt`/pre-fit table carry the same issue. Generator lead: `tmp/paper_bridge/figscripts/gen_FlxPin_allmodels_20261023.m`, lines 55–72. It interpolates `ForceStrainForFit.mat` after replacing NaNs in z with zero, rather than directly calling the equation. Cross-check against `Code/Matlab/Functions/HuntEq.m`, correct diameter, pressure units, strain versus normalized strain, maximum contraction, hysteresis setting, force direction, and joint moment sign. Compare the source equation with the lookup table before declaring either wrong. Export the valid-point mask and report both per-model-domain and common-domain scores if the domains differ.

Ben rejects the initial flexor paragraph (paper line 328 / dissertation line 375) and initial extensor paragraph (paper line 368 / dissertation line 425). Use `minimizeFlxPin10_results.mat` row 77 and `minimizeExt10mmX3_results.mat` row 32 as adopted records; distinguish them from explicitly historical single-bracket fits. Follow AGENTS.md's fold bookkeeping: two training/two validation/one unused flexor test, and three training/three validation/three unused extensor tests. Map nominal test IDs to the actual resting lengths printed in the manuscript. Recompute scores and memberships before replacing prose, not the reverse.

## P1 — Confirmed fig:FlxGrp layout

Ben explicitly confirmed in this chat: use a **2×2 tiled layout**, three data panels titled **10 mm**, **20 mm**, **20 mm**, and a dedicated legend in tile 4. Preserve A/B/C margin letters. The last panel is the redesigned two-BPA route; state that distinction in the caption despite the repeated 20 mm title. Keep measured points, corrected predictions, original/improved-force baselines, reduced-pressure points, and human targets correctly identified. Use shared axes where comparison is meaningful, but do not clip data merely to share a range. Deliver one coherent asset for each document, not the dissertation's tiny pasted C inset.

Generator lead: `tmp/paper_bridge/figscripts/gen_FlxGrp_3panel_20261023.m`, especially layout lines 184–277. Paper asset: `FlxGrp_20261023.eps`; dissertation currently combines `FlxGrp.eps` and `Flx20mm_meetsHuman.pdf`. Preserve all fifteen verified redesigned-flexor measurements and their actual source/pressure metadata. The 2×2 instruction was posted on the χ₄-table thread, but its target is confirmed as this figure.

## P2 — Routing geometry and supplementary candidate

Paper `fig:ExtensorRoute_redesigned`: reduce the oversized eight-pose diagram; add recognizable femur/tibia outlines; Ben permits fewer representative poses. Generator lead: `tmp/paper_bridge/figscripts/gen_ExtensorRoute_A4.m`. Select poses that explain contact/release behavior, rather than retaining eight merely to fill a grid. Maintain ≥10 pt final text, meaningful dimensions, and a clear legend. The dissertation's same-named figure uses a different composite/hardware photograph; reconcile purpose before replacing it.

Ben's conditional supplementary-material comment anchors to `fig:RouteVariants`, not the routing-geometry figure. Count the final floats against the actual target-journal rules before moving it. Current Main.tex has eleven figure floats and four table floats; no journal limit was presumed by this audit.

## P2 — Table fit, χ₄, and notation

The paper's fourteen-column `tab:FlxPinGof` leaves the page and has the wrong orientation. Reduce unnecessary repetition or split into readable tables before shrinking text. Check the final page in its intended reading orientation. Preserve values and memberships from the numerical audit.

Remove the paper's empty χ₄ column in `tab:AdoptedXi`; this model declares χ₀–χ₃. Explain why χ₃ is not fitted to pinned-flexor identification but is used in later flexor redesign, rather than reconciling the paper dash and dissertation 0.158 silently.

Ben's dissertation requests also apply to shared equations: uppercase torque M, and `\mathbf{k}` for knee-axis subscripts, distinguishing them from basis/summation indices. Apply consistently to axes, captions, equations, and generated figures. Preserve moment-arm variables.

## P2 — Dissertation mechanics schematic

At `fig:xiProjection`, split the existing compliance/wrapping composite into two figures as Ben requested and put alpha on the circle. The follow-up “Then we don't need to mention” has an incomplete instruction; inspect its saved source anchor before deleting an explanation. This request is dissertation-specific; the paper has a different mechanics-figure set.

## Return package

Use the established `Dissertation/Notes/advisor_revision_20261023/REVIEW_INBOX/<machine>/<new-batch>/` with `REVIEW.md`, bounded `patches/`, and `figures/`. Include original and replacement filenames, labels, numerical changes, and a series-level CSV for audited plots. Keep the current GoF prose edits from this audit. Use the dissertation, MATLAB, mechanics, and relevant plotting skills. No new optimization is needed for these repairs.
