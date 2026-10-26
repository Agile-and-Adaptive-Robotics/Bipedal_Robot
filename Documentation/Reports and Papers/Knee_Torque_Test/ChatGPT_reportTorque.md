# Torque-paper and dissertation reconciliation

Status: audit and requested GoF prose revision in progress, 26 October 2026.

Paper: [Design Factors for Isometric Robotic Knee Torque](https://www.overleaf.com/project/622e9536b8ea1e0a7f167cf5). Dissertation: [Bolen-Dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb), `30-torque.tex`.

Live source snapshots and review evidence are preserved in [review_audit_20261026](review_audit_20261026/). This report will distinguish Ben's pre-existing comments and tracked edits from changes made during this audit.

## Findings established so far

- The paper review overview contains 31 active items: 29 in `Main.tex`, one bibliography insertion, and one `latexmkrc` insertion. Nine active comment threads are Ben's; two are older advisor comments. The full items are preserved in `paper_review_entries_with_positions.json`.
- The paper has 272 resolved historical threads. All expandable text was opened, and the full archive, including 74 substantive messages by Ben, is preserved in `paper_resolved_comments.json`. Resolution events are retained separately from substantive comments/replies.
- The dissertation torque review overview contains 13 active items. These include Ben's uppercase torque-symbol request, bold knee-axis subscript request, and instruction to separate a composite figure and mark alpha on the circle.
- The paper explains RMSE, FVU, and maximum absolute residual in its Identification Procedure, but less fully than the dissertation. The dissertation's claim that the metrics separately map to offset, shape/amplitude, and phase is too strong. Ben has authorized an improved explanation in both documents that makes the role of multiobjective optimization clear.
- Ben confirmed during this audit that the 2×2 layout request applies to `fig:FlxGrp`: three data panels titled **10 mm**, **20 mm**, and **20 mm**, with the legend in the fourth panel. Its placement as a reply to the χ₄-table thread does not change the intended target.

## Current replot and correction targets

1. `fig:FlxGrp`: confirmed 2×2 layout, three data panels and a dedicated legend tile. Preserve distinctions between the original 10 mm/20 mm configurations and the redesigned two-BPA 20 mm route.
2. `fig:ExtensorRoute_redesigned`: reduce the oversized multi-pose geometry figure; add femur/tibia outlines so the routing is intelligible. Ben permits fewer pose panels.
3. `tab:FlxPinGof`: correct page fit and rotation, retaining data and training/validation labels.
4. `fig:FlxPin` / Hunt equation: audit equation, units, strain domain, valid-point mask, and scores before regenerating; Ben questioned the calculation.
5. Incorrect measured/hybrid point series: Ben's comment navigates to the `fig:Ext20mm` figure environment. Verify the intended data series and figure generator before changing marker colors or values.
6. `tab:AdoptedXi`: the paper has an extra χ₄ column although the declared correction family is χ₀–χ₃. The 2×2 figure request is a separate action on `fig:FlxGrp`.
7. Initial identification narratives: Ben marked the flexor paragraph's values and training/validation membership incorrect and marked the initial extensor paragraph “all incorrect.” Both have corresponding historical-fit material in the dissertation that requires source reconciliation.

No plots have been regenerated and no comments or tracked changes have been resolved, accepted, or rejected during this audit. Detailed per-item anchors, the document differences, and the final GoF change record follow as the audit completes.
