# Chapter 4 and SADb directives 29 October 2026

Local changes only. No Overleaf upload, comments resolved, tracked changes accepted, commits, or pushes. The annotated Word document has not been found: the supplied local path was checked after Ben closed Word and contains no comments part, tracked insertions/deletions, or new text. Office recent-file shortcuts point to that same path. No connected Word session was available. The original review document is preserved and not overwritten; annotation-dependent revisions remain pending.

## Directive compliance

| Directive | Evidence |
|---|---|
| Section 4.3.4 title | `30-torque.tex:727` |
| Biomimetic robot framing | `30-torque.tex:665` |
| Direct-route summary table | `30-torque.tex:667` |
| Cam geometry and capacity table | `30-torque.tex:687` |
| Schematic caption | `30-torque.tex:781` |
| Origin reaction | `30-torque.tex:774` |
| Views callout | `50-database.tex:8` |
| Map controls callout | `50-database.tex:26` |
| Citation connections callout | `50-database.tex:28` |

- Figure 4.11: removed source tile 4 (the redesigned-route panel C); retained source tiles 1 and 3 as a one-by-two layout. The actual XData and YData of retained objects are bit-identical. Legend is southwest within panel B. The removed C results were previously duplicated in the later physical-test section; the earlier discussion is removed and the later values retained under Exceeding Human Torque Magnitude.
- Figure 4.14 (flexor physical test, not extensor): all scatter objects have MarkerFaceAlpha 0.75 and MarkerEdgeColor none. The group carrying the real measured series is stacked top/last, above pressure-specific predictions. XData and YData are unchanged. All edits operate on saved figure objects; no model evaluation or optimization was run.
- Section 4.3.3: describes inherited biomimetic robot geometry with human-muscle torque references, rather than a redesign solely of Gait2392. Prose source words decrease from 356 to 125, excluding tables and figure blocks. Two summary tables retain the direct extensor/flexor predictions and the cam capacity, pressure, and geometric values. Redundant peak-value repetition and the claim that +5.02 percent is exactly 5 percent are removed; the underlying numbers are retained. No fit coefficients, numerical results, or pressure data were changed.
- Reverse pulley: caption starts Example of. The schematic labels F_BPA=2F_20, F_t=F_BPA/n, transmitted tendon travel n times BPA contraction, and origin reaction F_BPA. n counts cable parts supporting the moving block, not arbitrary redirect-pulley wheels. With the perpendicular exit shown, the fixed-mount magnitude is F_BPA sqrt(1+1/n^2), not F_BPA(n-1). The full vector equation and angle-dependent magnitude remain in the text. For n=2 this is sqrt(5)F_20. The origin and common fixed-tackle reactions are distinct; their combination on the entire femur is not labeled as the fixed-mount reaction. The old reverse_pulley_principle.pdf was locked against replacement, so the active local include is reverse_pulley_biomimetic.pdf; the prior graphic is preserved.
- SADb: A through F are callouts, not subpanels. Interface areas are called the record area and connection list. Explicit Figure reference plus callout letter is used for each discussion of a marked feature. The six-letter caption remains intact.
- FVU remains addressed by Ben's reply, with no new interpretation or recalculation.

## Formatting sources and bounded finding

PSU Graduate School, Electronic Thesis and Dissertation Formatting Requirements, https://www.pdx.edu/gradschool/etd-formatting-requirements (checked 29 October 2026): figures must be titled/labeled, sequentially numbered, legible, and consistently formatted. PSU permits disciplinary style provided its own requirements take precedence. This page does not prescribe the exact wording for lettered screenshot callouts.

PSU ETD and Format Review Resources, https://www.pdx.edu/gradschool/etd-and-format-review-resources: PSU does not require a separate style guide beyond its formatting requirements. Therefore “Figure ..., callout (B)” is a consistent editorial choice, not a quoted PSU mandate. Repeating the figure reference for each callout avoids mistaking the marked UI features for separate figure panels. This was a focused guidance check, not a complete submission-format audit.

## Verification

Saved MATLAB object checks are in Figures/30-torque/Components/directives_20261029/plot_verification.txt. Editable revised .fig files and the editing script are beside it. Original chapter sources and graphics are in Archive/ch45_directives_20261029. Published Chapter 20 and the separate torque paper were not edited. Final build and visual-check results are appended after completion.

Final verification: build_chapter45_directives_release_20261029, 306 pages, five successful passes, no errors, unresolved references, duplicate labels, or missing figures. One inherited 0.66867 pt overfull box remains. Visual checks covered the 1x2 flexor figure, both new tables, renamed subsection, measured-last flexor plot, schematic, and SADb callout references. Figure 4.11 fonts were increased to remain legible at final inclusion size. Build folders and QA images are archived in Archive/ch45_directives_20261029. The annotated Word-copy question remains pending; no Word document was overwritten in this turn.


## Saved annotation follow-up

The 16:14 saved Word file contains Ben’s plain-text instructions and anchors. The earlier pending question is superseded. See CH45_ANNOTATED_REVIEW_20261029.md and CH45_COMMENTS_REPLIES_AND_FIXES_REVISED_20261029.docx for the applied directives, newly read Overleaf reply, tracked-deletion reconciliation, build, and remaining items. The annotated input was preserved unchanged.
