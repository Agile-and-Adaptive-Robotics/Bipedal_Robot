# ChatGPT torque-paper / dissertation report

Completed 2026-10-26T15:09:07-07:00. Machine: EB475WS4.

**Main findings:** the two documents contain materially different figures, data counts, historical-fit descriptions, and presentation of results. Several of Ben's current paper comments also apply to the dissertation. A blind chapter replacement or automatic transfer from either document would propagate errors.

Paper: [Design Factors for Isometric Robotic Knee Torque](https://www.overleaf.com/project/622e9536b8ea1e0a7f167cf5). Dissertation: [Bolen-Dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb), `30-torque.tex`.

## Scope and evidence

- Captured live source from both Overleaf editors. Line numbers below refer to the saved **pre-GoF-edit snapshots**, not PDF line numbers or character positions. Use section/figure labels to locate the current live passages.
- Paper review overview: **31 active items** before this audit's GoF edit: 29 in Main.tex, one in BolenFrontiers22.bib, and one in latexmkrc. Main.tex has **11 comment threads** (nine Ben, two older advisor) and **18 tracked changes** (17 Ben, one advisor). The two other files each have a tracked insertion by Ben.
- Dissertation torque overview: **13 active items** before this audit's GoF edit: four comment threads and nine tracked changes, all attributed to Ben.
- Resolved paper history: **272 threads**, fully expanded and preserved, including **74 substantive messages by Ben**. Resolution events are not counted as substantive replies.
- Only the subsequently authorized GoF prose was edited during this audit. Figures, tables, numerical results, pending review dispositions, and comments were not changed. No commits, pushes, model runs, or messages to ZCode.

## GoF: does the paper explain it?

**Yes, but briefly.** The paper's Identification Procedure (snapshot line 300) defined RMSE, FVU, and maximum absolute residual, and named the solver. The dissertation had a longer account at lines 310–314. Its one-to-one mapping of metrics onto offset, shape/amplitude, and phase was too strong.

**Ben's clarification is applied in both documents:** the metrics are interpreted together, rather than assigning one number to each type of curve disagreement. The revised explanation states what each metric computes, how their combined interpretation assesses torque magnitude and variation, and why multiobjective optimization retains trade-offs between overall agreement, comparison across tests, and local error. The Pareto front and held-out selection are explained. Existing algorithm settings and selection procedures were retained.

For a fixed test, RMSE and FVU are transformations of the same squared residuals; their aggregation across tests of different scales can weight those tests differently. FVU does not exclude mean bias, and maximum residual is not insensitive to a uniform offset. This mathematical point is why the combined interpretation is appropriate.

Exact applied prose: [gof_prose_patches.json](review_audit_20261026/gof_prose_patches.json). Whole-file readbacks matched precisely one bounded replacement in each document. [Applied diff](review_audit_20261026/applied_GoF_prose.diff).

## Current paper comments and corresponding dissertation actions

| ID / Ben's comment | Paper anchor | Dissertation finding / required action |
|---|---|---|
| P-C01 | This figure got replaced with an incorrect one. The blue dots are where the back-calculated yellow dots should be. | Figure environment `fig:Ext20mm`, lines 506–510; `Figures/Ext20mm_cam_20261023.eps`. | Compare with `fig:extensor_physical_final`, lines 599–604, `extensor_20mm_physical.png`. The paper plots six blue measurements; the dissertation describes nine measurements plus yellow pressure-specific predictions. Audit data provenance before changing colors or transferring either result. |
| P-C02 | Don't use "holds". Call it data measurements. | Physical extensor results, especially line 504 and caption line 510. | Still applies: `sec:redesign_physical`, line 595, says nine/fifteen “holds,” and later text says “between holds.” Search both documents and captions for measurement uses; replace with data measurements. Do not alter the mechanical meaning of holding a load when describing the procedure. |
| P-C03 | Too big for page. Should probably have a femur and tibia shape in there to make it less confusing. We could reduce the amount of panels shown. | `fig:ExtensorRoute_redesigned`, lines 485–490; eight poses plus legend. | Dissertation uses a different composite under the same label, lines 551–567: cropped extensor_cam_route.pdf panels and a hardware photograph. Review leg outlines and final-page readability in each version; the paper's eight-pose replot should not silently replace the dissertation composite. |
| P-C04 | Candidate to move to supplementary material if we are over the Figure/Table limit. | Blank-line anchor directly before `fig:RouteVariants`, lines 474–480; confirmed by comment navigation and surrounding source. | No equivalent RouteVariants figure in the torque chapter. Its route/pulley material is distributed between the chapter and Appendix C. Treat this as a conditional paper-length decision, not an instruction to delete dissertation results. Main paper currently has 11 figure floats and four table floats; journal limits were not assumed. |
| P-C05 | There's no /Chi_{4} lol XD. Reply: Make a 2x2 panel figure with the legend in the last panel. Title panels 10 mm, 20 mm, 20 mm. | `tab:AdoptedXi` header, line 457, has a χ₄ column. | Dissertation table at lines 528–540 already has χ₀–χ₃ only. Remove the empty paper χ₄ column in a future table edit. **Ben explicitly confirmed that the reply's 2×2 instruction applies to fig:FlxGrp**, not this table. Preserve these as two separate actions. |
| P-C06 | Leaves the page and is rotated the wrong way. | `tab:FlxPinGof`, especially line 356; 14 columns in a sidewaystable. | Dissertation separates pre-optimization and adopted-term tables and already carries Ben's tracked table-layout changes (lines 350–407). Fit and orientation must be checked in the paper's article layout separately. Do not rotate until the final reader-facing page is inspected. |
| P-C07 | I don't think Hunt Equation was calculated correctly here. | `fig:FlxPin`, caption line 334; `FlxPin_allmodels_20261023.eps`. | Dissertation pre-optimization Figure `fig:FlxPin_preopt` and its GoF table, lines 341–368, use the same Hunt comparison. Audit the calculation and score mask in both documents. Existing low-strain undefined-point exclusions do not establish that the equation itself was evaluated correctly. |
| P-C08 | This is all incorrect | Initial extensor paragraph, line 368. | Corresponding historical identification paragraph remains in the dissertation at line 425. Reconcile with saved fit records, intended parameter stage, held-out fold, and the terms actually used to generate the displayed curves. Do not relabel old single-bracket numbers as adopted two-stage results. |
| P-C09 | This is incorrect in the values and which tests were training and validation. The whole paragraph is incorrect. | Initial flexor paragraph, line 328. | Near-verbatim corresponding paragraph remains in the dissertation at line 375, including 11.9 mm, 102/10.6 kN/m, initial-fit errors, and a single validation test. Both need the same record-level reconciliation. Adopted two-bracket fold is two training/two validation; identification history must be distinguished from current figure membership. |
| P-C10 | Advisor, 6 May: how does the shorter length result in this? | Old thread currently navigates to the conclusion near line 556. Original highlighted phrase no longer survives. | Anchor is historically uncertain; retain the thread. The force/resting-length and wrapping explanations may be relevant, but do not claim the original question answered without the original quoted context. |
| P-C11 | Advisor, 6 May: Can this be connected with Jon's work? | Old thread also navigates to conclusion line 556; original anchor no longer survives. | Identity of Jon's specific work is not established by this thread. Preserve as unresolved historical context; do not invent a reference or retrofit a claim. |

## Differences beyond the comment anchors

| Issue | Paper | Dissertation / disposition |
|---|---|---|
| Cam-extensor data count and result | Six measurements; only one exceeds the human target at the measured, lower pressures. | Nine measurements and a broader “most of the tested range” summary. Current `Knee_Extensor_20mm.m` contains nine entries. This needs a shared data export and recalculated comparisons, not a prose-only merge. |
| Biomimetic extensor GoF | Line 407 reports adopted RMSE 1.1702, FVU 0.6681, max 2.2088; line 468 also distinguishes older terms. | Line 465 reports 0.8920 / 0.3396 / 1.4188, while line 540 separately reports about 1.17 / 0.67. The dissertation is internally inconsistent about the parameter/evaluation version. |
| Flexor panel C | `FlxGrp_20261023.eps` is a single merged three-panel asset; caption fully describes the redesigned route. | A/B in FlxGrp.eps plus a small separate C from Flx20mm_meetsHuman.pdf. C caption ends at “The redesigned two-BPA flexor route.” Carry over a complete explanation with the confirmed 2×2 replot. |
| Original 10 mm extensor photograph/curve | One A/B figure, photo A and torque plot B. | Separate torque plot and pass-through photograph with two labels. A future merge must retain reference aliases and update panel references. |
| Adopted χ₃ scope | Flexor row is a dash, extensor 0.158; table describes identification. | Both rows show 0.158, caption says it was identified on extensor and used in flexor redesign, not fitted to pinned flexor. Explain identification versus deployment instead of silently harmonizing cells. |
| Pulley claims | Route procedure says mono is one tackle and bi one tackle per joint; says all record runs were unity transmission. | Chapter describes a later size-resolved transmission study; Appendix C has pulley torque/margin surfaces. Distinguish the older unity-routing runs from the later study. Earlier reviewed rulings describe mono/bi by the dominant DoF used in the route calculation; verify implementation before porting the paper explanation. |
| Deleted unity-gain paragraph | Ben deleted the long “The knee-flexor run quantifies…” paragraph at 2:20 pm. | The first paragraph of `sec:disc_pulley`, line 589, still leads with unity transmission and 81.23→85.74 N m / 0.803 N m. Mark for editorial reconciliation with Ben's deletion; do not restore the long paragraph from history. |
| Design-factor wording | Ben changed budget→allowance, added the biomimicry trade-off, and added isokinetic/isotonic future-testing discussion. | Chapter contains related anatomical-bounds and routing trade-offs, but not the same explicit design-factor synthesis. Selectively carry over meaning, correcting “biomimcry” and the repeated “in in”; do not copy tracked prose mechanically. |
| Removed caveat block | Ben deleted a long limitations paragraph, then scoped the remaining statement to this stand/configurations/actuator families. | Related limits are distributed in route and physical-test sections. Preserve the deletion as editorial history; retain necessary, accurately scoped limitations without importing the deleted block wholesale. |
| Whole-body capability survey | Section `sec:capability_survey` reports a 19-actuator 5×5-domain comparison (66.32→75.30% torque coverage; 97.05→95.23% active-force coverage). | Not a corresponding standalone torque-chapter section. Appendix C has separate route-search/capability-support material. Verify whether the same campaign is represented before cross-referencing; do not conflate these coverage statistics with route feasibility counts. |
| Torque/knee-axis notation | Existing paper equations retain several lowercase moment/torque forms and plain knee subscripts. | Ben explicitly requests uppercase torque M and bold knee-axis subscript k in dissertation comments. Apply consistently to both manuscripts in a separate notation pass, distinguishing torque from moment-arm variables and knee subscripts from summation indices. |
| Composite mechanics illustration | Paper has the frame figures but no `fig:xiProjection` composite. | Ben asks to split `fig:xiProjection` and put alpha on the circle. This is a dissertation-specific replot, not a paper change to force through without checking its different figure set. |
| GoF exposition | Shorter before this audit. | Longer but overassigned metrics to specific error causes. **Both revised live during this audit** to describe their combined effect and multiobjective selection. |
| Paper conclusion grammar | Ben's tracked opening currently reads “started with an improved … model and actually delivers about a joint.” | Preserve the recorded edit, but polish the resulting sentence before submission; it is not a suitable sentence to transplant into the dissertation. |
| Compile warnings | GoF revision compiles with 0 errors and 20 warnings: duplicate fig:extensor_misalign, small-font substitutions, and their summary. | Record/reconcile these remaining warnings separately from scientific revisions. No new labels or citation keys were introduced by the GoF edits. |

## Every active paper item: exact text and source location

The original UI character positions are preserved in JSON as anchor evidence. They are converted to snapshot line numbers here; they are not presented as line numbers themselves. The excerpt is the containing source line, not a claim to recover the complete highlighted range. Short tracked fragments are grouped by their actual location, so “is→was” is not treated as a free-standing scientific instruction.

### P-I01 — Insertion — keyboard_arrow_downBolenFrontiers22.bib1

```text
bbolen83
14 April, 3:40 pm
Added:Â 
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I02 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:32 pm
This is incorrect in the values and which tests were training and validation. The whole paragraph is incorrect.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I03 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:30 pm
I don't think Hunt Equation was calculated correctly here. 
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I04 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:29 pm
Leaves the page and is rotated the wrong way.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I05 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:31 pm
This is all incorrect
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I06 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:28 pm
There's no /Chi_{4} 
lol XD
bbolen83
26 October, 2:44 pm
Make a 2x2 panel figure with the legend in the last panel. Title panels 10 mm, 20 mm, 20 mm.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I07 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:27 pm
Candidate to move to supplementary material if we are over the Figure/Table limit.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I08 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:26 pm
Too big for page. Should probably have a femur and tibia shape in there to make it less confusing. We could reduce the amount of panels shown. 
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I09 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:21 pm
Changed: is
to was
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I10 — Deletion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:20 pm
Deleted:Â  The knee-flexor run quantifies what the mono-pulley freedom is worth at unity gain. As a uniarticular knee flexor, the biceps femoris short head is exactly the site the mono-pulley variant serves, and the run converged to a collision-feasible layout with zero bound violations, both actuators of the rig reeved over a single tackle of one line part, so the gain was 1 and no force multiplication was engaged; the run therefore measures the redirect baseline, not the transmission trade. Even the redirect alone deepened the peak predicted flexion torque from \qty{81.23}{\newton\meter} on the direct route to \qty{85.74}{\newton\meter} on the same knee-axis measure, \qty{85.56}{\newton\meter} on the insertion-side measure after the tackle's force loss, a gain of \qty{4.51}{\newton\meter} or \qty{5.6}{\percent}, while holding the off-axis torque to at most \qty{0.803}{\newton\meter}, under \qty{1}{\percent} of the peak. The direct and pulley runs optimize different objective functions, so their objective values, 74.47 and 0.08465, are not commensurable; the torque comparison, not the objective comparison, is the meaningful one. That comparison should still be read as the value of the added routing freedom in this pair of runs rather than as an isolated pulley effect, since no run at fixed objective separates the pulley's contribution. The predicted peaks of the direct and pulley routes exceed the human flexor scale of \qty{23.16}{\newton\meter} by factors of 3.51 and 3.70, and the torque-for-range exchange that a reeved tackle buys is quantified by the doctrine re-solves of the batch analysis that follows.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I11 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:18 pm
Don't use "holds". Call it data measurements. 
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I12 — Comment — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:18 pm
This figure got replaced with an incorrect one. The blue dots are where the back-calculated yellow dots should be.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I13 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:12 pm
Changed: The route should keep $\varepsilon^{*}$ high through the functional range, because the wrapping loss of Equation~(\ref{eq:deltaL}) scales with $(1 - \varepsilon^{*})^{2}$ and therefore grows exactly where an under-contracted actuator is weakest.
to For a given $l_{0}$,
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I14 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:12 pm
Changed: M
to m
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I15 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:13 pm
Changed: budget
to allowance
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I16 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:13 pm
Added:Â  This increased moment arm comes at a cost to biomimcry.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I17 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:15 pm
Changed: filament
to filament
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I18 — Deletion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 1:55 pm
Deleted:Â The redesign results carry caveats that bound their use. The flexor redesign solves its second actuator's endpoints from the first by the mirrored-endpoint rule rather than optimizing the two routes independently, so the mirror symmetry is built into the design; freeing the second route's endpoints is the next design refinement. All redesign solutions are optimizer-stage: they have not been propagated into the simulation models, and the physical validation rests on sparse isometric holds, six points on the extensor cam route and fifteen on the flexor route, rather than dense sweeps of either configuration. The stiffness pair itself was selected by an external design requirement rather than uniquely identified, because the pinned-flexor data admit a family of near-equivalent pairs; an independent stiffness measurement of the test stands would close that non-identifiability and either confirm or replace the adopted pair. The wrapping factor $\chi_{3}$ was identified on the $\diameter$\qty{10}{\mm} extensor family and applied unchanged to the $\diameter$\qty{20}{\mm} redesigns, so its transfer across actuator diameters is untested. On the characterization side, the fictive-strain formulation inherits the conclusion of our force characterization that the unloaded \qty{620}{\kPa} contraction $\varepsilon_{620}$ should be recorded for each actuator used on a robot \citep{bolen_2026}, since every correction term is normalized by it. The whole-body extension identified its own next steps: the route-builder assembly failure needs a fix before the twelve terminated routes can be rebuilt, and the two torque-short feasible routes need either larger actuators or relaxed corridors.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I19 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 1:56 pm
Changed: One scope statement, finally, belongs with these caveats: t
to T
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I20 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 1:56 pm
Changed: ese
to is
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I21 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 1:56 pm
Changed: s
to , these knee configurations,
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I22 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 1:57 pm
Changed: is
to ese
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I23 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 1:57 pm
Changed: y, and o
to ies. O
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I24 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:02 pm
Changed: sizes
to types
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I25 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:02 pm
Added:Â  Isometric torque testing has yielded important design factors for roboticists seeking to design biomimetic humanoid robots. Additional testing such as isokinetic and isotonic that are frequently used in in human dynamometer testing may yield equally as important insights.
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I26 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:09 pm
Changed: asked what a
to started with an improved
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I27 — Insertion — keyboard_arrow_downMain.tex29

```text
bbolen83
26 October, 2:09 pm
Added:Â force-length-pressure model and
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I28 — Deletion — keyboard_arrow_downMain.tex29

```text
ajh26
6 May, 12:01 pm
Deleted:Â s/a theX_0K_1K_2X_{0}KX_{0}
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I29 — Comment — keyboard_arrow_downMain.tex29

```text
ajh26
6 May, 12:06 pm
how does the shorter length result in this?
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I30 — Comment — keyboard_arrow_downMain.tex29

```text
ajh26
6 May, 12:08 pm
Can this be connected with Jon's work?
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

### P-I31 — Insertion — keyboard_arrow_downlatexmkrc1

```text
bbolen83
7 November, 2:38 pm
Added:Â # This shows how to use the xr package with latexmk. # John Collins 2023-03-29 # # The xr package ("a system for eXternal References") is used by a document # to make references to sections, equations, etc in other external # documents. # The definitions in this file enable latexmk to apply latexmk to # automatically update an external document whenever its .tex file changes, # so that the references in the main document stay up to date. # Notes: # 1. This version is defined to put the files from the compilations of # the external documents into a defined subdirectory, to segregate # potentially many generated files from the main document # directories. # 2. But for latexmk's custom dependency mechanism to be used, as here, # the aux file from compilation of a subdocument must be generated in # the same directory as the corresponding source .tex file. So the # .aux file is copied. # 3. It is assumed that the external documents are to be compiled by # pdflatex. This can be changed, of course, by changing the '-pdf' # option given to the invoked latexmk to whatever is needed. # 4. An ideal implementation would also ensure that recompilation of an # external document also happens whenever any of its other source # files changes. But this is not done in the present version, and # would probably entail either the use of internal latexmk variables # or extra enhancements to latexmk. # 5. The code uses subroutines copy and fileparse that are loaded by # latexmk from the Perl packages File::Copy and File::Basename. # 6. It also uses some not-yet-documented features of latexmk: an array # variable @file_not_found and subroutines popd, pushd, and # rdb_add_generated. #-------------------- # Configurable choices for compilation of external documents # Subdirectory for output files from compilation of external documents: $sub_doc_output = 'output-subdoc'; # Options to supply to latexmk for compilation of external documents: @sub_doc_options = (); push @sub_doc_options, '-pdf'; # Use pdflatex for compilation of external documents. # Replace '-pdf' by '-pdfdvi', 'pdfxe', or 'pdflua' if needed. #-------------------- # Add a pattern for xr's log-file message about missing files to latexmk's # list. Latexmk's variable @file_not_found is not yet documented. # This line isn't necessary for v. 4.80 or later of latexmk. push @file_not_found, '^No file\\s*(.+)\s*$'; add_cus_dep( 'tex', 'aux', 0, 'makeexternaldocument' ); sub makeexternaldocument { if ( $root_filename ne $_[0] ) { my ($base_name, $path) = fileparse( $_[0] ); pushd $path; my $return = system "latexmk", @sub_doc_options, "-aux-directory=$sub_doc_output", "-output-directory=$sub_doc_output", $base_name; if ( ($sub_doc_output ne '') && ($sub_doc_output ne '.') ) { # In this case, .aux file generated by pdflatex isn't in same # directory as the .tex file. # Therefore: # 1. Actual generated aux file must be listed as produced by this # rule, so that latexmk deals with dependencies correctly. # (Problem to overcome: If $sub_dir_output is same as $aux_dir # for the main document, xr may read the .aux file in the # aux_dir rather than the one the cus dep is assumed by latexmk # to produce, which is in the same directory as the .tex source # file for this custom dependency.) # Use not-yet-documented latexmk subroutine rdb_add_generated # to do this: # 2. A copy of the .aux file must be in same directory as .tex file # to satisfy latexmk's definition of a custom dependency. rdb_add_generated( "$sub_doc_output/$base_name.aux" ); copy "$sub_doc_output/$base_name.aux", "."; } popd; return $return; } }
```

File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.

## Every active dissertation torque item

### D-I01 — Insertion — line 25

```text
bbolen83
26 October, 1:16 pm
Added:Â International Society of Biomechanics (ISB)
```

Containing section: **Robot Leg Model and Knee Test Stands**.

> Details on the biomimetic robot knee used in this study can be found in \citep{steele_biomimetic_2018,steele_development_2017,steele_biomimetic_2018-1}. In brief, this 1-DoF joint (Fig. \ref{fig:ICR}) uses a four-bar linkage to change the effective moment arm during rotation. The intersection point of the internal crossed links define the instantaneous center of rotation (ICR). This causes the joint ICR location in the femur frame to translate posteriorly in the $X$ direction and distally in the $Y$ direction except for a brief proximal move near full flexion. Reference International Society of Biomechanics (ISB) coordinate frames are given by \citet{wu_for_1995,wu_isb_2002}. This is similar to how the human knee behaves \citep{bolen_determination_2019,morrow_optimization_2020,seth_opensim_2018,delp_opensim_2007,delp_interactive_1990}. The correct location of the ICR is important for both theoretical and experimental torque calculations.

### D-I02 — Insertion — line 34

```text
bbolen83
26 October, 1:17 pm
Changed: U
to In the pinned knee joint, u
```

Containing section: **Robot Leg Model and Knee Test Stands**.

> In the pinned knee joint, uniarticular knee extensor and knee flexor are used to actuate the joint. Each knee joint type was tested with a single $\diameter$\qty{10}{\mm} Festo BPA per muscle. The biomimetic knee was additionally tested with $\diameter$\qty{20}{\mm} Festo BPAs in a paired configuration: the extensor with two BPAs in parallel, and the flexor with two BPAs whose pathways are mirrored about the $Z$ axis. Each artificial muscle is pinned at the muscle origin location on the femur while the other end is either pinned or attached to the muscle insertion location on the tibia via an artificial tendon.

### D-I03 — Comment — line 36

```text
bbolen83
26 October, 1:20 pm
I don't like Torque as lower case m. Use upper case. 
```

Containing section: **Robot Leg Model and Knee Test Stands**.

> Torque about the knee is generated by the moment of the BPA force about the joint's instantaneous center of rotation (ICR). Two vectors enter that moment: the position vector $\mathbf{d}$, drawn from the ICR to the point where the actuator force is applied, and the force vector $\mathbf{F} = \norm{\mathbf{F}}\hat{\mathbf{u}}$, whose direction $\hat{\mathbf{u}}$ is the line of action of the actuator. The moment of the force about the joint is $\mathbf{m} = \mathbf{d} \times \mathbf{F}$, and its magnitude is the force magnitude multiplied by the moment arm $r$, the perpendicular distance from the joint axis to the line of action,

### D-I04 — Comment — line 43

```text
bbolen83
26 October, 1:24 pm
Make all subscript k's that refer to knee axis styled as \mathbf{k} to avoid confusion with i,j,k aka x,y,z, aka e_1, e_2, e_3 basis vectors.
```

Containing section: **Robot Leg Model and Knee Test Stands**.

> The moment arm is therefore a scalar, and it is the same quantity whether it is computed from the position vector with Eq.~\eqref{eq:momentarm_def} or measured directly with a tape measure. In this study it is calculated using the method specified by \citet{young_analyzing_2019}. The scalar moment arm $r_{k}$ about the joint axis $\hat{\mathbf{z}}_{k}$ (i.e., the $Z$ axis in the knee frame) is calculated as 

### D-I05 — Comment — line 270

```text
bbolen83
26 October, 1:38 pm
Should be two separate figures. alpha can be marked on the circle. 
```

Containing section: **Length Offset and System Stiffness**.

> \includegraphics[width=0.83\textwidth,height=0.20\textheight,keepaspectratio]{bpa_wrapping_loss.pdf}

### D-I06 — Comment — line 308

```text
bbolen83
26 October, 2:35 pm
Then we don't need to mention
```

Containing section: **Multiobjective Identification Algorithm**.

> This initial identification was subsequently extended to the full test sets in two stages. First, the flexor model was given a second compliant bracket at the origin-side attachment (a hip-frame cantilever carrying the BPA origin), so that the compliance of both ends of the actuator path is represented; the stiffness pair $(\chi_{1}, \chi_{2})$ is shared by the two brackets in fixed axis assignments (Appendix~\ref{app:model_geometry}, Table~\ref{tab:app_stiffness_arrays}), leaving a three-term fit $(\chi_{0}, \chi_{1}, \chi_{2})$ cross-validated over the four remaining flexor tests, two of them held out per fold; the encoder-corrected $l_{0}=\qty{47.9}{\cm}$ test was not used in the two-bracket identification. Second, the extensor refit locked $(\chi_{1}, \chi_{2})$ at the flexor-identified values, enforcing a single stiffness pair across the flexor and extensor configurations, and solved only the extensor-specific terms $(\chi_{0}, \chi_{3})$ on a three-test training set with three tests held out for validation, with the extensor length offset allowed to take negative values. The bounds and the adopted values of both stages are recorded in Appendix~\ref{app:model_geometry} (Tables~\ref{tab:app_opt_bounds} and~\ref{tab:app_xi_values}); the fit outcomes are summarized in Tables~\ref{tab:FlxPinGof} and~\ref{tab:ExtPinGof}.

### D-I07 — Insertion — line 352

```text
bbolen83
23 October, 10:56 am
Changed: \begin{sidewaystable}[p] \centering\small\setlength{\tabcolsep}{4pt} \caption[Pre-optimization pinned-flexor goodness of fit]{Pre-optimization pinned-flexor goodness of fit. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. Scoring-point restrictions are stated in the text.} \label{tab:FlxPinPreGof} \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}lccccccccc@{}} \toprule & \multicolumn{3}{c}{\textbf{Rigid body}} & \multicolumn{3}{c}{\textbf{Improved char.\ alone}} & \multicolumn{3}{c}{\textbf{Back-calculated (hybrid)}} \\ \cmidrule(lr){2-4} \cmidrule(lr){5-7} \cmidrule(l){8-10} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{48.5}{\cm}$ (training) & 5.779 & 1.1761 & 9.798 & 4.871 & 0.8354 & 9.865 & 1.635 & 0.0991 & 3.842 \\ $l_{0}=\qty{41.0}{\cm}$ (training) & 6.633 & 3.2313 & 11.409 & 6.074 & 2.7090 & 9.262 & 2.187 & 0.1184 & 4.243 \\ $l_{0}=\qty{45.7}{\cm}$ (validation) & 5.315 & 0.9408 & 9.991 & 4.465 & 0.6640 & 9.536 & 1.246 & 0.0409 & 2.275 \\ $l_{0}=\qty{40.6}{\cm}$, tendon (validation) & 3.945 & 0.6750 & 7.527 & 3.226 & 0.4512 & 7.376 & 2.085 & 0.0755 & 3.073 \\ \midrule Mean (training) & 6.206 & 2.2037 & 10.604 & 5.472 & 1.7722 & 9.563 & 1.911 & 0.1088 & 4.043 \\ Mean (validation) & 4.630 & 0.8079 & 8.759 & 3.845 & 0.5576 & 8.456 & 1.666 & 0.0582 & 2.674 \\ \bottomrule \end{tabular} \end{sidewaystable}
to \begin{table}[H] \caption{Pre-optimization pinned-knee flexor torque-prediction goodness of fit per test, comparing the three pre-optimization torque predictions of Figure~\ref{fig:FlxPin_preopt}: the Hunt-equation rigid-body prediction, the unoptimized prediction of the improved BPA characterization, and the back-calculated (hybrid) torque. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:FlxPinPreGof} \begingroup \renewcommand{\arraystretch}{1}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc ccccccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Hunt equation}} & \multicolumn{3}{c}{\textbf{Improved char.\ alone}} & \multicolumn{3}{c}{\textbf{Back-calculated (hybrid)}} \\ \cmidrule(lr){3-5} \cmidrule(lr){6-8} \cmidrule(l){9-11} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule 48.5 & -- & 5.779 & 1.1761 & 9.798 & 4.871 & 0.8354 & 9.865 & 1.635 & 0.0991 & 3.842 \\ 41.0 & -- & 6.633 & 3.2313 & 11.409 & 6.074 & 2.7090 & 9.262 & 2.187 & 0.1184 & 4.243 \\ 45.7 & -- & 5.315 & 0.9408 & 9.991 & 4.465 & 0.6640 & 9.536 & 1.246 & 0.0409 & 2.275 \\ 40.6 & 35 & 3.945 & 0.6750 & 7.527 & 3.226 & 0.4512 & 7.376 & 2.085 & 0.0755 & 3.073 \\ \bottomrule \end{tabular} \endgroup \end{table}
```

Containing section: **Flexor on Pin Joint**.

> \begin{table}[H]

### D-I08 — Insertion — line 384

```text
bbolen83
23 October, 10:59 am
Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-flexor goodness of fit]{Optimized pinned-flexor goodness of fit at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.} \label{tab:FlxPinGof} \begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}} \toprule & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){2-4} \cmidrule(l){5-7} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{48.5}{\cm}$ (training) & 7.617 & 2.1278 & 25.531 & 2.165 & 0.1718 & 5.467 \\ $l_{0}=\qty{41.0}{\cm}$ (training) & 9.694 & 2.1149 & 21.076 & 1.084 & 0.0277 & 1.766 \\ $l_{0}=\qty{45.7}{\cm}$ (validation) & 8.803 & 1.9380 & 26.264 & 1.368 & 0.0468 & 3.614 \\ $l_{0}=\qty{40.6}{\cm}$, tendon (validation) & 6.373 & 0.5935 & 16.859 & 1.611 & 0.0427 & 5.891 \\ \midrule Mean (training) & 8.656 & 2.1213 & 23.303 & 1.624 & 0.0998 & 3.616 \\ Mean (validation) & 7.588 & 1.2657 & 21.562 & 1.490 & 0.0447 & 4.752 \\ \bottomrule \end{tabularx} \end{table}
to \begin{table}[H] \caption{Pinned-knee flexor torque-prediction goodness of fit: the improved BPA characterization alone (the zero-correction prediction of the text) and the corrected model at the adopted correction terms (Table~\ref{tab:AdoptedXi}); training and validation tests follow the adopted cross-validation fold. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:FlxPinGof} \begingroup \renewcommand{\arraystretch}{1.4}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc cccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){3-5} \cmidrule(l){6-8} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule \multicolumn{8}{@{}l}{\textbf{Training}} \\ 48.5 & -- & 7.617 & 2.1278 & 25.531 & 2.165 & 0.1718 & 5.467 \\ 41.0 & -- & 9.694 & 2.1149 & 21.076 & 1.084 & 0.0277 & 1.766 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 8.656 & 2.1213 & 23.303 & 1.624 & 0.0998 & 3.617 \\ \midrule \multicolumn{8}{@{}l}{\textbf{Validation}} \\ 45.7 & -- & 8.803 & 1.9380 & 26.264 & 1.368 & 0.0468 & 3.614 \\ 40.6 & 35 & 6.373 & 0.5935 & 16.859 & 1.611 & 0.0427 & 5.891 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 7.588 & 1.2657 & 21.562 & 1.490 & 0.0447 & 4.752 \\ \bottomrule \end{tabular} \endgroup \end{table}
```

Containing section: **Flexor on Pin Joint**.

> \begin{table}[H]

### D-I09 — Insertion — line 431

```text
bbolen83
23 October, 11:01 am
Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-extensor goodness of fit]{Optimized pinned-extensor goodness of fit at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.} \label{tab:ExtPinGof} \begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}} \toprule & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){2-4} \cmidrule(l){5-7} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{43.6}{\cm}$ (training) & 2.437 & 0.9618 & 4.284 & 0.978 & 0.1634 & 1.818 \\ $l_{0}=\qty{43.2}{\cm}$, tendon (training) & 2.878 & 4.7304 & 4.817 & 1.525 & 1.3270 & 2.481 \\ $l_{0}=\qty{45.7}{\cm}$ (training) & 2.178 & 1.0586 & 4.974 & 0.490 & 0.0535 & 1.206 \\ $l_{0}=\qty{40.0}{\cm}$ (validation) & 3.092 & 0.5638 & 4.842 & 1.244 & 0.1303 & 2.196 \\ $l_{0}=\qty{40.0}{\cm}$, tendon (validation) & 3.622 & 6.3919 & 7.845 & 1.350 & 0.8878 & 1.890 \\ $l_{0}=\qty{46.5}{\cm}$ (validation) & 0.949 & 0.3152 & 2.272 & 0.699 & 0.1714 & 1.559 \\ \midrule Mean (training) & 2.498 & 2.2502 & 4.692 & 0.997 & 0.5146 & 1.835 \\ Mean (validation) & 2.554 & 2.4236 & 4.986 & 1.098 & 0.3965 & 1.882 \\ \bottomrule \end{tabularx} \end{table}
to \begin{table}[H] \caption{Pinned-knee extensor torque-prediction goodness of fit, at the adopted correction terms with the effective-stiffness pair locked to the flexor-identified values (Table~\ref{tab:AdoptedXi}): the improved BPA characterization alone (the zero-correction prediction of the text) and the corrected model. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:ExtPinGof} \begingroup \renewcommand{\arraystretch}{1}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc cccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){3-5} \cmidrule(l){6-8} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule \multicolumn{8}{@{}l}{\textbf{Training}} \\ 43.6 & -- & 2.437 & 0.9618 & 4.284 & 0.978 & 0.1634 & 1.818 \\ 43.2 & 22 & 2.878 & 4.7304 & 4.817 & 1.525 & 1.3270 & 2.481 \\ 45.7 & -- & 2.178 & 1.0586 & 4.974 & 0.490 & 0.0535 & 1.206 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 2.498 & 2.2503 & 4.692 & 0.997 & 0.5146 & 1.835 \\ \midrule \multicolumn{8}{@{}l}{\textbf{Validation}} \\ 40.0 & -- & 3.092 & 0.5638 & 4.842 & 1.244 & 0.1303 & 2.196 \\ 40.0 & 40 & 3.622 & 6.3919 & 7.845 & 1.350 & 0.8878 & 1.890 \\ 46.5 & -- & 0.949 & 0.3152 & 2.272 & 0.699 & 0.1714 & 1.559 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 2.554 & 2.4236 & 4.986 & 1.098 & 0.3965 & 1.882 \\ \bottomrule \end{tabular} \endgroup \end{table}
```

Containing section: **Extensor on Pin Joint**.

> \begin{table}[H]

### D-I10 — Insertion — line 486

```text
bbolen83
23 October, 9:48 am
Changed: addresses.
to ...addresses. Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of Section~\ref{sec:redesign_physical}: across its fifteen tested measurements from \ang{-119} to \ang{5}, the measured torque magnitude met or exceeded the human biceps-femoris short-head maximum at every point, with margins from \qty{0.42}{\newton\meter} at \ang{-105} to \qty{17.0}{\newton\meter} at \ang{-14}.
```

Containing section: **Flexor on Biomimetic Knee**.

> The torque of the single $\diameter$\qty{20}{\mm} flexor actuator of panel \textbf{B} did not exceed the magnitude of the maximum human isometric flexion torque between \ang{60} and \ang{100} of knee flexion, and within a $\pm 20\%$ band about the human mean it still fell short of human-level torque between \ang{78} and \ang{103} of flexion. Because operating at \qty{620}{\kPa} across the full desired range of motion would load the printed femur beyond its ultimate strength, single points were also recorded at the reduced pressures \qtylist{560;421;325;281}{\kPa}, which show that the desired range of motion is achievable at reduced pressure. This shortfall in flexion torque defines the target that the attachment-route redesign of Section~\ref{sec:ongoing} ...addresses. Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of Section~\ref{sec:redesign_physical}: across its fifteen tested measurements from \ang{-119} to \ang{5}, the measured torque magnitude met or exceeded the human biceps-femoris short-head maximum at every point, with margins from \qty{0.42}{\newton\meter} at \ang{-105} to \qty{17.0}{\newton\meter} at \ang{-14}.

### D-I11 — Insertion — line 489

```text
bbolen83
23 October, 9:50 am
Changed: \begin{center} \includegraphics[width=1\textwidth]{FlxGrp.eps} \end{center}
to \begin{center} \includegraphics[width=0.72\textwidth]{FlxGrp.eps} \par\smallskip \begin{minipage}[t]{0.30\textwidth}\vspace{0pt}\centering {\textbf{(C)}}\par\smallskip \includegraphics[width=\linewidth]{Flx20mm_meetsHuman.pdf} \end{minipage} \end{center}
```

Containing section: **Flexor on Biomimetic Knee**.

> \begin{center}

### D-I12 — Insertion — line 497

```text
bbolen83
23 October, 9:53 am
Added:Â  \textbf{(C)} The redesigned two-BPA flexor route
```

Containing section: **Flexor on Biomimetic Knee**.

> \caption{Biomimetic-knee torque versus knee angle for the flexor BPAs. Knee angle is negative in flexion, and flexion torque is negative. \textbf{(A)} $\diameter$\qty{10}{\mm}, $l_{0}=\qty{41.5}{\cm}$ actuator with a \qty{21}{\mm} artificial tendon at \qty{604}{\kPa}. \textbf{(B)} $\diameter$\qty{20}{\mm}, $l_{0}=\qty{41.5}{\cm}$ actuator with a \qty{14}{\mm} artificial tendon at \qty{620}{\kPa}. Indigo filled circles are measured torque; the magenta line is the optimized prediction of the corrected model using the flexor terms identified on the pinned knee; the dashed light orange line is the unoptimized prediction of the improved BPA characterization alone. In \textbf{(B)}, the black dotted line is the maximum human isometric torque, and the pink circles are measured points at the reduced pressures \qtylist{560;421;325;281}{\kPa} (from left to right, respectively), recorded because full pressure would load the printed femur beyond its ultimate strength. \textbf{(C)} The redesigned two-BPA flexor route}\label{fig:FlxGrp}

### D-I13 — Deletion — line 536

```text
bbolen83
23 October, 10:53 am
Deleted:Â +
```

Containing section: **Follow-Up Identification and Route Redesign**.

> Flexor & $0.00394$ & $3.998\times10^{4}$ & $1.473\times10^{4}$ & $0.158$\\

## Historical comments and Ben replies

The full, untruncated archive contains 272 resolved threads and 74 substantive messages by Ben: [historical_review_archive.md](review_audit_20261026/historical_review_archive.md). It preserves each surviving quoted anchor, author, date string, reply, and resolution event. The following index identifies every substantive Ben message; these are historical records rather than new requirements.

- **R001, 11 November, 1:55 pm:** Remove this or create another example using the 20 mm BPA?
- **R002, 13 November, 12:02 pm:** This paragraph initially had a comparison between Festo and experimental data for the 10 and 20 mm BPAs of a hypothetical 100 mm resting length. Using real data instead, I made it a comparison with just the 10 mm diameter BPA.
- **R003, 27 October, 4:37 pm:** The goal is to provide an aim towards future work on BPA characterization. Specifically, to find active, passive, and velocity effects that are a strong analogue to Thelen's Hill muscle model.
- **R008, 6 November, 6:10 pm:** Can you? I believe you created this file because I don't have it on my system. I have the original, which is 4-panel figure.
- **R009, 6 November, 6:10 pm:** Can you delete these?
- **R020, 10 November, 3:54 pm:** I commented this out. Let me know if what is in the supplementary materials is sufficient.
- **R028, 6 August 2025, 10:59 pm:** I think they will want to know what the 7.5 mm offset is from
- **R029, 6 August 2025, 3:35 pm:** I agree. The experimentally collected data has more outliers. Should I set aside 20% of the collected data and use it for validation, instead? Should I discard the equation rows for Festo?
- **R030, 4 August 2025, 7:23 pm:** I think it's confusing to add PSI here when everything else is in metric.
- **R030, 4 August 2025, 7:26 pm:** What if we say, "(originally specified as the near-equivalent 90 PSI)"?
- **R031, 2 August 2025, 7:49 pm:** I don't understand what this is a comment on.
- **R032, 4 August 2025, 7:56 pm:** What is this for?
- **R033, 2 August 2025, 7:42 pm:** model of what? -CM
- **R034, 2 August 2025, 7:42 pm:** model of what? -CM
- **R035, 2 August 2025, 7:54 pm:** embellish a bit. -CM
- **R035, 2 August 2025, 7:55 pm:** Maybe more descriptive? Like "environmental exploration".
- **R036, 2 August 2025, 7:41 pm:** present vs past tense. -CM
- **R036, 4 August 2025, 7:18 pm:** I changed to past tense. - BB
- **R037, 2 August 2025, 7:39 pm:** are there other characteristics? Might just say poor force and length characterization. - CM
- **R038, 16 July 2025, 9:49 am:** It is sometimes so repetitive to type "uninflated diameter" that I think we should keep phi.
- **R039, 16 July 2025, 9:50 am:** I see abstracts define variables all the time.
- **R040, 16 July 2025, 9:43 am:** This is the same knee joint. Steele has 2 papers in 2018. I didn't change the knee design. I don't see this figure called out in the text anymore but I believe that it did, and it referenced my previous papers and the OpenSim papers.
- **R042, 16 July 2025, 9:55 am:** All of the dots in Figure 3. It's too many to list in text.
- **R061, 13 June 2023, 11:32 pm:** Add citation for biomechanics society reference frames
- **R099, 13 December 2023, 10:04 am:** For both of these tables, the values have slightly shifted and I've also changed the gof measures to be R^2, RMSE, and largest error. This exists in an excel spreadsheet and I will upload it.
- **R100, 14 January 2024, 12:57 am:** I have to check but I don't think the yellow is your 2017 model. I think it's just what is predicted pre-optimization.
- **R103, 10 January 2024, 2:30 am:** Why isn't this important?
- **R105, 10 January 2024, 2:30 am:** Yes
- **R105, 20 January 2024, 11:33 pm:** Done with this section now.
- **R106, 10 January 2024, 2:28 am:** cite what? Is the way I cited the Festo Tool earlier acceptable?
- **R109, 10 January 2024, 2:34 am:** I thought that's why we put it in the table
- **R110, 10 January 2024, 2:24 am:** Do you mean, "our modeling work?"
- **R110, 10 January 2024, 2:24 am:** Didn't you write this? It looks like it's in green.
- **R110, 4 February 2024, 3:01 pm:** Ok I understand. Yes, I faced this same dilemma. I checked the length calculating aspect of our algorithm and the solidworks model. I think that hole and position tolerance play a factor. I also think that the endcap and bracket holes deform over repeated cycles.
- **R111, 4 February 2024, 3:02 pm:** Not in the pinned flexor group. The brackets created a large moment arm.
- **R112, 10 January 2024, 2:23 am:** Why? Doesn't Frontiers just want figures at the end, anyway?
- **R113, 27 August 2023, 6:30 pm:** But it is the validation that the optimization works with a different resting length
- **R118, 13 December 2024, 3:49 pm:** why?
- **R132, 10 January 2024, 2:30 am:** Isoclines?
- **R132, 13 January 2024, 7:30 pm:** I put the isolines in part (A) like you requested.
- **R133, 14 January 2024, 12:50 am:** Alex says not to use the 40mm data experiment since it was inconclusive. He might want to repeat the experiment using a 20mm to replace a vastus muscle. He might also want to add a patella for the extensors (this would require modifying the solid model so I can use a mate in solidworks that will allow me to know where the patella is as a function of knee angle).
- **R134, 4 February 2024, 2:53 pm:** Yes, now it is. It can all be {eq:exp} if we want and it'll make the second table messier but that's ok.
- **R135, 14 January 2024, 1:04 am:** Why are these coefficients changing? I need to know so I'm not spending time on creating latex tables that I have to go to lengthy efforts to revise.
- **R136, 1 September 2023, 3:19 am:** We need to determine if we want to use the word "contraction" or "strain."
- **R137, 13 December 2023, 9:56 am:** No longer true.
- **R141, 13 December 2023, 9:51 am:** I removed some of the description of points (in a translation matrix) to the Supplementary material section.
- **R143, 13 December 2023, 9:52 am:** Note to self: Cite Steele, and/or Bolen, and/or Delp (Opensim)
- **R145, 13 December 2023, 9:48 am:** Should I make a table of which resting actuator lengths were tested at each diameter and each load cell?
- **R146, 25 August 2023, 10:13 pm:** I don't exactly know what you mean here.
- **R148, 13 June 2023, 9:16 pm:** Add citations from Dr. Hunt or other NeuroNex researchers?
- **R148, 25 August 2023, 10:11 pm:** Do any of my references currently in Zotero stand out to you as ones I can use here, or do I need to spend time to re-read them?
- **R150, 25 August 2023, 10:33 pm:** Should we not note this?
- **R151, 25 August 2023, 10:32 pm:** How should I specify that the load cell was different for different resting lengths, as I note in the commented paragraph below?
- **R152, 13 June 2023, 11:09 pm:** One word? Should I include description and figure in supplementary section?
- **R154, 13 June 2023, 11:39 pm:** Add a reference to supplementary table with joint and muscle locations?
- **R159, 8 December 2023, 12:37 pm:** verify
- **R163, 13 June 2023, 10:57 pm:** Delete whole paragraph?
- **R165, 10 August 2023, 5:21 pm:** I strongly disagree.
- **R167, 18 August 2023, 1:51 pm:** Compare against two virtual models and use the new data to create an updated robot model.
- **R170, 13 June 2023, 11:05 pm:** include Cody's Master's thesis? Uses shorter BPAs I believe.
- **R171, 13 June 2023, 10:56 pm:** add citation or footnote for gait2392?
- **R172, 13 June 2023, 9:19 pm:** We use Festo Corp. BPAs because...
- **R174, 13 June 2023, 11:26 pm:** add citation or footnote with url?
- **R176, 13 June 2023, 9:02 pm:** Which should I use? Phi looks nicer, but is a little harder to code.
- **R182, 15 December 2022, 3:36 pm:** Isn't it important to transition from Intro to Methods and Materials?
- **R188, 15 December 2022, 4:00 pm:** But the 20mm and 40mm tests weren't done at full pressure. Only two data points for the 20mm were done using this pressure. The rest weren't, because we were breaking brackets and didn't want to miss the deadline. These larger BPAs aren't allowed (by Festo) to be at full pressure unless they have some contraction b/c their force will exceed the max theoretical force and permissible operation range.
- **R190, 15 December 2022, 3:41 pm:** We didn't really track which tests used which. It was more of a space consideration, because the analog ones come in various lengths, and the digital one is long. I can try to look through my pictures or take my best guess.
- **R191, 15 December 2022, 3:46 pm:** Almost all 10mm ID actuator lengths had at least one test recorded recorded with crane scale and the other with a load cell. Should I change my figures and add MarkerEdgeColor to the data points that used the crane scale, and then add a note here that says the difference is noted in the figures in the Results section?
- **R232, 25 March 2022, 2:31 pm:** Where they wrong in this publication as well?
- **R234, 27 March 2022, 5:23 pm:** yes
- **R252, 26 March 2022, 4:31 pm:** Done, except instead of saying flexibility I might want to say something similar but make it clear we mean curving/bending.
- **R254, 26 March 2022, 4:28 pm:** Yes it's how their template is set up and how one of Dr. Hunt's previous Frontiers papers has it too.
- **R255, 26 March 2022, 4:27 pm:** I think this code just makes it say "Bolen, et al."
- **R266, 25 March 2022, 2:34 pm:** So after I talk about this in discussion then shall I include extra graphs that compare theoretical calculations with and without these results taken into account?

## ZCode delivery and verification

Use [ZCODE_REPLOT_REQUESTS.md](review_audit_20261026/ZCODE_REPLOT_REQUESTS.md). Return bounded `.tex`/figure changes under the existing machine/batch review inbox. Include exact generating script, source-data identifiers, parameter row/fold, plotted series export, and a rendered preview. Do not overwrite live chapters or accept/reject tracked changes.

Evidence files: [raw document diff](review_audit_20261026/torque_documents_raw.diff), [paragraph alignment](review_audit_20261026/paragraph_alignment.json), [figure/table manifest](review_audit_20261026/figure_table_manifest.json), [paper review items](review_audit_20261026/paper_review_entries_with_positions.json), [dissertation review items](review_audit_20261026/dissertation_torque_review_entries.json).

The diff is a complete source diff; structural/template differences are expected because one file is a standalone article and the other a dissertation chapter. The curated differences above identify the actionable scientific and editorial issues. No current Overleaf asset hash comparison was obtained; matching local filenames are generator leads, not proof of byte-identical live images.

Final live compilation: dissertation Errors 0 / Warnings 0; torque paper Errors 0 / Warnings 20 (duplicate-label and font-substitution warnings listed in the saved compile evidence). Both full Overleaf documents compiled after the bounded prose edits. No local model or plotting run was performed.
