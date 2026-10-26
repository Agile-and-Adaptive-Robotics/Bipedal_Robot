# Review of the three incoming batches

The files in this folder are proposed, bounded patches. No Overleaf text, tracked change, or advisor comment was altered. No chapter was replaced. Original submissions remain in their machine folders; `incoming_manifest.json` records their hashes. The patches require anchor reconciliation against the live document before application. A clean fragment proof does not establish a clean dissertation build.

## Application order and anchors

1. **Tools inventory:** use `patches/P01_replacement.tex` and `P01_companion_table.tex` together, replacing the existing opening tools/contributions inventory and its capability table. Both use the approved framing: Researchers, Gait2392 second, tools created, and selected demonstrations. The table retains `tab:sns_capability`. Do not reinsert a second table or describe Rybak as a Simulink implementation. The SolidWorks assembly remains author-confirmed; its inaccessible saved model does not prevent listing the tool.
2. **Conversion:** `P07_conversion_summary.tex` replaces the conversion discussion within `sec:sns_existing_tools`; `P06_appendix_conversion_section.tex` supplies `sec:app_conversion_repairs` in Appendix D. These are the previously reviewed short narrative and detailed repair account, rather than the incoming single long paragraph. Preserve nearby unrelated material.
3. **Figure 6.1/history:** the copied `P04_figure61_block.tex` and `P05_controller_history.tex` remain available for the already agreed correction. Check whether these have already been applied. The original Connor–Ben diagram is the shared-RG architecture; later bilateral RG and sensory extensions are separate. Do not restore Ben's deleted suspended-stepping qualifier.
4. **W2L tests:** `EG_controlled_tests.tex` replaces only `sec:w2l_perturbation` through the start of the next subsection/section. Preserve the existing section and figure labels. Upload `w2l_perturbation_fully_applied.pdf` under that exact name. The mass-condition figures are 32–38 strides, minimum pelvis 0.8981 m, and nine tilt failures, rather than 33–38, at least 0.90 m, or robust balance. Continued stepping and survival are different outcomes.
5. **Ablation:** insert `EB_feedback_ablation.tex` after the W2L perturbation figure, if the live document does not already contain this comparison. It is contact and Ib-force ablation on the standard muscle-actuated W2L model, with Ia/II retained. The scored interval is 25.499 s within a 30 s simulation. Same harness configuration does not mean the same measured support force. Saved metrics support the comparison; no preserved raw trajectory was supplied for it.
6. **Deng:** `DENG_original_engine.tex` is an insertion in `sec:preliminary_sim_results`, immediately after the existing original-model/original-engine discussion, before any extension/boundary paragraph. Reconcile that adjacent paragraph to avoid duplicate qualifications. It concerns the archived antagonist model in AnimatLab, not quantitative validation of the Simulink port or the 2023 biarticular model. Use the repaired `deng_original_reviewed.pdf`. Merge the separate bibliography entry only if its key is absent. Existing bibliography bytes were not changed.
7. **Torque discussion:** `EG_pulley_discussion.tex` replaces the body of `sec:disc_pulley` through the next subsection/section. The force/travel ratio applies specifically to the two-part reversed tackle. The two optimized routes are not an isolated pulley A/B experiment.
8. **Merged extensor figure:** `EG_extensor_merged_figure.tex` replaces the two figure environments carrying `fig:Ext10mm_52cm` and `fig:Ext_passthru`. It actually draws A/B letters and preserves both labels as aliases. Change prose referring to the pass-through photograph to “Figure~\ref{fig:Ext10mm_52cm}B”; torque-plot references should identify panel A. Check both old-label references before applying.
9. **Captionless physical figures:** `EG_physical_figures_captionless_OPTION.tex` is conditional on retaining Ben's caption-removal instruction recorded by the incoming batch. It replaces only the two physical-final figure environments and explicitly advances the figure counter. Bare labels without captions do not advance it. These figures still need marker-shape differentiation before final presentation. The FlxGrp shared-axis asset is staged separately; do not confuse it with the distinct three-panel 20 mm validation figure `Flx20mm_56cm`.

## Decisions by incoming file

| Batch/file | Decision | Reason |
|---|---|---|
| EB475WS4 P01, knee reflex | Caption correction ready; figure layout repair pending | Simulink settles at 44.4 degrees, NumPy at 43.0, maximum deviation 3.6. Supplied “both near 43.5” is inaccurate. Circuit clipping remains. |
| EB475WS4 P01, knee CPG | Hold | The displayed interval rises without a visible complete cycle; the motion strip does not establish the caption's one-cycle claim. |
| EB475WS4 P01, loaded arm | Useful selected demo; layout repair pending | Use direct sag: 2.00 versus 7.85 degrees. Do not turn the angle ratio into an independently validated stiffness measurement. |
| EB475WS4 P01, ablation | Rewritten | Partial sensory ablation, not total deafferentation; standard muscle model, not BPA; scored-window and support distinctions added. |
| easteregg2 P01 conversion | Superseded by reviewed P07/P06 | Separates driver-joint limit restoration from removal of coupled path-point limits; avoids an unreadable single paragraph. |
| easteregg2 P02 capability | Do not apply | Old framing, unsupported aggregate claims, and malformed doubled TeX commands. Approved tools inventory retained instead. |
| easteregg2 P03 W2L tests | Rewritten | Commissural deletion retains 13/12 strides; Ib deletion alone stops both legs at this operating point. Old six-field search results cannot be silently relabeled as newly optimized all-field results. |
| easteregg2 P04 pulley | Accepted with specific ratio wording | Halving force/doubling travel is the two-part reversed-tackle case. |
| easteregg2 P05 merge | Corrected | Actual A/B letters and both reference labels retained. |
| easteregg2 P06 physical figures | Conditional, corrected | Counter advancement restored; marker repair still needed. |
| easteregg2_2 P01 original Deng | Rewritten | Removes “unmodified” ambiguity, scopes the original engine/model, and separates 0.1 s paper-timed from 0.5 s extended stimuli. |
| easteregg2_2 P02 bibliography | Completed separately | Publisher metadata verified against the DOI page. |
| easteregg2_2 P03 Deng figure | Rebuilt from saved raw charts | Removes printed TeX, overlaps, tiny labels, and misleading universal 2.0–2.5 s shading. Solid/dashed traces provide a non-color distinction. RG axes use their own voltage ranges. |

## Evidence and limits

- Original raw Deng chart files are preserved under `evidence/deng_raw`; the new plot uses these saved values and does not run the model. Incoming summary statistics, phase method, export, and supervisor verdict remain in the original easteregg2_2 batch. The quoted phase uncertainty is the saved analysis's post-stimulus standard deviation, not an experimental confidence interval.
- W2L numerical corrections are from the previously audited fully applied-field replay and body-change evidence preserved with Revision E. This review does not rerun controller optimization or simulations.
- Ablation source: `Code/Matlab/SNS_Simscape/results/ABLATION_PAIR_REPORT_20261026.md` and its paired JSON metrics. Standard W2L actuator definitions are in `Code/MuJoCo_SNS/spinal/w2l_mujoco/w2l_ground.xml`.
- Bibliography primary source: [Springer chapter record](https://link.springer.com/chapter/10.1007/978-3-031-39504-8_2).

## Repeated-error scorecard

These are review failures to correct in our workflow, not additional experiments demanded from Ben. The chapter's tools-and-selected-demonstrations framing remains intact.

| Repeated issue | Responsibility | Required check |
|---|---|---|
| Calling retained stepping robust walking/recovery | ZCode drafts; ChatGPT if accepted | Read stride, support, height, and tilt gates separately before writing the claim. |
| Reviving unsupported commissural necessity and zero-effect compensation | ZCode batch; ChatGPT review gate | Match prose and caption to corrected end-to-end parameter forwarding and actual deletion counts. |
| Replacing user-approved framing or original circuit history | Both systems | Compare against Ben's explicit rulings before proposing a patch. |
| A caption asserts content the plotted panel does not show | ZCode figure generation; ChatGPT visual review | Open the rendered figure; verify cycle, model identity, panel letters, units, and stimulus windows. |
| Treating cross-implementation agreement as published-experiment reproduction | Both systems | Name the model, implementation, stimulus, comparator, and scope in the result itself. |

The prior count-duplication, parameter-forwarding, and mislabeled joint-range bugs remain invalidated evidence, not compensated pathways. Their historical correction must be retained in the report.
