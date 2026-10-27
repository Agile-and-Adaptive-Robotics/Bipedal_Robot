<!-- TWO_SESSION_UPDATE_20261027_START -->
# ChatGPT sessions: 26–27 October 2026

Updated 27 October 2026, 14:32 Pacific, on **easteregg2**. Ben's final session numbering is **Session 1 = Complete MuJoCo Neural Tasks** and **Session 2 = Review dissertation and torque paper**. Session 2 prepared this consolidated entry from its own work, Session 1's saved recovery report, and Session 1's latest visible progress messages. Session 1 is still active; its progress statements below are a dated snapshot, not a declaration that its task is finished. Earlier 26 October entries below remain intact and retain their original attribution and verification records.

## Session 2 — dissertation and torque-paper review

### Completed on 27 October

- Restored the paper's bracket-deflection development: force components, compliance matrix, displacement, projection, effective stiffness, and brackets/tendon in series. Restored contact/release calculations in the supplement and their main-text pointer.
- Created separate vector schematics for correction-factor projection and BPA bending across a flexed knee. The latter shows the actuator, attachments, tendon, radius, and full alpha contact angle. These appear as dissertation Figures 4.5/4.6 and paper Figures 5/6. PDF, SVG, PNG, and alternative-text assets are saved in both figure trees, with a reproducible generator.
- Renamed the dissertation subsection to “Identification of the Torque-Model Correction Terms.” Removed obsolete 48.0/47.9 cm exception language from the adopted identification account while preserving the historical inventory in Appendix C. Corrected stiffness/frame descriptions against executable evaluators; adopted parameters and fold selections were preserved.
- Corrected the Morrow/Bolen placement objective's sign and normalization against `Code/Matlab/Mesh_Optimization/costFunctionKnee.m`, added its penalty/weight/meaning table, and moved its citation out of Background into the torque-methods account. Retained the subsequent stiffness-aware cost evolution. Corrected interpretation of identification metrics, root-finding tolerance, and empirical bending losses.
- Reconciled the capability comparison with the saved 27-site CSV: cell coverage 37.28%/65.19%, axis means 38.39%/66.17%, mean torque-ratio scores −1.93/−3.63, and 23 improvements/four ties. Clarified transmission/sizing changes, baseline diameter substitutions, and wrong-sign torque scores. Added the 27-site DoF/diameter inventory to Appendix C and the supplement.
- Corrected the actuator-count conclusion: no saved tested combination meets all five angles; the best saved three-BPA/two-part-tackle case meets four. Removed the unsupported causal +4.51 N·m redirect claim. Did not restore the rejected 40.2 N·m clause or apply the stale Revision_C attachment swap. **pA remains origin; pB remains insertion.** Measurements, fitted parameters, and result mats were not altered, and no optimization was rerun.
- Corrected abstract/conclusion and validation descriptions, including 15 flexor and nine extensor measured points. Retained **17 main-paper display items (12 figures, five tables)** under Ben's instruction to consolidate later if needed. Repaired affected table pagination, oversized figures/tables, and the paper's empty template ending.
- Built and visually checked the affected local outputs: **paper 41 pages**, **supplement eight pages**, and **dissertation 320 pages**. Final local builds have no errors, undefined references/citations, duplicate labels, or oversized floats; inherited font/typesetting warnings remain. This verifies the reviewed torque material and targeted layout repairs, not every claim in every dissertation chapter.
- Uploaded changed sources and schematic PDFs to both Overleaf projects and recompiled. The [torque paper](https://www.overleaf.com/project/622e9536b8ea1e0a7f167cf5) compiles with **zero errors and zero undefined-reference warnings**; 19 font-substitution warnings remain. Refreshed its uploaded auxiliary files so supplementary references resolve, and left `Main.tex` selected. The supplement's independent live build was not separately verified; its eight-page local build was verified.
- The [dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb) compiles with **zero errors** after adding the live `figs/30-torque/` search path to Appendix C. Existing live Appendix E content was preserved. Existing comments and tracked-change dispositions were preserved. Native screen control stopped at Ben's request; subsequent checks used the Chrome extension.

### Outstanding / deferred

- The live dissertation still defines `fig:sadb_explorer`, `tab:sadb_afferents`, and `tab:sadb_pathways` in both Chapter 5 and Appendix E. These produce three duplicate-label warnings plus the summary warning and can misdirect references. They were explained to Ben but **have not been fixed**. The local review build lacks that extra live appendix content, so its clean label check does not clear the live conflicts.
- Consolidating the paper's 17 display items is deferred. Remaining font/typesetting warnings have not all been eliminated. No full dissertation-wide scientific audit was claimed.
- No commits or pushes were made. Existing ZCode and unrelated session changes were preserved for Ben's review.

Evidence and commit-review map: `Documentation/Reports and Papers/Knee_Torque_Test/CHATGPT_REVIEW_20261027.md`. Per-source diffs, verification data, and upload archive: `Documentation/Reports and Papers/Knee_Torque_Test/Review_20261027/`. Local builds: `Knee_Torque_Test/build_review_20261027/` and `Dissertation/ProofFinal/build_r12_torque_review_20261027/`. Figure generator: `Knee_Torque_Test/Figures/Figure components/draw_torque_schematics_20261027.py`. Source snapshots and Overleaf confirmation images: `tmp/torque_review_20261027/`. All these paths are under `Documentation/Reports and Papers/` except the explicitly named root `tmp/` path.

## Session 1 — MuJoCo/neural-model recovery and CAD robot

### Completed and recorded in the recovery report

- Repaired W2L gain registration so saved effective settings reach the compiled network; unknown settings now fail explicitly. Audited 56 effective conductance dimensions and three runtime sensor dimensions, excluding two inert saved settings.
- Verified supported stepping in fresh processes. Trial 458 produced 28 left/21 right contact events, median phase 0.496063, 2.246231 m displacement, and 0.091786 m/s fitted speed over the report's analysis interval. Four independent replays matched; both reported candidates passed their separate 12-condition perturbation checks. These are **supported stepping**, with a vertical harness still active even at `wsup=0`; they do not establish harness-free or physiological walking.
- Added an exact W2L neural backend with full-trajectory parity. Its isolated neural step improved by 1.425×; this is not a whole-simulation speedup claim. Rejected the previous sparse backend's closed-loop parity. A separately compiled Di Russo-inspired adaptation matched its own equations, without establishing equivalence to the original author's controller.
- Recovered the original Di Russo repository, custom SCONE fork, model, and parameter structure. Corrected the previous source attribution and identified the original reduced nine-active-coordinate/18-muscle Thelen model. Rejected the surrogate search's somersault candidate as walking and stopped that optimization route.
- Regenerated W2L perturbation/feedback panels with effective settings and recorded their interpretation limits. Imported the AH001 CAD geometry, preserved both four-bar knees, added the adjustable lumbar hinge/spine/payload, and mechanically verified a 75 kg model. The source CAD was opened read-only and not saved during that recorded import.

Evidence: `Code/MuJoCo_SNS/spinal/reports_20261027_codex/REPORT.md`, its verification JSON files, and `ZCODE_HANDOFF.md`; CAD build and mechanical checks: `Code/MuJoCo_SNS/robot_torso/README.md` and its audit files. That saved recovery report describes an earlier 92-route scaffold; the newer work below supersedes that scaffold's route-count and frozen-point status.

### In progress at this update

Session 1's latest visible messages report correcting the coordinate frame to **x anterior, y superior, z right**, refreshing Ben's knee pose, and using the robotbody 27-actuator-per-side selection mirrored left (**54 selected actuators**). They report no reversed neutral actions after attachment corrections (largest shift about 21 mm), moving-point functions matching OpenSim numerically, and conditional gates passing source checks. These are Session 1's current progress reports; Session 2 has not independently rerun those checks.

The SNS network is now reported to drive actuator forces. Initial supported tests collapsed rather than balanced; Session 1 replaced the force allocator with a bounded solver and is testing standing before assessing stepping. Small length jumps at semitendinosus conditional-guide switches remain flagged in the full source-route audit; those routes are outside the selected 27-actuator map. **Supported balance and neural locomotion of the new CAD robot are not yet demonstrated.** Foot integration and the final routing/controller verification remain part of that active work.

### Still incomplete / blocked

- **Di Russo walking replication is incomplete.** Native original-controller replay requires the author's custom SegmentalController runtime, absent from the installed environment; the original-model baseline remains outstanding.
- **Actual CAD robot neural locomotion is incomplete**, despite verified mechanics and the newer active coupling work.
- **Task 3 has not been launched**, pending the requested completion of tasks 2 and 5. Earlier recovery notes and handoff prompts are evidence of preparation, not completion.
- Session 1's recorded recovery work made no commits, pushes, or dissertation edits. Subsequent changes are still active and should be reviewed against its next saved report before committing.

<!-- TWO_SESSION_UPDATE_20261027_END -->

---

<!-- ETHICS_CLOSING_20261026_START -->
# Follow-up: ethical-case paragraph closing

Ben requested a replacement for “The failure modes of these experiments damage hardware, not tissue.” Applied a one-sentence tracked replacement in live `15-background.tex`, Section `sec:bg_ethics`, local snapshot line 230:

> Robotic models therefore offer both experimental repeatability and a way to reduce reliance on the invasive animal preparations that established many of the neural mechanisms studied here.

The replacement closes the paragraph around repeatability and reduced reliance on invasive animal preparations, then leads into the next paragraph's discussion of how those preparations appear in circuit schematics. It makes no claim that animal use has already been quantitatively reduced. No numbers, citations, other prose, comments, or existing tracked-change dispositions were changed. Exact whole-file comparison passed. Overleaf Errors 0 / Warnings 0. Full local build `build_astra_05` passed with no errors, undefined references, or missing figures; existing overflow messages remain.

Before/after source and exact replacement: `LIVE_ASTRA_20261026/followup_ethics_closing/`.
<!-- ETHICS_CLOSING_20261026_END -->

---

<!-- CAT_TRANSITION_20261026_START -->
# Follow-up: transition from cat schematics to experimental procedures

Updated 2026-10-26T13:56:48-07:00. Ben asked that the paragraph transition explain how apparently normal or happy-looking cat schematics differ from the experimental preparation. Applied in live Overleaf, `15-background.tex`, Section `sec:bg_ethics`, local snapshot line 232.

Restructured the first paragraph around the schematic's appearance, the brainstem-transection preparation, the archival figure, and the specific meaning of replacement. It describes the impression a drawing can give without claiming to know the animal's emotional state. Retained the existing source keys for the preparation and footage; removed the unsupported “almost exclusively in schematic” generalization and the vague accounting/counterweight wording.

New first paragraph:

> Schematics of decerebrate-cat locomotion emphasize the neural circuits under study, but often leave the surgical preparation out of view. A simplified drawing of a walking cat can make the animal appear intact, alert, or even contented. That appearance should not be confused with the condition of the experimental animal: a brainstem transection has disconnected the forebrain, while neural circuits below the transection can still generate treadmill stepping \citep{brown_intrinsic_1911, severin_1967}. Figure~\ref{fig:decerebrate_cat}, taken from archival footage, shows the physical preparation behind the circuit diagrams \citep{decerebrate_cat_footage}. This distinction makes the ethical meaning of replacement concrete: the robotic experiments discussed here are intended to reproduce selected neural mechanisms without repeating these invasive preparations in living animals.

New opening of the next paragraph:

> The contrast between an intact animal and the experimental preparation is particularly clear in publicly available teaching footage.

The rest of the procedural paragraph, the figure, its caption, and the citations remain as they were in the freshly captured live source. Two exact text selections were replaced in Reviewing mode. Whole-file readback confirms these were the only changes. No comments resolved or tracked changes accepted/rejected. Overleaf recompile: Errors 0 / Warnings 0. Local full build `build_astra_04` completed with no errors, undefined references, or missing figures; existing overflow messages remain.

Before/after source and exact replacements: `LIVE_ASTRA_20261026/followup_cat_transition/`.
<!-- CAT_TRANSITION_20261026_END -->

---

# Follow-up clarity correction: feed-forward and feedback control

Applied in live Overleaf, `15-background.tex`, Section `sec:bg_neuromech`, local snapshot line 204. Ben identified the sentence beginning “Hybrid work makes the complementarity…” as unclear. Replaced that sentence with three explicit sentences:

> \citet{haeufle_benefit_2018} combined pre-recorded muscle activation patterns (feed-forward control) with sensory reflex responses (feedback control) in a simulated human walking model. They tested whether the model could keep walking after an unexpected drop in ground height. The benefit of combining the signals depended on the muscle being controlled and on the proportion of activation supplied by each signal.

The revision identifies the two signals, defines the step-down test, and states what the muscle-dependent finding means. It was verified against the [original paper](https://www.frontiersin.org/journals/computational-neuroscience/articles/10.3389/fncom.2018.00080/full). The same citation key is retained; no numerical result was added or altered. Whole-file readback confirms this was the only content change in the chapter. Applied in Reviewing mode; no comments resolved or tracked changes accepted/rejected.

Verification: Overleaf Errors 0 / Warnings 0. Full local build `LIVE_ASTRA_20261026/revised_project/build_astra_03/` passed all four compilation steps with no errors, undefined references, or missing figures. Existing overflow messages remain. Before/after source, the exact replacement, and a screenshot are in `LIVE_ASTRA_20261026/followup_hybrid/`.

---

<!-- LIVE_ASTRA_20261026_START -->
# Live Overleaf revision — completed 26 October 2026

Ben explicitly authorized live editing. Applied **24 bounded text replacements plus the detailed conversion-appendix insertion**, and uploaded two reviewed figures. Changed chapters: introduction, background (one consistency paragraph), SNS, overall discussion, conclusion, and Appendix D. **The peer-reviewed isometric-force chapter was not edited.**

The SNS chapter now states the knowledge gap, names the created tools, and explains what the selected demonstrations establish. “Researchers” and Gait2392-second ordering are applied. The tools table is explicitly referenced. Figure 6.1 has visible A/B labels and the correct Connor–Ben shared-RG diagram. The conversion account incorporates Ben's robotbody/small-muscle explanation. The appendix now explains the Python Rybak resetting/non-resetting mechanism.

All edits were made in Reviewing mode. Existing tracked changes were not accepted/rejected; comments were not resolved, and no replies were posted. Ben's deleted suspended-stepping qualifier remains deleted. Whole-file readback verified all six edited files against the prepared text.

**Verification:** full local build passed with zero errors, undefined references, and missing figures. Final Overleaf compile showed Errors 0 and Warnings 0. Its 268-page PDF is saved. Box-overflow messages remain; the record does not certify every page or every advisor thread as polished.

Delivery folder: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026/`
- `APPLIED_REVISIONS.md`: exact patch purposes, current local line anchors, advisor/Ben response matrix, evidence, and remaining limits.
- `Dissertation_Overleaf_revised.pdf`: final live compiled PDF.
- `LIVE_CHANGES.diff`, `before/`, `after/`, `FINAL_VERIFICATION.json`: source and verification record.
- `overleaf_revised_gap.png`: live browser confirmation.

Do not reapply the old incoming batches or replace whole chapters. The new Deng and additional W2L perturbation proposals remain separate; the held CPG/other figure repairs were not imported. No simulations, commits, pushes, or outgoing ZCode messages were performed.
<!-- LIVE_ASTRA_20261026_END -->

---

# Three-batch review: reconciled packet, 26 October 2026

Reviewed EB475WS4, easteregg2, and easteregg2_2. Delivery: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/REVIEW_INBOX/REVIEWED_20261026/REVIEW_DECISIONS.md`. Use its bounded patches and application map; do not apply the incoming batches verbatim. `ZCODE_REPAIR_NOTES.md` identifies held figures and exact repairs. `REVIEWED_FRAGMENTS.pdf` is the standalone proof, not a chapter replacement.

Approved tools framing retained (Researchers; Gait2392 second). Corrected stale W2L necessity/mass claims, scoped the contact/Ib ablation, repaired the original-engine Deng figure from preserved raw charts, completed a separate publisher-verified bibliography fragment, and corrected torque panel lettering/label handling. Knee CPG figure is held because the plotted interval does not show the claimed cycle. Knee-reflex and loaded-arm graphics need clipping/annotation repairs. Repeated failures and ownership are recorded in the packet's self-scoring matrix. Prior invalidated counts/no-op deletions remain invalidated.

Fragment proof compiles cleanly; external chapter references are explicitly mocked. No Overleaf edits, active chapter replacements, bibliography rewrites, model runs, commits, pushes, or messages to ZCode were performed. Application must reconcile the live text, advisor comments, and Ben's tracked deletions first.

---

# P01 author confirmation — coupled CAD assembly

Ben confirms the SolidWorks-derived robot assembly with a neural controller and BPAs using his force model; its saved model is on the temporarily unavailable laptop. Include it as an author-confirmed contribution. P01 uses representative proof-of-concept examples, not exhaustive results for every tool. Do not block the prose on further experiments or ask Ben again for this artifact now. Rybak is Python. Details: Revision_E/review/targeted_reconciliation/tools_inventory/AUTHOR_CONFIRMATION.md.

---

# P01 reframed: tools created and representative demonstrations

Ben rejected the incomplete-work framing. The current P01 fragments list delivered model/tool contributions and functional-verification/proof-of-concept demonstrations, with credit to the reused packages. Evidence and naming boundaries are in Revision_E/review/targeted_reconciliation/tools_inventory/TOOLS_AND_DEMONSTRATIONS.md. The paragraph and capability table are paired targeted replacements, not chapter uploads. No Overleaf edits or new model runs were made.

---

# Application rule — targeted patches, not chapter replacement

The live Overleaf document is the baseline. Do not upload/paste Revision E chapter files wholesale: they lack current tracked deletions, changes, and comment anchors. Use `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/targeted_reconciliation/OVERLEAF_TARGETED_PATCHES.md`. It contains exact-text replacements, comment/date evidence, and the current table/deletion decisions. Ben’s deleted sentence must stay deleted. Supply modeling evidence separately from prose patches; do not launch parallel chapter rewrites. No Overleaf changes or inter-agent messages were made by this reconciliation.

---

# Figure 6.1 and MyoConverter detail correction — 26 October 2026

Revision E now contains Ben’s supplied original Connor–Ben shared-RG diagram as labeled Figure 6.1(B), with the mechanical walker labeled (A). Later per-side RG and RG/PF/contact-feedback changes are attributed to Ben and distinguished from the original model. The heading is exactly “Block Programming Neural Pathways.” New Appendix D.1 consolidates the conversion repairs; the main summary and existing Appendix A force comparisons remain. Details and line evidence: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/figure61_correction/CORRECTIONS.md`. Editable chapter and appendix sources remain in Revision_E/chapters. No Overleaf edits, commits, or pushes.

---

# Revision E submission clarity pass — 26 October 2026

The local PDF now leads with the concrete SNS framework: Gait2392 and separate BPA models, neural-controller alternatives, adjustable simulated harnesses, staged tuning, neural stimuli, external perturbations, recorded responses, and bounded Optuna searches. It makes no solver-superiority or global-optimum claim. Existing results, including posture failures and the scoped Ib finding, remain qualified.

Seventeen principal prose changes across six chapters align the Introduction, Background, Database, SNS, Discussion, and Conclusion. The recorded clarity-related comments have line-based treatments in `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/clarity_pass/comment_compliance.json`; `CLARITY_PASS.md` is the readable review. Earlier texts are preserved. Acceptance remains for Ben/advisor to judge. No Overleaf changes, model runs, commits, or pushes.

---

# Revision E follow-up — corrected remote evidence, 26 October 2026

The latest E supersedes the earlier “corrected evidence unavailable” status below. Read-only easteregg2 retrieval and independent raw-log recount confirm corrected control 17/16, commissural deletion 13/12, and Ib deletion 0/0. The runner forwarding repair changes the operating point; it does not retroactively make the old six-field search a per-synapse optimization. Nine of ten corrected mass-change cases fail tilt, and both push conditions exceed the 45° fall threshold. E now states those limits and restores the corrected figure.

Rybak's selected neural-behavior tests (13/13) are included with activity-based-reduction caveats. Deng and Di Russo experiment-level reproduction remains open. The human wean removes added hand drive, not all descending drive or all feedback. Ben does not want to request the never-supplied Hunt/Fletcher libraries; that is no longer an outstanding request to him.

Updated source, PDF, report, handoff, snapshot hashes, independent count audit, and self-scoring matrix are under `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/`. Earlier versions are preserved in its review folder. No Overleaf edits, remote edits, new model runs, commits, pushes, or external messages.

---

# Revision E — SNS review, 26 October 2026

Built the local Revision E SNS chapter against fetched commit `a988508a`. **Critical finding:** 55 of 61 per-target winner fields are discarded by the committed evaluator, and the runner does not forward them. The perturbation figure also duplicates left-heel counts into its right column; several nominal deletions are not applied. These are implementation/evidence issues, not just wording issues. A newer easteregg2 source may differ and needs recovery with its outputs.

Revision E replaces the confusing five-task capability table with gaps, reasons, tests, results, and limits in prose. It incorporates Ben's 26 October per-leg PF / shared four-neuron RG clarification and MN-feedback wording, recent Simulink tests, recorded mass-change outcomes, failed conductance robustness probes, and the negative Di Russo stability search. The defective perturbation figure is explicitly excluded and its source retained. Historical review scores and Ben's replies are preserved. Clarity is still pending Ben/advisor acceptance.

Files: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/Revision_E_local.pdf`, `REPORT.md`, `ZCODE_HANDOFF.md`, and `review/self_scoring_matrix.json`. The report includes exact section/line/text evidence, source hashes, numerical changes, and build results. Overleaf was read only; no commit, push, new modeling run, or external message was made. Other chapters are the Revision D baseline.

---

# Priority 1.1 corrected and expanded with Ben's replies — 23 October 2026

Ben's 2:03, 2:05, and 2:07 p.m. Overleaf replies were captured verbatim in Notes/advisor_review_20261023/ben_conversion_replies_20261023.json. They supply actual conversion problems omitted from the earlier generic prose pass: small-muscle pruning (3 per side, model/count reconciliation pending), OpenSim/ISB axis convention, patella-path representation, and the initial BPA-routed robotbody import with human muscle parameters triggering impossible-force checks. Ben additionally confirms restoration of the OpenSim/ISB frame and a ground-contact issue; the exact fixes, outcomes, and machine provenance require recovery.

The self-contained **ZCODE_HANDOFF_20261023.md at repository root** and Revision_D/ZCODE_HANDOFF.md now include these replies, the added context, exact evidence-recovery questions, and required dissertation/animatlab/mujoco/matlab/neuroscience/sns-toolbox skill use and relevant plugins. Use repository-relative paths and hostname-specific AGENTS.md instructions. GLM-5.3/Flash remain the execution preference; ChatGPT Astra is the final prose preference. These files are local/uncommitted: they must be copied or synchronized by Ben to reach the other machine. No commit, push, message, or new run was performed.

Correction to the earlier handoff: the advisor's examples were hypothetical, but Ben's subsequent issue list was not. It is now explicitly included rather than left as an unspecified request to discover whether any problems existed. Matrix M-126 and M-127 carry the verbatim follow-up evidence. No .tex source was changed in this handoff update, so the prior build is unaffected. Overleaf remains unchanged.

---

# Local prose revision D — 23 October 2026

**18 targeted edits across Introduction, Background, Database, and SNS are complete locally. Overleaf is unchanged.** Ben clarified that ChatGPT Astra is his intended final prose model. No model was launched, no provider was configured, and no scheduled monitoring was restarted.

Working copy: [Revision_D_local.pdf](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/Revision_D_local.pdf>). Actionable code/model handoff: [ZCODE_HANDOFF.md](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/ZCODE_HANDOFF.md>). Exact before/after changes: [changes.diff](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/review/changes.diff>); per-edit record: [changes.json](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/review/changes.json>).

## Source and scope

The five captured Overleaf chapter snapshots supplied the reviewed source. Uncaptured chapters, main.tex, bibliography, and initial assets came from the existing Revision_C copy solely to provide local build support. This is **not a fully reconciled current Overleaf project** and is not an upload bundle. ZCode must merge the targeted changes into fresh source, not replace live chapters wholesale. ProofFinal and Revision_C were not edited by this pass.

The existing flexor panel Flx20mm_meetsHuman.pdf was missing from the Revision_C support assets. It was copied unchanged from tmp/overleaf_live_20261023 for the build. No figure was generated or replaced. The database pivot caption's “animal preparations” was explicitly changed to “experimental preparations”; all other figure blocks were preserved, and all figure inclusions, labels, and citation commands remain conserved.

## Changes and self-check

| Edit | File and final source line | Change / reason | Local self-score |
|---|---|---|---|
| D-01 | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/10-introduction.tex:12>) | Correct subject/verb agreement and tense. | 3/4; acceptance unverified |
| D-02 | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/10-introduction.tex:12>) | Make the dispersed findings the explicit plural subject. | 3/4; acceptance unverified |
| D-03 | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/10-introduction.tex:12>) | Remove duplicate text and misspelling; scope the missing combination to this project instead of asserting an unverified field-wide first. | 3/4; acceptance unverified |
| D-04 | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:10>) | Define archived by its actual app-only status; source sadb_app.html archive control/header. | 3/4; acceptance unverified |
| D-05 | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:10>) | Explain citation direction; verified build_citation_graph.py:8 and build_app.py:320-322,725-744. | 3/4; acceptance unverified |
| D-06 | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:12>) | Explain record types and relationships in reading order; preserve existing example and explicitly avoid claiming complete links. | 2/4; acceptance unverified |
| D-07 | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:25>) | Explain the purpose of controlled vocabulary without inventing category definitions. | 2/4; acceptance unverified |
| D-08 | [50-database.tex:40](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:40>) | Use the existing preparation category without implying that every category is an animal species. | 2/4; acceptance unverified |
| D-09 | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:52>) | Answer the explicit future-work request without claiming all requested links are absent. | 2/4; acceptance unverified |
| D-10 | [60-sns.tex:14](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:14>) | Replace abstract missing-step language with observed outcome, limitation, and specific requirement; the existing 92-actuator value remains later in this chapter. | 2/4; acceptance unverified |
| D-11 | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:27>) | Keep the advisor five-part structure but distinguish tasks from completed systems, remove tested-it ellipsis, and define graphical construction. | 2/4; acceptance unverified |
| D-12 | [60-sns.tex:54](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:54>) | Remove literal ellipsis from table cell; preserve all capabilities and limits. | 3/4; acceptance unverified |
| D-13 | [60-sns.tex:61](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:61>) | Define checks rather than implying unspecified audit discoveries; actual defects and repairs remain a ZCode evidence request. | 2/4; acceptance unverified |
| D-14 | [60-sns.tex:63](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/60-sns.tex:63>) | Define BPA interface inputs, force output, feedback, and limits at first mention without changing model functionality. | 2/4; acceptance unverified |
| D-15 | [15-background.tex:224](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/15-background.tex:224>) | Complete a sentence fragment with the existing scoped attribution. | 3/4; acceptance unverified |
| D-16 | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/15-background.tex:226>) | Remove unsupported overnight throughput/statistical-power assertions; explicitly logged, not silently harmonized. No numerical results or figures removed. | 2/4; acceptance unverified |
| D-17 | [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/15-background.tex:230>) | Clarify repeatability versus biological replication; retain Three Rs definition/citations and the existing figure. | 2/4; acceptance unverified |
| D-18 | [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D/chapters/50-database.tex:8>) | State the project-specific evidence gap in a direct sentence. | 2/4; acceptance unverified |

Scores of 3 concern a directly checked local clarification; scores of 2 retain an unresolved contextual, reader-understanding, or implementation-evidence dependency. No score of 4 was assigned. Earlier 0/4 clarity failures remain in the historical matrix; these local revisions do not erase them. Context IDs in changes.json link nearby matrix entries but do not claim every nearby short comment is addressed.

## Verification

- Fresh pdflatex → bibtex → pdflatex ×2: **257 pages, zero errors, zero undefined references, zero missing figures.** Full logs and build_result.json are in review/ and the named build directory.
- Figure paths/counts, label lists, and citation command lists are unchanged in all four edited files. Bibliography was copied byte-for-byte and never re-encoded.
- All numerical measured/model results are preserved. In SNS, duplicated mentions of the 92-actuator count and Gait2392 name were shortened; their substantive values remain. Literal list ordinals 1–5 became an enumerate list. Numeric-token differences are logged in conservation.json; no result value was silently changed.
- Background's unsupported claims about thousands of overnight combinations, hundreds of cycles per condition, unmatched throughput, and statistical power were explicitly replaced by a qualitative statement of repeatable simulation tests. These were unverified performance claims, not benchmark results. The assertion that identical models are independent biological subjects was explicitly removed. Old wording is retained in changes.json/diff.
- Citation directions checked against SADb_audit/build_citation_graph.py:8 and SADb_audit/app/build_app.py:320–322,725–744. App-only archive status checked against the app header/archive control. These are local code checks, not claims of an exhaustive live-base schema audit.
- Database opening and SNS task-list page renders were inspected: text was legible, with no clipping on those sampled pages. **26 overfull boxes remain in the full build.** The document is not declared layout-clean.

## Rubric disposition and unresolved work

Gate 1: no new external prose copied; existing citation commands and publication notes preserved. Gate 2: changed prose checked for banned boilerplate, unexplained labels, and reviewer-comment leakage; no new reviewer text or agent comments placed in .tex. Gate 3: numerical results conserved, scoped claims checked as above; detailed implementation claims still require ZCode's artifact audit. Gate 4: direct clarifications made; broad clarity concerns remain provisional until Ben/advisor reads the result. **No overall scientific-validation GO is claimed.**

ZCode's handoff specifies evidence to retrieve for actual conversion defects, five-system integration status, bilateral organization/authorship, source-model replication, earlier-library comparison, and database links. It requires existing evidence before new runs, source-backed edge tables before optimization, and one baseline plus one published stimulus/deletion protocol before expanding a replication campaign. GLM-5.3 / Flash are Ben's preferred execution models there. ChatGPT Astra is reserved for bounded final prose polishing after evidence is settled.

All source edits are local and reviewable. No Overleaf changes, uploads, replies, accept/reject actions, resolutions, commits, or pushes were performed.

---

# Advisor review matrix and execution contract

Updated 2026-10-23 12:49:56 Pacific. This is the primary navigation and scoring view. Older captures and accountability findings remain verbatim below. **Overleaf is unchanged, and scheduled monitoring remains canceled.**

Ben's controlling instruction: section, source line, highlighted/changed text, author, and comment must be recorded; advisor edits without comments and Ben's edits count too. This matrix contains 137 current review entries: {'Comment thread': 63, 'Advisor tracked edit': 63, 'Ben tracked edit': 11}. All messages and replies within each thread are retained.

## How to use the references

Line links refer to **frozen read-only copies of the current Overleaf editor source** in Notes/advisor_review_20261023/matrix_source, not stale ProofFinal chapters or PDF pages. Section numbers are calculated from that snapshot's LaTeX headings. Pending text is included as shown by the editor; deleted wording remains in the exact change record. Source hashes are in matrix_manifest.json. Raw character offsets are retained only in JSON for reproducibility.

**Correction to my earlier explanation:** the values were review-card start offsets in JavaScript UTF-16 units, not reliable human navigation references. They give no highlight length. Some old threads collapsed to the end of a replaced chapter. Those rows explicitly say the original current location is unavailable; an end-of-file line is not presented as the original passage. Historical exact highlights are labeled historical. For other comments, the source line is context, not a fabricated exact selection. Exact before/after tracked edits are recorded independently of this limitation.

## Self-scoring rules, used before claiming an item is handled

Ben already explained the unclear prose and gaps. He does not need to restate them more forcefully. The obligation is on the editing agents to demonstrate compliance.

| Score | Required evidence |
|---|---|
| Unassessed | No scored revision/evidence review yet. Not a zero and not a pass. This capture does not award itself quality credit. |
| 0 | Revised passage still fails the specific request, contradicts evidence, or repeats the same documented defect. |
| 1 | Some wording changed, but the referent, gap, causal connection, attribution, or result remains unclear. |
| 2 | Direct answer is present, but context or factual verification is incomplete. |
| 3 | Comment intent, surrounding argument, attribution, and cited evidence checked; the local revision is reviewable. Reader acceptance remains unverified. |
| 4 | Ben or the advisor confirms that the concern is resolved. Agents cannot award themselves this score. |

Every future revision must populate: proposed/actual before and after, evidence path, exact checks, score, remaining limitation, and reviewer disposition. A build checks LaTeX; it does not establish clarity or model validity. A repeated substantive concern reopens the row and caps its prior editorial score at 1 pending reassessment. Do not deduct points from Ben for his edits. His edits and the advisor's edits provide examples for evaluating the agents' prose; authorship in Overleaf is an account label, not proof of which agent wrote it.

**Initial scored failures:** rows with “organization is still unclear”, “This all still feels a bit awkward”, and “This paragraph confuses me” receive 0/4 for the earlier agent clarity revision. This is scored from explicit reader feedback, not inferred from the mere presence of an unresolved thread. Other rows remain unassessed until their actual proposed revision and evidence are checked. This makes no claim of a whole-dissertation quality percentage.

Gap check for each affected section: (1) what was unknown or unavailable, (2) why it mattered, (3) what Ben did, (4) what observed result establishes, and (5) what remains unresolved. A list of software names or missing combinations does not pass. Keep knowledge gaps, engineering capability gaps, completed results, and proposed integration distinct.

## Work allocation and minimum additional evidence

Ben's preference is recorded: this chat specifies minor local edits and review criteria; GLM-5.3 and GLM-5.3-Flash handle additional testing, modeling, and optimization through his ZCode setup; the requested final polishing model is “6.1 Astra”, whose exact available identifier must be confirmed before any launch. No models or other chats were invoked in this pass, and no paid work was launched.

| Work | Smallest useful output | Gate before expansion |
|---|---|---|
| Definitions, pronouns, citation direction, paragraph order | One targeted replacement with before/after and preserved citations. No new simulation. | Reader can identify subject, contribution, and link to the argument. |
| “Results? Lessons?” / conversion “audit” | Existing experiment, specific observation, actual defect/fix if one existed, and limitation, with trace/log evidence. | Search existing records before requesting a rerun. The advisor's hypothetical wrong origins/handedness are not findings. |
| Accurate connectome transcription | Source-backed node and edge tables: source, target, sign, weight/convention, delay, cell model, source page/figure/table, and ambiguity. Overlay/render for visual checking plus machine-readable comparison against the built SNS network. | Resolve ambiguous edges; verify every edge and neuron, including omissions and extra edges. A matching census or a plausible diagram alone cannot pass. Never tune weights to hide uncertain wiring. |
| Reference neural replication | One source model and one published stimulus/deletion protocol first, matching equations, inputs, initial state, units, solver/time step, and scoring window; reference/replica traces and quantitative error. | Predefine observable metrics and tolerances justified by the source. Report transcription fidelity separately from model reduction or adaptation. Rhythm alone is not a complete replication. |
| Coupled neuromechanical replication | Original controller AND source mechanical/actuator/contact assumptions; baseline and the chosen perturbation/deletion result. | A controller port on different mechanics is an adaptation until the difference is characterized. Frozen/NaN runs fail. Support-assisted stepping is not unsupported walking. |
| Optimization | Only after wiring and baseline protocol pass; state which parameters may vary, objective, constraints, fixed controls, seeds, and independent validation. | Do not use optimization as the first response to a transcription or explanation failure. |
| Final polish | Only approved affected paragraphs plus their evidence matrix, one chapter at a time. Preserve numbers, labels, citations, and factual limits. | Fix evidence first; constrain context and output size. No automatic whole-dissertation rewrite or model switching. |

This queue is a proposed execution specification, not a claim that new experiments have been completed. The local replication README still labels its contents drafts; an existing verification script is not execution evidence. CROSS_PLATFORM_REPORT.md contains prior matched tests and explicit limits; those should be examined before commissioning duplicate runs. Replication model status needs a separate artifact audit before being upgraded.

## Matrix index

“Anchor at EOF” means the old location was lost in chapter replacement. Open the card for historical highlights and messages. All 137 rows have an author, type, source context or explicit lost-anchor status, and an unassessed score field.

| ID | Author / type | Section | Source line | Text / action | Self-score |
|---|---|---|---|---|---|
| [M-001](#m-001) | ajh26 / Comment thread | 1.1 Motivation | [10-introduction.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:8>) | This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too. | Unassessed |
| [M-002](#m-002) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:8>) | Changed: repair the discrepancy to fill in the gaps | Unassessed |
| [M-003](#m-003) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:8>) | Changed: rational to accurate | Unassessed |
| [M-004](#m-004) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) | Changed: N to Attempts to implement n | Unassessed |
| [M-005](#m-005) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) | Added: models | Unassessed |
| [M-006](#m-006) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) | Added: , and was difficult to sift through for new researchers to the field | Unassessed |
| [M-007](#m-007) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) | Changed: tested to validated | Unassessed |
| [M-008](#m-008) | ajh26 / Advisor tracked edit | 1.1 Motivation | [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) | Changed: their availability did not establish controller behavior after model conversion, actuator substitution, or the addition of contact feedbac... | Unassessed |
| [M-009](#m-009) | ajh26 / Comment thread | 1.2 The Biomimetic Humanoid Robot Project | [10-introduction.tex:20](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:20>) | double check | Unassessed |
| [M-010](#m-010) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | expand | Unassessed |
| [M-011](#m-011) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | You should have a good idea of how many now that you have gathered them. | Unassessed |
| [M-012](#m-012) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | use mirror language here and for your chapter order. Your chapters should be in the order you present the gaps/you should present the gaps in the o... | Unassessed |
| [M-013](#m-013) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | It is not clear how this work connects to the Biomimetic Humanoid Project and how your work builds on it. | Unassessed |
| [M-014](#m-014) | ajh26 / Comment thread | 2.2.1 Robot Actuation Technologies / Overview | [15-background.tex:36](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:36>) | Motor center is a little off | Unassessed |
| [M-015](#m-015) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | This figure is confusing and probably not needed, | Unassessed |
| [M-016](#m-016) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | not on your collected data though, right? Be more specific. | Unassessed |
| [M-017](#m-017) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | huh? | Unassessed |
| [M-018](#m-018) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | huh? | Unassessed |
| [M-019](#m-019) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | I'm now a little lost...which method does 2392 use? | Unassessed |
| [M-020](#m-020) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | I assume the 23 corresponds to the DoF. What about the 92? | Unassessed |
| [M-021](#m-021) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | non seq | Unassessed |
| [M-022](#m-022) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | I feel like this is jumping all over the place | Unassessed |
| [M-023](#m-023) | ajh26 / Comment thread | 2.5 Neural Control of Locomotion: Central Pattern Generators | [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) | what comparison was made? You don't say | Unassessed |
| [M-024](#m-024) | ajh26 / Comment thread | 2.6 Sensory Feedback in Locomotion | [15-background.tex:165](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:165>) | Chapter 5 attempts to collect and organize all known studies on sensory feedback. | Unassessed |
| [M-025](#m-025) | ajh26 / Comment thread | 2.8 Synthetic Nervous Systems | [15-background.tex:196](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:196>) | non seq | Unassessed |
| [M-026](#m-026) | ajh26 / Comment thread | 2.9 Neuromechanical Simulation | [15-background.tex:204](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:204>) | huh? | Unassessed |
| [M-027](#m-027) | ajh26 / Advisor tracked edit | 2.9 Neuromechanical Simulation | [15-background.tex:210](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:210>) | Changed: ing to s | Unassessed |
| [M-028](#m-028) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) | Added: elaborate | Unassessed |
| [M-029](#m-029) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) | Changed: too elaborate to integrate by hand, joined to to and complex | Unassessed |
| [M-030](#m-030) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) | Changed: too large to wire by hand. to . | Unassessed |
| [M-031](#m-031) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) | Deleted: , and the ones adjacent to it in the same ecosystem, | Unassessed |
| [M-032](#m-032) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:214](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:214>) | Deleted: , and near the ends of a joint's range the routed paths of sparsely wrapped muscles can bunch up or pass through bone, an interference doc... | Unassessed |
| [M-033](#m-033) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:214](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:214>) | Added: | Unassessed |
| [M-034](#m-034) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:216](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:216>) | Added: | Unassessed |
| [M-035](#m-035) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:218](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:218>) | Added: | Unassessed |
| [M-036](#m-036) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:220](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:220>) | Added: | Unassessed |
| [M-037](#m-037) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:222](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:222>) | Changed: Combined with the SNS-Toolbox's MuJoCo backend \citep{nourse_sns_2023}, to SNS Toolbox \citep{nourse_sns_2023} | Unassessed |
| [M-038](#m-038) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: convert the laboratory's to develop a | Unassessed |
| [M-039](#m-039) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: modified to robotic | Unassessed |
| [M-040](#m-040) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Added: , convert it | Unassessed |
| [M-041](#m-041) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Added: , | Unassessed |
| [M-042](#m-042) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: Simulation speed is the practical enabler to This new neuromechanical model can simulate very quickly | Unassessed |
| [M-043](#m-043) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: , to . | Unassessed |
| [M-044](#m-044) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: t to T | Unassessed |
| [M-045](#m-045) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: , to . | Unassessed |
| [M-046](#m-046) | ajh26 / Advisor tracked edit | 2.9.1 Neuromechanical Simulation / Simulation Tools | [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) | Changed: and that to This new model | Unassessed |
| [M-047](#m-047) | bbolen83 / Ben tracked edit | 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case | [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) | Changed: form of reduction to n experimental repeatability | Unassessed |
| [M-048](#m-048) | ajh26 / Comment thread | 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case | [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) | huh? | Unassessed |
| [M-049](#m-049) | bbolen83 / Ben tracked edit | 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case | [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) | Changed: Refinement follows as well: t to T | Unassessed |
| [M-050](#m-050) | ajh26 / Comment thread | 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case | [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) | huh? | Unassessed |
| [M-051](#m-051) | ajh26 / Comment thread | 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case | [15-background.tex:232](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:232>) | Non seq | Unassessed |
| [M-052](#m-052) | bbolen83 / Ben tracked edit | 4.2.1 Results / Flexor on Pin Joint | [30-torque.tex:352](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:352>) | Changed: \begin{sidewaystable}[p] \centering\small\setlength{\tabcolsep}{4pt} \caption[Pre-optimization pinned-flexor goodness of fit]{Pre-optimiza... | Unassessed |
| [M-053](#m-053) | bbolen83 / Ben tracked edit | 4.2.1 Results / Flexor on Pin Joint | [30-torque.tex:384](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:384>) | Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-flexor goodness of fit]{Optimized pinned-flexor goodness of fit at the adopted... | Unassessed |
| [M-054](#m-054) | bbolen83 / Ben tracked edit | 4.2.2 Results / Extensor on Pin Joint | [30-torque.tex:431](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:431>) | Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-extensor goodness of fit]{Optimized pinned-extensor goodness of fit at the ado... | Unassessed |
| [M-055](#m-055) | bbolen83 / Ben tracked edit | 4.2.4 Results / Flexor on Biomimetic Knee | [30-torque.tex:486](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:486>) | Changed: addresses. to ...addresses. Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of Section~\ref{sec:redesign_physic... | Unassessed |
| [M-056](#m-056) | bbolen83 / Ben tracked edit | 4.2.4 Results / Flexor on Biomimetic Knee | [30-torque.tex:489](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:489>) | Changed: \begin{center} \includegraphics[width=1\textwidth]{FlxGrp.eps} \end{center} to \begin{center} \includegraphics[width=0.72\textwidth]{FlxGr... | Unassessed |
| [M-057](#m-057) | bbolen83 / Ben tracked edit | 4.2.4 Results / Flexor on Biomimetic Knee | [30-torque.tex:497](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:497>) | Added: \textbf{(C)} The redesigned two-BPA flexor route | Unassessed |
| [M-058](#m-058) | bbolen83 / Ben tracked edit | 4.3.2 Discussion / Follow-Up Identification and Route Redesign | [30-torque.tex:536](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:536>) | Deleted: + | Unassessed |
| [M-059](#m-059) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) | Changed: The knowledge gap was the lack of a to When I began this work, there was no | Unassessed |
| [M-060](#m-060) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) | Deleted: , available to this project, | Unassessed |
| [M-061](#m-061) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) | Changed: support to be used to design | Unassessed |
| [M-062](#m-062) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) | Added: the | Unassessed |
| [M-063](#m-063) | ajh26 / Comment thread | 5.1 The Resource | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:10>) | This is unclear | Unassessed |
| [M-064](#m-064) | ajh26 / Comment thread | 5.1 The Resource | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:10>) | Unclear why there would be an 'archived' record. | Unassessed |
| [M-065](#m-065) | ajh26 / Comment thread | 5.1 The Resource | [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:10>) | organization is still unclear | 0/4: repeated clarity failure |
| [M-066](#m-066) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Changed: in to is contained in | Unassessed |
| [M-067](#m-067) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Added: , an online database program. It | Unassessed |
| [M-068](#m-068) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Added: tables for | Unassessed |
| [M-069](#m-069) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Added: Pathways | Unassessed |
| [M-070](#m-070) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Added: and | Unassessed |
| [M-071](#m-071) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Deleted: , and Recorders tables | Unassessed |
| [M-072](#m-072) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Changed: The other tables organize review literature and recording methods to Review papers contain a short summary | Unassessed |
| [M-073](#m-073) | ajh26 / Advisor tracked edit | 5.1 The Resource | [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) | Deleted: I curate the linked Airtable base and rebuild the read-only application from its exported records. | Unassessed |
| [M-074](#m-074) | ajh26 / Comment thread | 5.1 The Resource | [50-database.tex:17](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:17>) | I don't really understand this. Also, be clear about what a paper cites, vs what cites the paper | Unassessed |
| [M-075](#m-075) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) | Added: This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. | Unassessed |
| [M-076](#m-076) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) | Changed: Zotero libraries to libraries | Unassessed |
| [M-077](#m-077) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) | Changed: of the AARL and myself, and papers already curated in Airtable to used to provide background and understanding for the previous publicatio... | Unassessed |
| [M-078](#m-078) | ajh26 / Comment thread | 5.2 Curation Methods | [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) | expand | Unassessed |
| [M-079](#m-079) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) | Deleted: and the laboratory's PDF collection | Unassessed |
| [M-080](#m-080) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) | Deleted: This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. | Unassessed |
| [M-081](#m-081) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) | Added: and placed in an archive | Unassessed |
| [M-082](#m-082) | ajh26 / Comment thread | 5.2 Curation Methods | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) | huh? | Unassessed |
| [M-083](#m-083) | ajh26 / Comment thread | 5.2 Curation Methods | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) | huh? | Unassessed |
| [M-084](#m-084) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) | Deleted: ; missing text was flagged rather than supplied from recollection | Unassessed |
| [M-085](#m-085) | ajh26 / Advisor tracked edit | 5.2 Curation Methods | [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) | Deleted: The retained and archived sets distinguish the working collection from the broader audit. | Unassessed |
| [M-086](#m-086) | ajh26 / Advisor tracked edit | 5.3 Coverage and Use | [50-database.tex:35](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:35>) | Added: of the database | Unassessed |
| [M-087](#m-087) | ajh26 / Comment thread | 5.3 Coverage and Use | [50-database.tex:35](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:35>) | animal? | Unassessed |
| [M-088](#m-088) | ajh26 / Advisor tracked edit | 5.3 Coverage and Use | [50-database.tex:43](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:43>) | Added: \section{Discussion} | Unassessed |
| [M-089](#m-089) | ajh26 / Advisor tracked edit | 5.4.1 Discussion / Use in Neural-Controller Design | [50-database.tex:49](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:49>) | Added: \subsection{Further improvement and work} | Unassessed |
| [M-090](#m-090) | ajh26 / Comment thread | 5.4.1 Discussion / Use in Neural-Controller Design | [50-database.tex:50](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:50>) | Expand what all still needs to be done. Some thoughts. Connecting feedback pathways and papers. Connecting reviews with papers and feedback pathways. | Unassessed |
| [M-091](#m-091) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | There are still lots of unknowns | Unassessed |
| [M-092](#m-092) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | After this sentence, define the gap in knowledge and what you did to fill it. | Unassessed |
| [M-093](#m-093) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | I thought it was presented 2x, once without the AI help, and the other with. | Unassessed |
| [M-094](#m-094) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | I thought you migrated it to a different system | Unassessed |
| [M-095](#m-095) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | You don't have to follow this exactly, but you should generally put this in the Chochrane form. https://www.cochrane.org/sites/default/files/upload... | Unassessed |
| [M-096](#m-096) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | what are notes? | Unassessed |
| [M-097](#m-097) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized? | Unassessed |
| [M-098](#m-098) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | approved by who? | Unassessed |
| [M-099](#m-099) | ajh26 / Comment thread | 5.4.2 Discussion / Further improvement and work | [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) | how was this decided? | Unassessed |
| [M-100](#m-100) | ajh26 / Advisor tracked edit | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:8>) | Added: For a course project in ME 610: Neuromechanical Modeling, | Unassessed |
| [M-101](#m-101) | ajh26 / Comment thread | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:8>) | for each leg? | Unassessed |
| [M-102](#m-102) | ajh26 / Comment thread | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:8>) | We? | Unassessed |
| [M-103](#m-103) | ajh26 / Comment thread | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:10>) | We? | Unassessed |
| [M-104](#m-104) | ajh26 / Comment thread | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:10>) | Huh? Just delete all this? | Unassessed |
| [M-105](#m-105) | ajh26 / Comment thread | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:12>) | Results? Lessons? | Unassessed |
| [M-106](#m-106) | ajh26 / Comment thread | 6.1 The AnimatLab Bipedal Walker | [60-sns.tex:14](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:14>) | huh? | 0/4: repeated clarity failure |
| [M-107](#m-107) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: Within this project, five model combinations were missing: a to I developed several software tools to aid in the development of a human-mi... | Unassessed |
| [M-108](#m-108) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: coupled to a CPG in MuJoCo--SNS-Toolbox; its to in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to ... | Unassessed |
| [M-109](#m-109) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Added: to the Gait2392 Model, also | Unassessed |
| [M-110](#m-110) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: ; the to controller. 3) I | Unassessed |
| [M-111](#m-111) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Added: the | Unassessed |
| [M-112](#m-112) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: controller to neural model | Unassessed |
| [M-113](#m-113) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: ; to and tested it... 4) I built a | Unassessed |
| [M-114](#m-114) | ajh26 / Comment thread | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | huh? | Unassessed |
| [M-115](#m-115) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: and to 5) I coupled | Unassessed |
| [M-116](#m-116) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Deleted: coupled | Unassessed |
| [M-117](#m-117) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Changed: These missing combinations prevented a to These tools enable | Unassessed |
| [M-118](#m-118) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Added: s | Unassessed |
| [M-119](#m-119) | ajh26 / Advisor tracked edit | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | Deleted: the same | Unassessed |
| [M-120](#m-120) | ajh26 / Comment thread | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | need a better word | Unassessed |
| [M-121](#m-121) | ajh26 / Comment thread | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) | * | Unassessed |
| [M-122](#m-122) | ajh26 / Comment thread | 6.2 Knowledge Gaps and the New Tools | [60-sns.tex:29](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:29>) | This all still feels a bit awkward. | 0/4: repeated clarity failure |
| [M-123](#m-123) | bbolen83 / Ben tracked edit | 6.3 Software Tools | [60-sns.tex:46](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:46>) | Changed: & Mechanism import achieved & Standalone reflex demonstrations do not establish control of imported models to ...& Mechanism import achiev... | Unassessed |
| [M-124](#m-124) | ajh26 / Comment thread | 6.3.1 Software Tools / Model Conversion and Actuator Coupling | [60-sns.tex:53](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:53>) | I don't really know what this means | Unassessed |
| [M-125](#m-125) | ajh26 / Comment thread | 6.3.1 Software Tools / Model Conversion and Actuator Coupling | [60-sns.tex:53](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:53>) | huh? | Unassessed |
| [M-126](#m-126) | ajh26 / Comment thread | 6.3.1 Software Tools / Model Conversion and Actuator Coupling | [60-sns.tex:53](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:53>) | I feel like this all seems to imply that there were problems. Perhaps if you are more specific here, that would clarify. Example You ran MyoConvert... | Unassessed |
| [M-127](#m-127) | ajh26 / Comment thread | 6.3.1 Software Tools / Model Conversion and Actuator Coupling | [60-sns.tex:55](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:55>) | What BPA interface? Non seq | Unassessed |
| [M-128](#m-128) | bbolen83 / Ben tracked edit | 6.4.3 Verification of the Implemented Interfaces / Closed-Loop Model Coupling | [60-sns.tex:109](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:109>) | Changed: The MuJoCo--Simulink bridge was tested against an independent Python integration using single-muscle activation, different control rates, ... | Unassessed |
| [M-129](#m-129) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | So many details and I don't understand what the goal is. What is the gap in knowledge? What did you do to fill it? | Unassessed |
| [M-130](#m-130) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | huh? | Unassessed |
| [M-131](#m-131) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | why would we need several? | Unassessed |
| [M-132](#m-132) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | This is all super awkward. You need to identify the gap in knowledge and how you fill that gap in knowledge. Much of what you wrote here should be ... | Unassessed |
| [M-133](#m-133) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | I don't understand the connection. | Unassessed |
| [M-134](#m-134) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | Overlapping text in figure | Unassessed |
| [M-135](#m-135) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | huh? | Unassessed |
| [M-136](#m-136) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | Wait, that is all you are saying about the other tools? It is hard for the reader to understand the limitations and gaps without a deeper descripti... | Unassessed |
| [M-137](#m-137) | ajh26 / Comment thread | Original current section unavailable | Anchor at EOF; see card | Did you compare and discuss against the libraries myself and Fletcher made here? | Unassessed |

## Full review cards

<a id="m-001"></a>

### M-001 · 10-introduction.tex · Comment thread

**Author/date:** ajh26 · 14 October, 10:59 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Artificial

**Earlier locator:** INTRO-01, source line 8 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 14 October, 10:59 am:**

~~~text
This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too.
~~~

**Ben (bbolen83) · 23 October, 6:20 am:**

~~~text
I revised the opening to state the sensory-evidence and model-integration gaps directly and explain how the dissertation addresses them. The counts are now in to Chapter 5.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Artificial muscles are a key enabling technology for this class of robot. Among the available actuation technologies, McKibben-style braided pneumatic actuators (BPAs) combine low mass, high force-to-weight ratio, and intrinsic compliance with a force--length curve that qualitatively resembles that of skeletal muscle \citep{liang_comparative_2020, hunt_modeling_2017}. BPAs have been used in quadrupedal robots \citep{aschenbeck_design_2006, hunt_modeling_2017, scharzenberger_design_2019}, musculoskeletal humanoids \citep{asano_design_2017, hitzmann_anthropomorphic_2018}, and dexterous manipulation systems \citep{chen_realizing_2019}. However, BPAs also differ from biological muscle in important ways: their maximum tensile force occurs at resting length rather than at an optimal fiber length, their maximum contraction ratio is substantially lower than that of human muscle, and their force output is a highly nonlinear function of internal pressure, contraction, and loading state \citep{hunt_modeling_2017, sarosi_comparative_2017, liang_comparative_2020}. My experience with \qty{10}{\mm} actuators made this gap concrete: manufacturer specifications overstate the force and contraction available at the short resting lengths a leg design requires, and designs based on those specifications delivered less range of motion and joint torque than intended. Nor could the published empirical models fill in the gaps, because each is fitted to measurements at a single resting length \citep{sarosi_comparative_2017, martens_modeling_2017, hunt_modeling_2017}; none accurately predicts how force changes when a designer chooses a different actuator length. Because muscle resting length is one of the primary degrees of freedom available to the robot designer, the absence of a force model valid at arbitrary resting length blocked accurate actuator selection. Designing a biomimetic robot around BPAs therefore requires quantitative, predictive models of the actuator and of the joint torques it can produce.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-002"></a>

### M-002 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:02 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: repair the discrepancy
to fill in the gaps
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Artificial muscles are a key enabling technology for this class of robot. Among the available actuation technologies, McKibben-style braided pneumatic actuators (BPAs) combine low mass, high force-to-weight ratio, and intrinsic compliance with a force--length curve that qualitatively resembles that of skeletal muscle \citep{liang_comparative_2020, hunt_modeling_2017}. BPAs have been used in quadrupedal robots \citep{aschenbeck_design_2006, hunt_modeling_2017, scharzenberger_design_2019}, musculoskeletal humanoids \citep{asano_design_2017, hitzmann_anthropomorphic_2018}, and dexterous manipulation systems \citep{chen_realizing_2019}. However, BPAs also differ from biological muscle in important ways: their maximum tensile force occurs at resting length rather than at an optimal fiber length, their maximum contraction ratio is substantially lower than that of human muscle, and their force output is a highly nonlinear function of internal pressure, contraction, and loading state \citep{hunt_modeling_2017, sarosi_comparative_2017, liang_comparative_2020}. My experience with \qty{10}{\mm} actuators made this gap concrete: manufacturer specifications overstate the force and contraction available at the short resting lengths a leg design requires, and designs based on those specifications delivered less range of motion and joint torque than intended. Nor could the published empirical models fill in the gaps, because each is fitted to measurements at a single resting length \citep{sarosi_comparative_2017, martens_modeling_2017, hunt_modeling_2017}; none accurately predicts how force changes when a designer chooses a different actuator length. Because muscle resting length is one of the primary degrees of freedom available to the robot designer, the absence of a force model valid at arbitrary resting length blocked accurate actuator selection. Designing a biomimetic robot around BPAs therefore requires quantitative, predictive models of the actuator and of the joint torques it can produce.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-003"></a>

### M-003 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:03 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: rational
to accurate
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Artificial muscles are a key enabling technology for this class of robot. Among the available actuation technologies, McKibben-style braided pneumatic actuators (BPAs) combine low mass, high force-to-weight ratio, and intrinsic compliance with a force--length curve that qualitatively resembles that of skeletal muscle \citep{liang_comparative_2020, hunt_modeling_2017}. BPAs have been used in quadrupedal robots \citep{aschenbeck_design_2006, hunt_modeling_2017, scharzenberger_design_2019}, musculoskeletal humanoids \citep{asano_design_2017, hitzmann_anthropomorphic_2018}, and dexterous manipulation systems \citep{chen_realizing_2019}. However, BPAs also differ from biological muscle in important ways: their maximum tensile force occurs at resting length rather than at an optimal fiber length, their maximum contraction ratio is substantially lower than that of human muscle, and their force output is a highly nonlinear function of internal pressure, contraction, and loading state \citep{hunt_modeling_2017, sarosi_comparative_2017, liang_comparative_2020}. My experience with \qty{10}{\mm} actuators made this gap concrete: manufacturer specifications overstate the force and contraction available at the short resting lengths a leg design requires, and designs based on those specifications delivered less range of motion and joint torque than intended. Nor could the published empirical models fill in the gaps, because each is fitted to measurements at a single resting length \citep{sarosi_comparative_2017, martens_modeling_2017, hunt_modeling_2017}; none accurately predicts how force changes when a designer chooses a different actuator length. Because muscle resting length is one of the primary degrees of freedom available to the robot designer, the absence of a force model valid at arbitrary resting length blocked accurate actuator selection. Designing a biomimetic robot around BPAs therefore requires quantitative, predictive models of the actuator and of the joint torques it can produce.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-004"></a>

### M-004 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:04 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: N
to Attempts to implement n
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Attempts to implement neural control models exposes two further gaps. First, the project lacked a structured, searchable resource connecting measured sensory-afferent properties to the species, preparations, and pathway roles on which controller design depends. Those findings were distributed across primary papers and reviews \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, and was difficult to sift through for new researchers to the field. Second, the project lacked validated neuromechanical models linking a spinal CPG to both a human-scale musculoskeletal model and its BPA-routed robot counterpart. Published CPG architectures and simulation tools supplied components \citep{rybak_modelling_2006, mccrea_organization_2008, markin_afferent_2010, nourse_sns_2023}, but they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses. The SADb and simulation chapters address these evidence and implementation gaps, respectively.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-005"></a>

### M-005 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:04 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: models
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Attempts to implement neural control models exposes two further gaps. First, the project lacked a structured, searchable resource connecting measured sensory-afferent properties to the species, preparations, and pathway roles on which controller design depends. Those findings were distributed across primary papers and reviews \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, and was difficult to sift through for new researchers to the field. Second, the project lacked validated neuromechanical models linking a spinal CPG to both a human-scale musculoskeletal model and its BPA-routed robot counterpart. Published CPG architectures and simulation tools supplied components \citep{rybak_modelling_2006, mccrea_organization_2008, markin_afferent_2010, nourse_sns_2023}, but they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses. The SADb and simulation chapters address these evidence and implementation gaps, respectively.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-006"></a>

### M-006 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:06 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: , and was difficult to sift through for new researchers to the field
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Attempts to implement neural control models exposes two further gaps. First, the project lacked a structured, searchable resource connecting measured sensory-afferent properties to the species, preparations, and pathway roles on which controller design depends. Those findings were distributed across primary papers and reviews \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, and was difficult to sift through for new researchers to the field. Second, the project lacked validated neuromechanical models linking a spinal CPG to both a human-scale musculoskeletal model and its BPA-routed robot counterpart. Published CPG architectures and simulation tools supplied components \citep{rybak_modelling_2006, mccrea_organization_2008, markin_afferent_2010, nourse_sns_2023}, but they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses. The SADb and simulation chapters address these evidence and implementation gaps, respectively.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-007"></a>

### M-007 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:06 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: tested
to validated
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Attempts to implement neural control models exposes two further gaps. First, the project lacked a structured, searchable resource connecting measured sensory-afferent properties to the species, preparations, and pathway roles on which controller design depends. Those findings were distributed across primary papers and reviews \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, and was difficult to sift through for new researchers to the field. Second, the project lacked validated neuromechanical models linking a spinal CPG to both a human-scale musculoskeletal model and its BPA-routed robot counterpart. Published CPG architectures and simulation tools supplied components \citep{rybak_modelling_2006, mccrea_organization_2008, markin_afferent_2010, nourse_sns_2023}, but they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses. The SADb and simulation chapters address these evidence and implementation gaps, respectively.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-008"></a>

### M-008 · 10-introduction.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:07 am  
**Section:** 1.1 Motivation  
**Frozen source line:** [10-introduction.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: their availability did not establish controller behavior after model conversion, actuator substitution, or the addition of contact feedback
to they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Attempts to implement neural control models exposes two further gaps. First, the project lacked a structured, searchable resource connecting measured sensory-afferent properties to the species, preparations, and pathway roles on which controller design depends. Those findings were distributed across primary papers and reviews \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, and was difficult to sift through for new researchers to the field. Second, the project lacked validated neuromechanical models linking a spinal CPG to both a human-scale musculoskeletal model and its BPA-routed robot counterpart. Published CPG architectures and simulation tools supplied components \citep{rybak_modelling_2006, mccrea_organization_2008, markin_afferent_2010, nourse_sns_2023}, but they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses. The SADb and simulation chapters address these evidence and implementation gaps, respectively.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-009"></a>

### M-009 · 10-introduction.tex · Comment thread

**Author/date:** ajh26 · 23 October, 9:10 am  
**Section:** 1.2 The Biomimetic Humanoid Robot Project  
**Frozen source line:** [10-introduction.tex:20](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:20>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 9:10 am:**

~~~text
double check
~~~

**Ben (bbolen83) · 23 October, 10:10 am:**

~~~text
I didn't redesign the linkages
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present work builds directly on a series of studies by me and my collaborators. In my master's project, the kinematics and dynamics of a biomimetic robot leg were matched against the Gait2392 OpenSim model, artificial muscle attachment and wrapping locations were determined for the robot legs, and MATLAB software was developed for evaluating muscle lengths, moment arms, and isometric joint torques over the full range of motion \citep{bolen_bipedal_2019, bolen_determination_2019}. That work showed that simply transplanting human muscle attachment locations onto a BPA-actuated robot does not reproduce human torque profiles; the placement optimization that closes this gap was completed during this dissertation and is summarized in Section~\ref{sec:foundation}. The rest of the platform advanced in parallel, and this dissertation builds on its hardware and demonstrations: the biomimetic four-bar knee linkage with a migrating instantaneous center of rotation, designed and validated by Steele \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}, is the knee of the robot leg whose torque transmission Chapter~\ref{ch:torque} measures, and the canine-inspired quadruped platform Muscle Mutt, formerly called DoggyDeux, demonstrated antagonistic BPA pairs under synthetic nervous system control on a physical robot \citep{scharzenberger_design_2019}, extending an earlier AARL demonstration of central pattern generator control of a pneumatically actuated robot \citep{hunt_development_2017}. The humanoid leg of this dissertation carries that line forward.
~~~

**Work class:** Attribution/scope. **Next check:** Verify contributors and unilateral/bilateral scope from the source record; make the smallest wording correction.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-010"></a>

### M-010 · 10-introduction.tex · Comment thread

**Author/date:** ajh26 · 22 October, 9:40 am  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [10-introduction.tex:47](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:47>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): but it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing.

**Earlier locator:** INTRO-07, source line 12 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 9:40 am:**

~~~text
expand
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I expanded the explanation of the gap and its connection to controller design, rather than retaining the undefined label.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-011"></a>

### M-011 · 10-introduction.tex · Comment thread

**Author/date:** ajh26 · 22 October, 9:41 am  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [10-introduction.tex:47](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:47>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): hundreds if not thousands

**Earlier locator:** INTRO-11, source line 14 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 9:41 am:**

~~~text
You should have a good idea of how many now that you have gathered them.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
Chapter 5 now reports the audited corpus size and retained-paper coverage, with explicit limits on what those counts establish.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-012"></a>

### M-012 · 10-introduction.tex · Comment thread

**Author/date:** ajh26 · 22 October, 9:44 am  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [10-introduction.tex:47](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:47>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): no closed-loop simulation pipeline joining the robot's musculoskeletal model to a biologically realistic neural controller, and no consolidated resource for the afferent data that such a controller requires.

**Earlier locator:** INTRO-13, source line 16 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 9:44 am:**

~~~text
use mirror language here and for your chapter order. Your chapters should be in the order you present the gaps/you should present the gaps in the order your chapters are.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I aligned the gap statements, contribution order, and chapter overview so they follow the same sequence.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-013"></a>

### M-013 · 10-introduction.tex · Comment thread

**Author/date:** ajh26 · 22 October, 9:49 am  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [10-introduction.tex:47](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/10-introduction.tex:47>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): In parallel, a biomimetic four-bar knee linkage with a migrating instantaneous center of rotation was designed and validated by Steele \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}, and a canine-inspired quadruped platform, Muscle Mutt, formerly called DoggyDeux, was developed with antagonistic BPA pairs and synthetic nervous system control \citep{scharzenberger_design_2019}. That platform extended an earlier AARL demonstration of central pattern generator control of a pneumatically actuated robot \citep{hunt_development_2017}.

**Earlier locator:** INTRO-14, source line 22 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 9:49 am:**

~~~text
It is not clear how this work connects to the Biomimetic Humanoid Project and how your work builds on it.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I rewrote this transition to explain how the work extends the lab’s biomimetic humanoid project through actuator characterization, corrected joint-torque predictions, and neural/mechanical model interfaces.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-014"></a>

### M-014 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:42 am  
**Section:** 2.2.1 Robot Actuation Technologies / Overview  
**Frozen source line:** [15-background.tex:36](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:36>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): \begin{figure}[htbp]

**Earlier locator:** BG-03, source line 36 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:42 am:**

~~~text
Motor center is a little off
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I checked the replacement motor schematic: the rotor and stator now share a center.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\begin{figure}[htbp]
~~~

**Work class:** Figure/layout. **Next check:** Check the actual rendered figure and the reader-facing purpose; no model experiment required.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-015"></a>

### M-015 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:45 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): \begin{figure}[htbp]

**Earlier locator:** BG-04, source line 100 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:45 am:**

~~~text
This figure is confusing and probably not needed,
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I removed the confusing graphical abstract and retained the figures that support specific comparisons in the text.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Figure/layout. **Next check:** Check the actual rendered figure and the reader-facing purpose; no model experiment required.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-016"></a>

### M-016 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:50 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): used throughout this dissertation

**Earlier locator:** BG-09, source line 113 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:50 am:**

~~~text
not on your collected data though, right? Be more specific.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
These are OpenSim model-derived moment arms, not measurements from my robot tests. I made that distinction explicit.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-017"></a>

### M-017 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:50 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): stressed

**Earlier locator:** BG-10, source line 115 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:50 am:**

~~~text
huh?
~~~

**Advisor (ajh26) · 22 October, 10:50 am:**

~~~text
definitions are stressed?
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I replaced the unclear terminology with a signed fixed-axis moment-arm definition and a separate explanation of OpenSim’s coordinate-length sensitivities.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-018"></a>

### M-018 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:52 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): written against

**Earlier locator:** BG-13, source line 115 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:52 am:**

~~~text
huh?
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I removed the repository digression and rewrote this explanation around the mechanics and the Gait2392 model.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-019"></a>

### M-019 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:53 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): judged

**Earlier locator:** BG-14, source line 115 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:53 am:**

~~~text
I'm now a little lost...which method does 2392 use?
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I clarified the Gait2392 knee kinematics and distinguished OpenSim’s coordinate-length sensitivity calculation from the fixed-axis geometric moment-arm expression.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-020"></a>

### M-020 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 10:53 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Gait2392

**Earlier locator:** BG-05, source line 111 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 10:53 am:**

~~~text
I assume the 23 corresponds to the DoF. What about the 92?
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I now define Gait2392 as a model with 23 degrees of freedom and 92 musculotendon actuators.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-021"></a>

### M-021 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 11:49 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): igure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed for it

**Earlier locator:** BG-19, source line 131 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 11:49 am:**

~~~text
non seq
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I shortened the geometry explanation and stated what the accompanying figures compare.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-022"></a>

### M-022 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 11:56 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfering with one another, and the result was saved as the gait2392\_robotbody model that the later chapters of this dissertation continue to refine \citep{bolen_bipedal_2019}. 

**Earlier locator:** BG-28, source line 152 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 11:56 am:**

~~~text
I feel like this is jumping all over the place
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I reorganized this passage around the prior model, the human/robot comparison, and the subsequent placement optimization, removing the file-saving chronology.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-023"></a>

### M-023 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 11:57 am  
**Section:** 2.5 Neural Control of Locomotion: Central Pattern Generators  
**Frozen source line:** [15-background.tex:160](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:160>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): The biceps femoris long head shows what the comparison exposed:

**Earlier locator:** BG-30, source line 152 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 22 October, 11:57 am:**

~~~text
what comparison was made? You don't say
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I now state the human-versus-robot comparison and the role of the subsequent attachment-placement optimization explicitly.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The two-level organization is the field's current working architecture, and its models account for specific observations beyond the deletions that motivated the split: the genetically identified interneuron classes mapped onto the model \citep{rybak_organization_2015} and the organization of circuits controlling two-joint muscles \citep{shevtsova_modeling_2016}, a convergence problem a robot muscle routed across two joints also presents. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}. How this architecture is expressed as simulable circuit models is the subject of Section~\ref{sec:bg_compneuro}; the locomotor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-024"></a>

### M-024 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 23 October, 8:39 am  
**Section:** 2.6 Sensory Feedback in Locomotion  
**Frozen source line:** [15-background.tex:165](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:165>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 8:39 am:**

~~~text
Chapter 5 attempts to collect and organize all known studies on sensory feedback.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Sensory feedback is not a passive correction layer but an integral part of locomotor pattern generation. Reviews of afferent control of the CPG \citep{van_de_crommert_review_1998} and of muscle proprioceptive feedback and spinal networks \citep{windhorst_muscle_2007, windhorst_spinal_2022} summarize the classical picture: group Ia afferents from muscle spindles signal muscle length and velocity, group Ib afferents from Golgi tendon organs signal force, and group II afferents signal length, with load-related feedback playing a role distinct from length-related feedback during locomotion \citep{dietz_significance_2000}. Removing proprioceptive feedback degrades locomotor patterns in mice \citep{akay_degradation_2014}; proprioceptive input regulates and adapts locomotor activity \citep{lam_role_2002}; and the relative contributions of proprioceptive and vestibular systems to locomotion are the subject of ongoing study \citep{akay_relative_2021, akay_sensory_2020}. Tonic central and sensory stimuli have been shown to facilitate involuntary air-stepping in humans \citep{selionov_tonic_2009}, demonstrating that simple graded stimulation of motoneuronal and rhythmic circuits can elicit stepping-like output, an experiment readily reproduced in simulation by tonic stimulation of pattern-formation and rhythm-generator populations. The proprioceptive and contact feedback channels designed in Chapter~\ref{ch:sns} follow this afferent organization.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-025"></a>

### M-025 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 22 October, 12:13 pm  
**Section:** 2.8 Synthetic Nervous Systems  
**Frozen source line:** [15-background.tex:196](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:196>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 22 October, 12:13 pm:**

~~~text
non seq
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I removed the misplaced dissertation-controller sentence and kept this section focused on computational neuroscience generally. CPG-specific material is in the CPG discussion.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
For robots, these biological architectures can be implemented as synthetic nervous systems (SNSs): networks of conductance-based spiking and non-spiking neurons that are small enough to run in real time yet biologically grounded enough to be compared against animal data. The approach has been pursued most persistently by the insect-biorobotics lineage of the Quinn laboratory at Case Western Reserve University and by laboratories closely connected to it, my own included. Works in this lineage span several platforms: CPG-based neural control of a pneumatically actuated dog robot \citep{hunt_development_2017}, analytically designed dynamic neural networks for human balance control \citep{hilts_dynamic_2019}, the full-muscle rat hindlimb model and its moment-arm analysis \citep{young_analyzing_2019}, the two-layer CPG neuromechanical walker built on it \citep{deng_neuromechanical_2019}, which supplies this dissertation's controller architecture, and the Muscle Mutt quadruped, a canine-inspired BPA robot designed at my laboratory, the Agile and Adaptive Robotics Laboratory (AARL) at Portland State University, for synthetic nervous system control \citep{scharzenberger_design_2019}. Related controllers from the wider neuromechanical community include sensory-entrained CPG control of bipedal walking \citep{li_bipedal_2017}, inter-leg coordination analyses in stick-insect-inspired controllers \citep{nourse_analyzing_2019}, reinforcement-reflex hybrids in insects \citep{goldsmith_investigating_2021}, and adaptive hindlimb walking controllers \citep{schilling_adaptive_2022}.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-026"></a>

### M-026 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 23 October, 8:42 am  
**Section:** 2.9 Neuromechanical Simulation  
**Frozen source line:** [15-background.tex:204](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:204>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 8:42 am:**

~~~text
huh?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Neuromechanical simulation, the coupling of a musculoskeletal model to a model of the neural controller that drives it, has advanced along three partially overlapping lines. The first line embeds reflex physiology directly into the controller. The muscle-reflex walking model of \citet{geyer_musclereflex_2010} demonstrated that positive force, length, and velocity feedback around each muscle, organized with simple interlimb coordination, is sufficient to produce self-stabilizing human-like walking dynamics and realistic muscle activations without any central pattern generator. That model became the reference point for a generation of reflex-based neuromechanical controllers, including the cat hindlimb models of Section~\ref{sec:bg_compneuro} \citep{markin_afferent_2010, shevtsova_modeling_2016} and their extensions. Hybrid work makes the complementarity of the two strategies explicit: \citet{haeufle_benefit_2018} combined feed-forward muscle activations with neuronal reflex feedback in simulated human walking and found that the robustness benefit of each depends on which muscle carries it.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-027"></a>

### M-027 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:43 am  
**Section:** 2.9 Neuromechanical Simulation  
**Frozen source line:** [15-background.tex:210](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:210>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: ing
to s
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\subsection{Simulation Tools}\label{sec:bg_tools}
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-028"></a>

### M-028 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:44 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: elaborate
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
All three lines depend on a software substrate that is itself part of the method: elaborate musculoskeletal geometry and complex neural models. The tools used in this dissertation are the following.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-029"></a>

### M-029 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:44 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: too elaborate to integrate by hand, joined to
to and complex
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
All three lines depend on a software substrate that is itself part of the method: elaborate musculoskeletal geometry and complex neural models. The tools used in this dissertation are the following.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-030"></a>

### M-030 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:44 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: too large to wire by hand.
to .
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
All three lines depend on a software substrate that is itself part of the method: elaborate musculoskeletal geometry and complex neural models. The tools used in this dissertation are the following.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-031"></a>

### M-031 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:45 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:212](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:212>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: , and the ones adjacent to it in the same ecosystem,
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
All three lines depend on a software substrate that is itself part of the method: elaborate musculoskeletal geometry and complex neural models. The tools used in this dissertation are the following.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-032"></a>

### M-032 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:46 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:214](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:214>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: , and near the ends of a joint's range the routed paths of sparsely wrapped muscles can bunch up or pass through bone, an interference documented within Gait2392-derived routings themselves
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
OpenSim \citep{delp_opensim_2007, seth_opensim_2018}, distributed through the SimTK (Simulation Toolkit) repository, provides the biomechanical modeling side: musculoskeletal models with detailed bone geometry, muscle routing, and force-generating properties, along with tools for analyzing moment arms and isometric torque. These models are approximations with known limits rather than exact anatomies: muscle paths are interpolations over discrete attachment and wrapping data \citep{bolen_bipedal_2019}. Their geometry is accurate for a purpose, namely moment arms and torque predictions over the functional range of motion, and that is the standard against which the robot model of this program is compared. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-033"></a>

### M-033 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:47 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:214](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:214>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: 
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
OpenSim \citep{delp_opensim_2007, seth_opensim_2018}, distributed through the SimTK (Simulation Toolkit) repository, provides the biomechanical modeling side: musculoskeletal models with detailed bone geometry, muscle routing, and force-generating properties, along with tools for analyzing moment arms and isometric torque. These models are approximations with known limits rather than exact anatomies: muscle paths are interpolations over discrete attachment and wrapping data \citep{bolen_bipedal_2019}. Their geometry is accurate for a purpose, namely moment arms and torque predictions over the functional range of motion, and that is the standard against which the robot model of this program is compared. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-034"></a>

### M-034 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:47 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:216](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:216>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: 
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The same SimTK infrastructure hosts SCONE, an open-source package for predictive simulation and optimization of neuromuscular control of biological motion \citep{geijtenbeek_scone_2019}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-035"></a>

### M-035 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:47 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:218](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:218>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: 
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
AnimatLab \citep{cofer_animatlab_2010} provides a graphical environment in which physical bodies and synthetic neural systems are constructed and co-simulated; the AARL's prior walker (Section~\ref{sec:prior_walker}) was built in it. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-036"></a>

### M-036 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:47 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:220](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:220>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: 
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
MuJoCo \citep{todorov_mujoco_2012} provides the physics-engine side: fast, accurate simulation of rigid-body dynamics with contacts, along with muscle-like actuators that include activation dynamics. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-037"></a>

### M-037 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:47 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:222](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:222>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: Combined with the SNS-Toolbox's MuJoCo backend \citep{nourse_sns_2023},
to SNS Toolbox \citep{nourse_sns_2023}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
MyoConverter bridges the biomechanical and physics-engine worlds by converting OpenSim musculoskeletal models to MuJoCo format with optimized muscle kinematics and kinetics \citep{caggiano_myosuite_2022, wang_myosim_2022}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-038"></a>

### M-038 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:48 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: convert the laboratory's
to develop a
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-039"></a>

### M-039 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:49 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: modified
to robotic
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-040"></a>

### M-040 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:49 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: , convert it
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-041"></a>

### M-041 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:49 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: ,
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-042"></a>

### M-042 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:49 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: Simulation speed is the practical enabler
to This new neuromechanical model can simulate very quickly
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-043"></a>

### M-043 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:50 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: ,
to .
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-044"></a>

### M-044 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:50 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: t
to T
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-045"></a>

### M-045 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:50 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: ,
to .
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-046"></a>

### M-046 · 15-background.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 8:50 am  
**Section:** 2.9.1 Neuromechanical Simulation / Simulation Tools  
**Frozen source line:** [15-background.tex:226](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:226>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: and that
to This new model
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system, a pipeline assembled and evaluated in Chapter~\ref{ch:sns}. This new neuromechanical model can simulate very quickly: a single overnight run can sweep thousands of parameter combinations or execute hundreds of simulated gait cycles per condition. Throughput that no human-subject or animal protocol could match. This new model transforms lesion and stimulation studies from bespoke demonstrations into statistically powered experiments.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-047"></a>

### M-047 · 15-background.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:27 am  
**Section:** 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case  
**Frozen source line:** [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: form of reduction
to n experimental repeatability
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The experiments this research program is designed to enable, including deleting a muscle from a limb, blocking a sensory pathway, lesioning part of a rhythm generator, and applying tonic stimulation to spinal circuits, cannot be performed on healthy humans. Invasive intervention on the human neuromuscular system is limited to clinical populations, is tightly constrained by ethics review, and permits neither the within-subject control nor the reversibility that neuromechanical hypotheses demand. The same experiments in animals also carry real moral cost, and their justification is governed by the Three Rs framework: replacement of animal subjects with non-animal methods, reduction of the numbers required, and refinement of the remaining procedures \citep{russell_principles_1959, tannenbaum_russell_2015}. A robot that implements a validated neuromechanical hypothesis is precisely such a replacement method: every deletion, lesion, and stimulation experiment that runs on the synthetic nervous system is an experiment that does not require a living subject, and the robot population can be duplicated exactly, so that ``subjects'' are statistically identical across conditions, an experimental repeatability no breeding colony can offer. The failure modes of these experiments damage hardware, not tissue.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-048"></a>

### M-048 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 23 October, 8:56 am  
**Section:** 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case  
**Frozen source line:** [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 8:56 am:**

~~~text
huh?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The experiments this research program is designed to enable, including deleting a muscle from a limb, blocking a sensory pathway, lesioning part of a rhythm generator, and applying tonic stimulation to spinal circuits, cannot be performed on healthy humans. Invasive intervention on the human neuromuscular system is limited to clinical populations, is tightly constrained by ethics review, and permits neither the within-subject control nor the reversibility that neuromechanical hypotheses demand. The same experiments in animals also carry real moral cost, and their justification is governed by the Three Rs framework: replacement of animal subjects with non-animal methods, reduction of the numbers required, and refinement of the remaining procedures \citep{russell_principles_1959, tannenbaum_russell_2015}. A robot that implements a validated neuromechanical hypothesis is precisely such a replacement method: every deletion, lesion, and stimulation experiment that runs on the synthetic nervous system is an experiment that does not require a living subject, and the robot population can be duplicated exactly, so that ``subjects'' are statistically identical across conditions, an experimental repeatability no breeding colony can offer. The failure modes of these experiments damage hardware, not tissue.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-049"></a>

### M-049 · 15-background.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:28 am  
**Section:** 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case  
**Frozen source line:** [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: Refinement follows as well: t
to T
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The experiments this research program is designed to enable, including deleting a muscle from a limb, blocking a sensory pathway, lesioning part of a rhythm generator, and applying tonic stimulation to spinal circuits, cannot be performed on healthy humans. Invasive intervention on the human neuromuscular system is limited to clinical populations, is tightly constrained by ethics review, and permits neither the within-subject control nor the reversibility that neuromechanical hypotheses demand. The same experiments in animals also carry real moral cost, and their justification is governed by the Three Rs framework: replacement of animal subjects with non-animal methods, reduction of the numbers required, and refinement of the remaining procedures \citep{russell_principles_1959, tannenbaum_russell_2015}. A robot that implements a validated neuromechanical hypothesis is precisely such a replacement method: every deletion, lesion, and stimulation experiment that runs on the synthetic nervous system is an experiment that does not require a living subject, and the robot population can be duplicated exactly, so that ``subjects'' are statistically identical across conditions, an experimental repeatability no breeding colony can offer. The failure modes of these experiments damage hardware, not tissue.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-050"></a>

### M-050 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 23 October, 8:56 am  
**Section:** 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case  
**Frozen source line:** [15-background.tex:230](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:230>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 8:56 am:**

~~~text
huh?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The experiments this research program is designed to enable, including deleting a muscle from a limb, blocking a sensory pathway, lesioning part of a rhythm generator, and applying tonic stimulation to spinal circuits, cannot be performed on healthy humans. Invasive intervention on the human neuromuscular system is limited to clinical populations, is tightly constrained by ethics review, and permits neither the within-subject control nor the reversibility that neuromechanical hypotheses demand. The same experiments in animals also carry real moral cost, and their justification is governed by the Three Rs framework: replacement of animal subjects with non-animal methods, reduction of the numbers required, and refinement of the remaining procedures \citep{russell_principles_1959, tannenbaum_russell_2015}. A robot that implements a validated neuromechanical hypothesis is precisely such a replacement method: every deletion, lesion, and stimulation experiment that runs on the synthetic nervous system is an experiment that does not require a living subject, and the robot population can be duplicated exactly, so that ``subjects'' are statistically identical across conditions, an experimental repeatability no breeding colony can offer. The failure modes of these experiments damage hardware, not tissue.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-051"></a>

### M-051 · 15-background.tex · Comment thread

**Author/date:** ajh26 · 23 October, 8:56 am  
**Section:** 2.10 Testing Neuromechanical Hypotheses on Robots: The Ethical Case  
**Frozen source line:** [15-background.tex:232](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/15-background.tex:232>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 8:56 am:**

~~~text
Non seq
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Figure~\ref{fig:decerebrate_cat} shows what that accounting leaves out. The image is a frame from footage of a real decerebrate-cat experiment, the class of preparation on which the field's understanding of spinal locomotion has been built since its beginnings \citep{brown_intrinsic_1911, decerebrate_cat_footage}: a living cat, its forebrain disconnected by a brainstem transection, walks on the treadmill because the neurons distal to the laceration are still living, still connected, and still firing action potentials. The published record of this work is drawn almost exclusively in schematic, labeled boxes for half-centers and color-coded reflex arcs that are honest about the circuitry and silent about the animal. The figure is included as the counterweight, documenting what the experiments looked like, and therefore what a replacement method stands in for.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-052"></a>

### M-052 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 10:56 am  
**Section:** 4.2.1 Results / Flexor on Pin Joint  
**Frozen source line:** [30-torque.tex:352](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:352>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: \begin{sidewaystable}[p] \centering\small\setlength{\tabcolsep}{4pt} \caption[Pre-optimization pinned-flexor goodness of fit]{Pre-optimization pinned-flexor goodness of fit. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. Scoring-point restrictions are stated in the text.} \label{tab:FlxPinPreGof} \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}lccccccccc@{}} \toprule & \multicolumn{3}{c}{\textbf{Rigid body}} & \multicolumn{3}{c}{\textbf{Improved char.\ alone}} & \multicolumn{3}{c}{\textbf{Back-calculated (hybrid)}} \\ \cmidrule(lr){2-4} \cmidrule(lr){5-7} \cmidrule(l){8-10} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{48.5}{\cm}$ (training) & 5.779 & 1.1761 & 9.798 & 4.871 & 0.8354 & 9.865 & 1.635 & 0.0991 & 3.842 \\ $l_{0}=\qty{41.0}{\cm}$ (training) & 6.633 & 3.2313 & 11.409 & 6.074 & 2.7090 & 9.262 & 2.187 & 0.1184 & 4.243 \\ $l_{0}=\qty{45.7}{\cm}$ (validation) & 5.315 & 0.9408 & 9.991 & 4.465 & 0.6640 & 9.536 & 1.246 & 0.0409 & 2.275 \\ $l_{0}=\qty{40.6}{\cm}$, tendon (validation) & 3.945 & 0.6750 & 7.527 & 3.226 & 0.4512 & 7.376 & 2.085 & 0.0755 & 3.073 \\ \midrule Mean (training) & 6.206 & 2.2037 & 10.604 & 5.472 & 1.7722 & 9.563 & 1.911 & 0.1088 & 4.043 \\ Mean (validation) & 4.630 & 0.8079 & 8.759 & 3.845 & 0.5576 & 8.456 & 1.666 & 0.0582 & 2.674 \\ \bottomrule \end{tabular} \end{sidewaystable}
to \begin{table}[H] \caption{Pre-optimization pinned-knee flexor torque-prediction goodness of fit per test, comparing the three pre-optimization torque predictions of Figure~\ref{fig:FlxPin_preopt}: the Hunt-equation rigid-body prediction, the unoptimized prediction of the improved BPA characterization, and the back-calculated (hybrid) torque. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:FlxPinPreGof} \begingroup \renewcommand{\arraystretch}{1}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc ccccccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Hunt equation}} & \multicolumn{3}{c}{\textbf{Improved char.\ alone}} & \multicolumn{3}{c}{\textbf{Back-calculated (hybrid)}} \\ \cmidrule(lr){3-5} \cmidrule(lr){6-8} \cmidrule(l){9-11} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule 48.5 & -- & 5.779 & 1.1761 & 9.798 & 4.871 & 0.8354 & 9.865 & 1.635 & 0.0991 & 3.842 \\ 41.0 & -- & 6.633 & 3.2313 & 11.409 & 6.074 & 2.7090 & 9.262 & 2.187 & 0.1184 & 4.243 \\ 45.7 & -- & 5.315 & 0.9408 & 9.991 & 4.465 & 0.6640 & 9.536 & 1.246 & 0.0409 & 2.275 \\ 40.6 & 35 & 3.945 & 0.6750 & 7.527 & 3.226 & 0.4512 & 7.376 & 2.085 & 0.0755 & 3.073 \\ \bottomrule \end{tabular} \endgroup \end{table}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\begin{table}[H]
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-053"></a>

### M-053 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 10:59 am  
**Section:** 4.2.1 Results / Flexor on Pin Joint  
**Frozen source line:** [30-torque.tex:384](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:384>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-flexor goodness of fit]{Optimized pinned-flexor goodness of fit at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.} \label{tab:FlxPinGof} \begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}} \toprule & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){2-4} \cmidrule(l){5-7} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{48.5}{\cm}$ (training) & 7.617 & 2.1278 & 25.531 & 2.165 & 0.1718 & 5.467 \\ $l_{0}=\qty{41.0}{\cm}$ (training) & 9.694 & 2.1149 & 21.076 & 1.084 & 0.0277 & 1.766 \\ $l_{0}=\qty{45.7}{\cm}$ (validation) & 8.803 & 1.9380 & 26.264 & 1.368 & 0.0468 & 3.614 \\ $l_{0}=\qty{40.6}{\cm}$, tendon (validation) & 6.373 & 0.5935 & 16.859 & 1.611 & 0.0427 & 5.891 \\ \midrule Mean (training) & 8.656 & 2.1213 & 23.303 & 1.624 & 0.0998 & 3.616 \\ Mean (validation) & 7.588 & 1.2657 & 21.562 & 1.490 & 0.0447 & 4.752 \\ \bottomrule \end{tabularx} \end{table}
to \begin{table}[H] \caption{Pinned-knee flexor torque-prediction goodness of fit: the improved BPA characterization alone (the zero-correction prediction of the text) and the corrected model at the adopted correction terms (Table~\ref{tab:AdoptedXi}); training and validation tests follow the adopted cross-validation fold. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:FlxPinGof} \begingroup \renewcommand{\arraystretch}{1.4}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc cccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){3-5} \cmidrule(l){6-8} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule \multicolumn{8}{@{}l}{\textbf{Training}} \\ 48.5 & -- & 7.617 & 2.1278 & 25.531 & 2.165 & 0.1718 & 5.467 \\ 41.0 & -- & 9.694 & 2.1149 & 21.076 & 1.084 & 0.0277 & 1.766 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 8.656 & 2.1213 & 23.303 & 1.624 & 0.0998 & 3.617 \\ \midrule \multicolumn{8}{@{}l}{\textbf{Validation}} \\ 45.7 & -- & 8.803 & 1.9380 & 26.264 & 1.368 & 0.0468 & 3.614 \\ 40.6 & 35 & 6.373 & 0.5935 & 16.859 & 1.611 & 0.0427 & 5.891 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 7.588 & 1.2657 & 21.562 & 1.490 & 0.0447 & 4.752 \\ \bottomrule \end{tabular} \endgroup \end{table}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\begin{table}[H]
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-054"></a>

### M-054 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 11:01 am  
**Section:** 4.2.2 Results / Extensor on Pin Joint  
**Frozen source line:** [30-torque.tex:431](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:431>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-extensor goodness of fit]{Optimized pinned-extensor goodness of fit at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.} \label{tab:ExtPinGof} \begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}} \toprule & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){2-4} \cmidrule(l){5-7} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{43.6}{\cm}$ (training) & 2.437 & 0.9618 & 4.284 & 0.978 & 0.1634 & 1.818 \\ $l_{0}=\qty{43.2}{\cm}$, tendon (training) & 2.878 & 4.7304 & 4.817 & 1.525 & 1.3270 & 2.481 \\ $l_{0}=\qty{45.7}{\cm}$ (training) & 2.178 & 1.0586 & 4.974 & 0.490 & 0.0535 & 1.206 \\ $l_{0}=\qty{40.0}{\cm}$ (validation) & 3.092 & 0.5638 & 4.842 & 1.244 & 0.1303 & 2.196 \\ $l_{0}=\qty{40.0}{\cm}$, tendon (validation) & 3.622 & 6.3919 & 7.845 & 1.350 & 0.8878 & 1.890 \\ $l_{0}=\qty{46.5}{\cm}$ (validation) & 0.949 & 0.3152 & 2.272 & 0.699 & 0.1714 & 1.559 \\ \midrule Mean (training) & 2.498 & 2.2502 & 4.692 & 0.997 & 0.5146 & 1.835 \\ Mean (validation) & 2.554 & 2.4236 & 4.986 & 1.098 & 0.3965 & 1.882 \\ \bottomrule \end{tabularx} \end{table}
to \begin{table}[H] \caption{Pinned-knee extensor torque-prediction goodness of fit, at the adopted correction terms with the effective-stiffness pair locked to the flexor-identified values (Table~\ref{tab:AdoptedXi}): the improved BPA characterization alone (the zero-correction prediction of the text) and the corrected model. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:ExtPinGof} \begingroup \renewcommand{\arraystretch}{1}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc cccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){3-5} \cmidrule(l){6-8} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule \multicolumn{8}{@{}l}{\textbf{Training}} \\ 43.6 & -- & 2.437 & 0.9618 & 4.284 & 0.978 & 0.1634 & 1.818 \\ 43.2 & 22 & 2.878 & 4.7304 & 4.817 & 1.525 & 1.3270 & 2.481 \\ 45.7 & -- & 2.178 & 1.0586 & 4.974 & 0.490 & 0.0535 & 1.206 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 2.498 & 2.2503 & 4.692 & 0.997 & 0.5146 & 1.835 \\ \midrule \multicolumn{8}{@{}l}{\textbf{Validation}} \\ 40.0 & -- & 3.092 & 0.5638 & 4.842 & 1.244 & 0.1303 & 2.196 \\ 40.0 & 40 & 3.622 & 6.3919 & 7.845 & 1.350 & 0.8878 & 1.890 \\ 46.5 & -- & 0.949 & 0.3152 & 2.272 & 0.699 & 0.1714 & 1.559 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 2.554 & 2.4236 & 4.986 & 1.098 & 0.3965 & 1.882 \\ \bottomrule \end{tabular} \endgroup \end{table}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\begin{table}[H]
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-055"></a>

### M-055 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:48 am  
**Section:** 4.2.4 Results / Flexor on Biomimetic Knee  
**Frozen source line:** [30-torque.tex:486](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:486>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: addresses.
to ...addresses. Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of Section~\ref{sec:redesign_physical}: across its fifteen tested measurements from \ang{-119} to \ang{5}, the measured torque magnitude met or exceeded the human biceps-femoris short-head maximum at every point, with margins from \qty{0.42}{\newton\meter} at \ang{-105} to \qty{17.0}{\newton\meter} at \ang{-14}.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The torque of the single $\diameter$\qty{20}{\mm} flexor actuator of panel \textbf{B} did not exceed the magnitude of the maximum human isometric flexion torque between \ang{60} and \ang{100} of knee flexion, and within a $\pm 20\%$ band about the human mean it still fell short of human-level torque between \ang{78} and \ang{103} of flexion. Because operating at \qty{620}{\kPa} across the full desired range of motion would load the printed femur beyond its ultimate strength, single points were also recorded at the reduced pressures \qtylist{560;421;325;281}{\kPa}, which show that the desired range of motion is achievable at reduced pressure. This shortfall in flexion torque defines the target that the attachment-route redesign of Section~\ref{sec:ongoing} ...addresses. Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of Section~\ref{sec:redesign_physical}: across its fifteen tested measurements from \ang{-119} to \ang{5}, the measured torque magnitude met or exceeded the human biceps-femoris short-head maximum at every point, with margins from \qty{0.42}{\newton\meter} at \ang{-105} to \qty{17.0}{\newton\meter} at \ang{-14}.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-056"></a>

### M-056 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:50 am  
**Section:** 4.2.4 Results / Flexor on Biomimetic Knee  
**Frozen source line:** [30-torque.tex:489](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:489>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: \begin{center} \includegraphics[width=1\textwidth]{FlxGrp.eps} \end{center}
to \begin{center} \includegraphics[width=0.72\textwidth]{FlxGrp.eps} \par\smallskip \begin{minipage}[t]{0.30\textwidth}\vspace{0pt}\centering {\textbf{(C)}}\par\smallskip \includegraphics[width=\linewidth]{Flx20mm_meetsHuman.pdf} \end{minipage} \end{center}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\begin{center}
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-057"></a>

### M-057 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:53 am  
**Section:** 4.2.4 Results / Flexor on Biomimetic Knee  
**Frozen source line:** [30-torque.tex:497](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:497>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added:  \textbf{(C)} The redesigned two-BPA flexor route
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\caption{Biomimetic-knee torque versus knee angle for the flexor BPAs. Knee angle is negative in flexion, and flexion torque is negative. \textbf{(A)} $\diameter$\qty{10}{\mm}, $l_{0}=\qty{41.5}{\cm}$ actuator with a \qty{21}{\mm} artificial tendon at \qty{604}{\kPa}. \textbf{(B)} $\diameter$\qty{20}{\mm}, $l_{0}=\qty{41.5}{\cm}$ actuator with a \qty{14}{\mm} artificial tendon at \qty{620}{\kPa}. Indigo filled circles are measured torque; the magenta line is the optimized prediction of the corrected model using the flexor terms identified on the pinned knee; the dashed light orange line is the unoptimized prediction of the improved BPA characterization alone. In \textbf{(B)}, the black dotted line is the maximum human isometric torque, and the pink circles are measured points at the reduced pressures \qtylist{560;421;325;281}{\kPa} (from left to right, respectively), recorded because full pressure would load the printed femur beyond its ultimate strength. \textbf{(C)} The redesigned two-BPA flexor route}\label{fig:FlxGrp}
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-058"></a>

### M-058 · 30-torque.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 10:53 am  
**Section:** 4.3.2 Discussion / Follow-Up Identification and Route Redesign  
**Frozen source line:** [30-torque.tex:536](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/30-torque.tex:536>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: +
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Flexor & $0.00394$ & $3.998\times10^{4}$ & $1.473\times10^{4}$ & $0.158$\\
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-059"></a>

### M-059 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:48 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: The knowledge gap was the lack of a
to When I began this work, there was no
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
When I began this work, there was no structured database of which measured afferent properties and feedback pathways could be used to design the robot's controller. Relevant findings were dispersed across species, preparations, tasks, and primary papers. Reviews organized parts of this evidence \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, but selecting a feedback rule still required recovering its experimental conditions and separating measurements from modeling assumptions. I built the SADb to connect those findings to their sources, afferent classes, and pathway roles so that controller choices could be traced to the evidence and remaining unknowns could be identified.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-060"></a>

### M-060 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:48 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: , available to this project,
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
When I began this work, there was no structured database of which measured afferent properties and feedback pathways could be used to design the robot's controller. Relevant findings were dispersed across species, preparations, tasks, and primary papers. Reviews organized parts of this evidence \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, but selecting a feedback rule still required recovering its experimental conditions and separating measurements from modeling assumptions. I built the SADb to connect those findings to their sources, afferent classes, and pathway roles so that controller choices could be traced to the evidence and remaining unknowns could be identified.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-061"></a>

### M-061 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:49 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: support
to be used to design
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
When I began this work, there was no structured database of which measured afferent properties and feedback pathways could be used to design the robot's controller. Relevant findings were dispersed across species, preparations, tasks, and primary papers. Reviews organized parts of this evidence \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, but selecting a feedback rule still required recovering its experimental conditions and separating measurements from modeling assumptions. I built the SADb to connect those findings to their sources, afferent classes, and pathway roles so that controller choices could be traced to the evidence and remaining unknowns could be identified.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-062"></a>

### M-062 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:49 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added:  the
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
When I began this work, there was no structured database of which measured afferent properties and feedback pathways could be used to design the robot's controller. Relevant findings were dispersed across species, preparations, tasks, and primary papers. Reviews organized parts of this evidence \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, but selecting a feedback rule still required recovering its experimental conditions and separating measurements from modeling assumptions. I built the SADb to connect those findings to their sources, afferent classes, and pathway roles so that controller choices could be traced to the evidence and remaining unknowns could be identified.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-063"></a>

### M-063 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 9:51 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:10>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 9:51 am:**

~~~text
This is unclear
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The primary deliverable is a browser application for searching paper records, counting them by selected classifications, and exploring citation neighborhoods (Figure~\ref{fig:sadb_clusters}). Its embedded data support offline use; its separate online-search tab requires a network connection. The distributed snapshot contains 1,393 records: 943 papers and 450 linked model, review, feedback-pathway, and recorder records. Archived records remain available alongside the retained collection. The pivot and citation views operate on papers; the table view searches across record types.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-064"></a>

### M-064 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 9:52 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:10>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 9:52 am:**

~~~text
Unclear why there would be an 'archived' record.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The primary deliverable is a browser application for searching paper records, counting them by selected classifications, and exploring citation neighborhoods (Figure~\ref{fig:sadb_clusters}). Its embedded data support offline use; its separate online-search tab requires a network connection. The distributed snapshot contains 1,393 records: 943 papers and 450 linked model, review, feedback-pathway, and recorder records. Archived records remain available alongside the retained collection. The pivot and citation views operate on papers; the table view searches across record types.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-065"></a>

### M-065 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 9:52 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:10>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 9:52 am:**

~~~text
organization is still unclear
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The primary deliverable is a browser application for searching paper records, counting them by selected classifications, and exploring citation neighborhoods (Figure~\ref{fig:sadb_clusters}). Its embedded data support offline use; its separate online-search tab requires a network connection. The distributed snapshot contains 1,393 records: 943 papers and 450 linked model, review, feedback-pathway, and recorder records. Archived records remain available alongside the retained collection. The pivot and citation views operate on papers; the table view searches across record types.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** 0/4 on the prior agent clarity revision: explicit repeated reader criticism; not a claim that every sentence is factually wrong. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-066"></a>

### M-066 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:54 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: in
to is contained in
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-067"></a>

### M-067 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:54 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: , an online database program. It
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-068"></a>

### M-068 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:57 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: tables for
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-069"></a>

### M-069 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:56 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added:  Pathways
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-070"></a>

### M-070 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:56 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: and
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-071"></a>

### M-071 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:57 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: , and Recorders tables
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-072"></a>

### M-072 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 9:57 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: The other tables organize review literature and recording methods
to Review papers contain a short summary
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-073"></a>

### M-073 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:00 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted:  I curate the linked Airtable base and rebuild the read-only application from its exported records.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The working base is contained in Airtable, an online database program. It contains linked tables for Papers, Feedback Pathways, Models, and Review Papers. Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text. Feedback records describe a sensory pathway, such as Ib feedback influencing the stance-to-swing transition, and link it to supporting papers; Models link computational implementations to their source literature. Review papers contain a short summary. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-074"></a>

### M-074 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:26 am  
**Section:** 5.1 The Resource  
**Frozen source line:** [50-database.tex:17](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:17>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:26 am:**

~~~text
I don't really understand this. Also, be clear about what a paper cites, vs what cites the paper
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\caption[SADb citation-neighborhood view]{Paper-only citation view focused on Grillner (1985). Soma size represents citation count; axes show year and citation count. Triangles mark 57 in-corpus citers, and circles 13 cited papers. Second-degree papers are gray. The panels identify the selected paper and its connections.}
~~~

**Work class:** Clarification; scope to verify. **Next check:** Define both directions: references are works this paper cites; citing papers are works that cite it. Then verify which relationship the app displays. No new model experiment.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-075"></a>

### M-075 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:01 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. The collection began with the libraries used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab. Records were reconciled by DOI and title, then expanded through citation searches. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-076"></a>

### M-076 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:19 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: Zotero libraries
to libraries
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. The collection began with the libraries used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab. Records were reconciled by DOI and title, then expanded through citation searches. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-077"></a>

### M-077 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:19 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: of the AARL and myself, and papers already curated in Airtable
to used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. The collection began with the libraries used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab. Records were reconciled by DOI and title, then expanded through citation searches. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-078"></a>

### M-078 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:21 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:21 am:**

~~~text
expand
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. The collection began with the libraries used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab. Records were reconciled by DOI and title, then expanded through citation searches. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. 
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-079"></a>

### M-079 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:21 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted:  and the laboratory's PDF collection
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. The collection began with the libraries used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab. Records were reconciled by DOI and title, then expanded through citation searches. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-080"></a>

### M-080 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:01 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:23](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:23>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication. The collection began with the libraries used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab. Records were reconciled by DOI and title, then expanded through citation searches. Reviews of spinal and proprioceptive feedback supplied an initial organization of the literature \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-081"></a>

### M-081 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:24 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added:  and placed in an archive
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Records were retained when they reported sensory-afferent measurements, investigated feedback in locomotion or posture, or modeled a relevant sensory pathway. Reviews were retained as sources for organizing and locating that evidence. Duplicate records and material outside this feedback scope were removed from the working base and placed in an archive. Classification used controlled vocabularies for preparations and afferent classes. Each mechanism summary was grounded in available source text, using the abstract when full text was unavailable. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-082"></a>

### M-082 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:22 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:22 am:**

~~~text
huh?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Records were retained when they reported sensory-afferent measurements, investigated feedback in locomotion or posture, or modeled a relevant sensory pathway. Reviews were retained as sources for organizing and locating that evidence. Duplicate records and material outside this feedback scope were removed from the working base and placed in an archive. Classification used controlled vocabularies for preparations and afferent classes. Each mechanism summary was grounded in available source text, using the abstract when full text was unavailable. 
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-083"></a>

### M-083 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:22 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:22 am:**

~~~text
huh?
~~~

**Advisor (ajh26) · 23 October, 10:22 am:**

~~~text
sensory feedback pathway?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Records were retained when they reported sensory-afferent measurements, investigated feedback in locomotion or posture, or modeled a relevant sensory pathway. Reviews were retained as sources for organizing and locating that evidence. Duplicate records and material outside this feedback scope were removed from the working base and placed in an archive. Classification used controlled vocabularies for preparations and afferent classes. Each mechanism summary was grounded in available source text, using the abstract when full text was unavailable. 
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-084"></a>

### M-084 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:23 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: ; missing text was flagged rather than supplied from recollection
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Records were retained when they reported sensory-afferent measurements, investigated feedback in locomotion or posture, or modeled a relevant sensory pathway. Reviews were retained as sources for organizing and locating that evidence. Duplicate records and material outside this feedback scope were removed from the working base and placed in an archive. Classification used controlled vocabularies for preparations and afferent classes. Each mechanism summary was grounded in available source text, using the abstract when full text was unavailable. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-085"></a>

### M-085 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:24 am  
**Section:** 5.2 Curation Methods  
**Frozen source line:** [50-database.tex:25](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:25>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted: The retained and archived sets distinguish the working collection from the broader audit.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Records were retained when they reported sensory-afferent measurements, investigated feedback in locomotion or posture, or modeled a relevant sensory pathway. Reviews were retained as sources for organizing and locating that evidence. Duplicate records and material outside this feedback scope were removed from the working base and placed in an archive. Classification used controlled vocabularies for preparations and afferent classes. Each mechanism summary was grounded in available source text, using the abstract when full text was unavailable. 
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-086"></a>

### M-086 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:29 am  
**Section:** 5.3 Coverage and Use  
**Frozen source line:** [50-database.tex:35](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:35>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: of the database
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
A user of the database can begin with a preparation or afferent class, retrieve the matching papers, and inspect their notes and source links. The pivot view counts papers along two selected dimensions (Figure~\ref{fig:sadb_pivot}); the citation view provides a second route from an anchor paper to related literature (Figure~\ref{fig:sadb_clusters}). These views help locate evidence, while the original experiment determines the conditions under which a finding applies.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-087"></a>

### M-087 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:29 am  
**Section:** 5.3 Coverage and Use  
**Frozen source line:** [50-database.tex:35](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:35>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:29 am:**

~~~text
animal?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
A user of the database can begin with a preparation or afferent class, retrieve the matching papers, and inspect their notes and source links. The pivot view counts papers along two selected dimensions (Figure~\ref{fig:sadb_pivot}); the citation view provides a second route from an anchor paper to related literature (Figure~\ref{fig:sadb_clusters}). These views help locate evidence, while the original experiment determines the conditions under which a finding applies.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-088"></a>

### M-088 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:31 am  
**Section:** 5.3 Coverage and Use  
**Frozen source line:** [50-database.tex:43](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:43>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added:  \section{Discussion}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex

~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-089"></a>

### M-089 · 50-database.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:33 am  
**Section:** 5.4.1 Discussion / Use in Neural-Controller Design  
**Frozen source line:** [50-database.tex:49](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:49>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added:  \subsection{Further improvement and work}
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex

~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-090"></a>

### M-090 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:34 am  
**Section:** 5.4.1 Discussion / Use in Neural-Controller Design  
**Frozen source line:** [50-database.tex:50](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:50>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:34 am:**

~~~text
Expand what all still needs to be done. Some thoughts. Connecting feedback pathways and papers. Connecting reviews with papers and feedback pathways.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
\subsection{Further improvement and work}
~~~

**Work class:** Future-work scope. **Next check:** Check existing schema links; specify the missing paper/pathway/review links without presenting finished links as future work.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-091"></a>

### M-091 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:14 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Every element of that specification was established experimentally

**Earlier locator:** DB-04, source line 6 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:14 am:**

~~~text
There are still lots of unknowns
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I removed the suggestion that the physiology is fully specified. The chapter now states the uneven coverage and the interpretation still needed before using a pathway in a controller.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-092"></a>

### M-092 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:16 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): literature

**Earlier locator:** DB-05, source line 6 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:16 am:**

~~~text
After this sentence, define the gap in knowledge and what you did to fill it.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I put the knowledge gap immediately after introducing SADb: the project lacked a structured, source-linked account of afferent measurements and pathways for selecting feedback rules. The next sentences explain how SADb addresses it.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-093"></a>

### M-093 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:20 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Neuroscience 2023 \citep{bolen_sensory_2023}

**Earlier locator:** DB-09, source line 10 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:20 am:**

~~~text
I thought it was presented 2x, once without the AI help, and the other with.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
Both presentations are now identified and cited separately: the database in 2023 and AI-assisted triage in 2024.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-094"></a>

### M-094 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:25 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): The database lives in a five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders)

**Earlier locator:** DB-13, source line 10 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:25 am:**

~~~text
I thought you migrated it to a different system
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I clarified the two roles: Airtable remains the editable linked working base, and the offline browser application is its read-only distribution. The app now carries all 1,393 pre-prune records; pivot and citation views remain paper-only.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-095"></a>

### M-095 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:33 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Methods

**Earlier locator:** DB-18, source line 14 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:33 am:**

~~~text
You don't have to follow this exactly, but you should generally put this in the Chochrane form. https://www.cochrane.org/sites/default/files/uploads/PDFs/guide_to_the_contents_of_a_cochrane_methodology_protocol_and_review.pdf
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I reorganized Curation Methods around sources, identity matching, inclusion criteria, source-text checking, and my decision authority. I also state that this was a curated research collection rather than an exhaustive or registered systematic review.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-096"></a>

### M-096 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:37 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): note

**Earlier locator:** DB-21, source line 18 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:37 am:**

~~~text
what are notes?
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I now define curation notes as short, source-grounded summaries of the mechanism, experimental conditions, and relevance to feedback design.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-097"></a>

### M-097 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:39 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): A note is written from the paper's own text, from the abstract at minimum and from the full text when available, never from the curator's recollection of the paper;

**Earlier locator:** DB-22, source line 18 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:39 am:**

~~~text
You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized?
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
The chapter now opens with the finished resource, its purpose, record types, organization, and a screenshot. Implementation history and supporting counts no longer precede that overview.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-098"></a>

### M-098 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:40 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): approved

**Earlier locator:** DB-23, source line 22 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:40 am:**

~~~text
approved by who?
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I made my responsibility explicit: I decide inclusion, classification, and deletion after checking the source.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-099"></a>

### M-099 · 50-database.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:40 am  
**Section:** 5.4.2 Discussion / Further improvement and work  
**Frozen source line:** [50-database.tex:52](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/50-database.tex:52>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): left as off-topic or orphaned

**Earlier locator:** DB-24, source line 22 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:40 am:**

~~~text
how was this decided?
~~~

**Advisor (ajh26) · 21 October, 11:41 am:**

~~~text
how were they identified?
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I added the relevance criteria and described DOI/title matching and source-text verification. Missing source text is flagged rather than reconstructed from recollection.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The present corpus is uneven. Its afferent classifications cover only part of the retained literature, and the depth of the mammalian locomotor material should not be taken as equal coverage of vestibular, visual, or other sensory systems. The application makes the curated papers easier to inspect and compare; additional curation and experimental interpretation remain necessary for extending the robot's sensory system.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-100"></a>

### M-100 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 10:48 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: For a course project in ME 610: Neuromechanical Modeling,
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
For a course project in ME 610: Neuromechanical Modeling, Morrow and I combined the bipedal mechanics of \citet{li_bipedal_2017} with the two-level central pattern generator (CPG) of \citet{deng_neuromechanical_2019} in AnimatLab \citep{cofer_animatlab_2010}. The mechanical model had a pelvis and bilateral thigh, shank, foot, and toe segments, with sagittal-plane joints and antagonist Hill-type muscles at the hip, knee, and ankle. I added a passive toe joint. The initial controller used one rhythm generator (RG) to control both sides, with a separate pattern-formation layer. I added afferent feedback and heel- and toe-contact signals during the subsequent development. Figure~\ref{fig:animatlab_model_architecture} shows the model and controller organization.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-101"></a>

### M-101 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:49 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:49 am:**

~~~text
for each leg?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
For a course project in ME 610: Neuromechanical Modeling, Morrow and I combined the bipedal mechanics of \citet{li_bipedal_2017} with the two-level central pattern generator (CPG) of \citet{deng_neuromechanical_2019} in AnimatLab \citep{cofer_animatlab_2010}. The mechanical model had a pelvis and bilateral thigh, shank, foot, and toe segments, with sagittal-plane joints and antagonist Hill-type muscles at the hip, knee, and ankle. I added a passive toe joint. The initial controller used one rhythm generator (RG) to control both sides, with a separate pattern-formation layer. I added afferent feedback and heel- and toe-contact signals during the subsequent development. Figure~\ref{fig:animatlab_model_architecture} shows the model and controller organization.
~~~

**Work class:** Attribution/scope. **Next check:** Verify contributors and unilateral/bilateral scope from the source record; make the smallest wording correction.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-102"></a>

### M-102 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:49 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:8](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:8>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:49 am:**

~~~text
We?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
For a course project in ME 610: Neuromechanical Modeling, Morrow and I combined the bipedal mechanics of \citet{li_bipedal_2017} with the two-level central pattern generator (CPG) of \citet{deng_neuromechanical_2019} in AnimatLab \citep{cofer_animatlab_2010}. The mechanical model had a pelvis and bilateral thigh, shank, foot, and toe segments, with sagittal-plane joints and antagonist Hill-type muscles at the hip, knee, and ankle. I added a passive toe joint. The initial controller used one rhythm generator (RG) to control both sides, with a separate pattern-formation layer. I added afferent feedback and heel- and toe-contact signals during the subsequent development. Figure~\ref{fig:animatlab_model_architecture} shows the model and controller organization.
~~~

**Work class:** Attribution/scope. **Next check:** Verify contributors and unilateral/bilateral scope from the source record; make the smallest wording correction.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-103"></a>

### M-103 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:50 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:10>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:50 am:**

~~~text
We?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I adjusted applied currents, synaptic weights, and muscle parameters by hand. The suspended model developed alternating muscle activity and stepping-like motion, presented at the 2021 NIH BRAIN Initiative Investigators Meeting \citep{bolenTorqueValuesArtificialKnee2021June15-17}. It did not achieve weight-bearing walking. This suspended stepping result preceded the later feedback and bilateral-RG extensions; it does not demonstrate walking with those additions.
~~~

**Work class:** Attribution/scope. **Next check:** Verify contributors and unilateral/bilateral scope from the source record; make the smallest wording correction.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-104"></a>

### M-104 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:51 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:10](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:10>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:51 am:**

~~~text
Huh? Just delete all this?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I adjusted applied currents, synaptic weights, and muscle parameters by hand. The suspended model developed alternating muscle activity and stepping-like motion, presented at the 2021 NIH BRAIN Initiative Investigators Meeting \citep{bolenTorqueValuesArtificialKnee2021June15-17}. It did not achieve weight-bearing walking. This suspended stepping result preceded the later feedback and bilateral-RG extensions; it does not demonstrate walking with those additions.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-105"></a>

### M-105 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:52 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:12](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:12>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:52 am:**

~~~text
Results? Lessons?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I subsequently expanded the neural model into left- and right-side circuits with interconnected rhythm generators (RGs). I constructed the right-side subsystems from the left-side organization and linked the two RGs through commissural pathways. This later model is referred to as the Walker Two-Layer (W2L) model. Further circuit schematics and implementation records are provided in Appendix~\ref{app:simulation_verification}.
~~~

**Work class:** Evidence/account. **Next check:** State the observed air-stepping result, the failure to achieve weight-bearing walking, and the specific next capability needed. Check existing traces before requesting any new simulation.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-106"></a>

### M-106 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 10:52 am  
**Section:** 6.1 The AnimatLab Bipedal Walker  
**Frozen source line:** [60-sns.tex:14](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:14>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 10:52 am:**

~~~text
huh?
~~~

**Advisor (ajh26) · 23 October, 11:11 am:**

~~~text
This paragraph confuses me.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The missing step was a tested link between spinal neural circuits and the human-scale and BPA-actuated mechanical models needed to study this robot, including the many joint degrees of freedom and 92 musculotendon actuators represented in Gait2392. Air stepping in the combined AnimatLab model did not establish how the same controller would behave with human musculotendon mechanics, robot BPA routing, or ground-contact feedback. I therefore needed models and tests that could separate the effects of the circuit, actuator mechanics, and numerical implementation.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** 0/4 on the prior agent clarity revision: explicit repeated reader criticism; not a claim that every sentence is factually wrong. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-107"></a>

### M-107 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:13 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: Within this project, five model combinations were missing: a
to I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-108"></a>

### M-108 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:15 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: coupled to a CPG in MuJoCo--SNS-Toolbox; its
to in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-109"></a>

### M-109 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:17 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: to the Gait2392 Model, also
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-110"></a>

### M-110 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:17 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: ; the
to controller. 3) I
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-111"></a>

### M-111 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:18 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: the
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-112"></a>

### M-112 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:18 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: controller
to neural model
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-113"></a>

### M-113 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:18 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: ;
to and tested it... 4) I built a
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-114"></a>

### M-114 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:19 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:19 am:**

~~~text
huh?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-115"></a>

### M-115 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:19 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: and
to 5) I coupled
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-116"></a>

### M-116 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:19 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted:  coupled
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-117"></a>

### M-117 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:19 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: These missing combinations prevented a
to These tools enable
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-118"></a>

### M-118 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:20 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Added: s
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-119"></a>

### M-119 · 60-sns.tex · Advisor tracked edit

**Author/date:** ajh26 · 23 October, 11:20 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Deleted:  the same
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-120"></a>

### M-120 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:21 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:21 am:**

~~~text
need a better word
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-121"></a>

### M-121 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:21 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:27](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:27>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:21 am:**

~~~text
*
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the Gait2392 human model in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated BPA-routed robot counterpart to the Gait2392 Model, also coupled to a CPG controller. 3) I combined the Li walker and Deng neural model with foot-contact feedback and tested it... 4) I built a graphical neural construction coupled to a maintained mechanical simulator; 5) I coupled the SolidWorks robot assembly to MATLAB BPA models and a neural controller. These tools enable controlled comparisons of neural organization across biological and artificial actuation. I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-122"></a>

### M-122 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:26 am  
**Section:** 6.2 Knowledge Gaps and the New Tools  
**Frozen source line:** [60-sns.tex:29](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:29>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:26 am:**

~~~text
This all still feels a bit awkward.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Existing software supplied separate parts of the solution. OpenSim supplied musculoskeletal geometry, MyoConverter translated it to MuJoCo, and SNS-Toolbox supplied neuron and synapse models \citep{delp_opensim_2007, caggiano_myosuite_2022, todorov_mujoco_2012, nourse_sns_2023}. The MuJoCo version used here supports scripted simulation and parameter searches but uses rigid tendons; its Python interface does not provide the block-based circuit assembly needed for this workflow. Simulink supplies that graphical construction, while Simscape Multibody imports mechanical geometry, including SolidWorks-derived mechanisms. AnimatLab retains the reference controllers but is unsupported, relies on an older physics engine, and ran these simulations slowly. The tools connect these complementary environments and use matched tests to check what transfers between them.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** 0/4 on the prior agent clarity revision: explicit repeated reader criticism; not a claim that every sentence is factually wrong. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-123"></a>

### M-123 · 60-sns.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:43 am  
**Section:** 6.3 Software Tools  
**Frozen source line:** [60-sns.tex:46](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:46>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: & Mechanism import achieved & Standalone reflex demonstrations do not establish control of imported models
to ...& Mechanism import achieved & The imported lower-body mechanism's ankle joints remain rigid, and standalone reflex demonstrations do not establish control of imported models
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
Mechanical imports & OpenSim-derived and SolidWorks-derived mechanisms in Simscape Multibody ...& Mechanism import achieved & The imported lower-body mechanism's ankle joints remain rigid, and standalone reflex demonstrations do not establish control of imported models  \\
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-124"></a>

### M-124 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:28 am  
**Section:** 6.3.1 Software Tools / Model Conversion and Actuator Coupling  
**Frozen source line:** [60-sns.tex:53](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:53>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:28 am:**

~~~text
I don't really know what this means
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
MyoConverter supplies the initial OpenSim-to-MuJoCo translation. I added repair and audit steps for the custom muscle paths and joint conventions, then coupled the converted model to an SNS-Toolbox controller. The converted model carries 92 muscle-tendon actuators. Geometry comparisons against the source and individual-actuator tests check that the conversion preserves the intended attachment paths and torque directions. The neural and mechanical models exchange signals at a \qty{2}{\milli\second} step.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-125"></a>

### M-125 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:28 am  
**Section:** 6.3.1 Software Tools / Model Conversion and Actuator Coupling  
**Frozen source line:** [60-sns.tex:53](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:53>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:28 am:**

~~~text
huh?
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
MyoConverter supplies the initial OpenSim-to-MuJoCo translation. I added repair and audit steps for the custom muscle paths and joint conventions, then coupled the converted model to an SNS-Toolbox controller. The converted model carries 92 muscle-tendon actuators. Geometry comparisons against the source and individual-actuator tests check that the conversion preserves the intended attachment paths and torque directions. The neural and mechanical models exchange signals at a \qty{2}{\milli\second} step.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-126"></a>

### M-126 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:31 am  
**Section:** 6.3.1 Software Tools / Model Conversion and Actuator Coupling  
**Frozen source line:** [60-sns.tex:53](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:53>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:31 am:**

~~~text
I feel like this all seems to imply that there were problems. Perhaps if you are more specific here, that would clarify. Example You ran MyoConverter and half hte origin points were wrong and needed to be corrected. The joint orientations were left hand instead of right hand. I dunno. Either there were problems, and it is a good thing you checked and fixed them. Or there weren't and it feels weird that you are calling out you created an 'audit' system.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
MyoConverter supplies the initial OpenSim-to-MuJoCo translation. I added repair and audit steps for the custom muscle paths and joint conventions, then coupled the converted model to an SNS-Toolbox controller. The converted model carries 92 muscle-tendon actuators. Geometry comparisons against the source and individual-actuator tests check that the conversion preserves the intended attachment paths and torque directions. The neural and mechanical models exchange signals at a \qty{2}{\milli\second} step.
~~~

**Work class:** Evidence/account. **Next check:** Replace generic audit language with actual checks and findings from artifacts. Hypothetical wrong origins or handedness in the advisor comment must not become claimed findings.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-127"></a>

### M-127 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 23 October, 11:32 am  
**Section:** 6.3.1 Software Tools / Model Conversion and Actuator Coupling  
**Frozen source line:** [60-sns.tex:55](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:55>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Highlighted text:** Not exposed by current review UI; source line below is context, not an exact highlight.

**Advisor (ajh26) · 23 October, 11:32 am:**

~~~text
What BPA interface? Non seq
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The BPA interface evaluates the quasi-static actuator force model from Chapter~\ref{ch:muscles} using the simulated muscle length and pressure. It provides an actuator model for controller development, rather than a simulation of valve and pressure dynamics. The controller receives length, velocity, force, and contact signals through the corresponding sensory channels. The correction terms of Chapter~\ref{ch:torque} do not by themselves verify the transferred routing or compliance of a dynamic model; these require separate checks.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Introduce inputs (simulated length and pressure), evaluated force model, output/coupling, and the valve-dynamics limitation in one local passage; verify this account against actual implementation.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-128"></a>

### M-128 · 60-sns.tex · Ben tracked edit

**Author/date:** bbolen83 · 23 October, 9:44 am  
**Section:** 6.4.3 Verification of the Implemented Interfaces / Closed-Loop Model Coupling  
**Frozen source line:** [60-sns.tex:109](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:109>) (review start; not a confirmed highlight range).  
**Anchor confidence:** Current review start mapped to frozen source; highlight end unavailable.

**Exact tracked change, including insertions/deletions without comments:**

~~~text
Changed: The MuJoCo--Simulink bridge was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} closed-loop run. Commands and measurements passed through the interface. However, under tonic drive alone the Simulink realization reached a fixed point while the NumPy reference oscillated. The network-level discrepancy remains a limit on transferring tuned controller parameters between implementations.
to The MuJoCo--Simulink bridge, built on the official MathWorks blockset that embeds the MuJoCo model inside Simulink as the mechanical model, was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} open-loop run, without afferent return; the closed-loop reflex experiment exercised commands and sensory measurements through the same interface. However, under tonic drive alone the Simulink realization reached a fixed point while the NumPy reference oscillated. The network-level discrepancy remains a limit on transferring tuned controller parameters between implementations.
~~~

**Current anchored source line (context, not asserted exact highlight):**

~~~latex
The MuJoCo--Simulink bridge, built on the official MathWorks blockset that embeds the MuJoCo model inside Simulink as the mechanical model, was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} open-loop run, without afferent return; the closed-loop reflex experiment exercised commands and sensory measurements through the same interface. However, under tonic drive alone the Simulink realization reached a fixed point while the NumPy reference oscillated. The network-level discrepancy remains a limit on transferring tuned controller parameters between implementations.
~~~

**Work class:** Tracked edit review. **Next check:** Check the exact before/after text, meaning, grammar, attribution, and factual scope. Recording is not acceptance.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-129"></a>

### M-129 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 11:48 am  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Two definitions cover most of the neuroscience vocabulary used below. First, the locomotor controller is a two-level central pattern generator: a rhythm generator sets the timing of the step cycle, and a pattern-formation layer shapes that timing into phase-appropriate activation of antagonist muscle groups \citep{rybak_modelling_2006}. Each leg's rhythm generator is a pair of mutually inhibiting antagonist populations, one extensor and one flexor, so that only one is strongly active at a time. Second, in the non-spiking (graded) formulation used for most of this work, each neuron is electrically a nonlinear RC circuit: synaptic current charges a membrane capacitance through a leak conductance, and the resulting membrane voltage, tens of millivolts, is both the neuron's state and its output \citep{nourse_sns_2023}. Graded chemical synapses couple these units: a depolarized presynaptic cell opens a conductance in the postsynaptic membrane whose sign makes the connection excitatory or inhibitory. The sensory channels that close the loop are the muscle afferents, Ia signaling muscle length and velocity, Ib signaling muscle force, and group II signaling length, together with heel and toe contact sensors.

**Earlier locator:** SNS-01, source line 6 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 11:48 am:**

~~~text
So many details and I don't understand what the goal is. What is the gap in knowledge? What did you do to fill it?
~~~

**Ben (bbolen83) · 21 October, 2:09 pm:**

~~~text
Gaps: Gait2392_simbody with CPG in Mujoco-SNS Toolbox
- Gait2392_robotbody (BPAs and BPA routing) with CPG
- Block-programming on a fast physics simulator (Animatlab no longer supported)
- Robot solidworks assembly with BPAs from Matlab and a neural controller
- Li walker with Deng neural controller and foot contact feedback.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I rebuilt the opening around the Li/Deng AnimatLab walker and my interconnected left/right RG extension. It then identifies the five missing model/controller combinations in my earlier reply and explains the tools and verification needed to investigate them. Detailed supporting results are in Appendix D.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-130"></a>

### M-130 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 12:45 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): the reflex-based and learned alternatives \citep{geyer_musclereflex_2010, song_deep_2021} are treated there.

**Earlier locator:** SNS-02, source line 8 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 12:45 pm:**

~~~text
huh?
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I removed the build-order shorthand. The opening now connects the bilateral AnimatLab construction to the missing interfaces and the verification each new tool provides.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-131"></a>

### M-131 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 12:49 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): built in any of several simulation engines

**Earlier locator:** SNS-06, source line 10 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 12:49 pm:**

~~~text
why would we need several?
~~~

**Ben (bbolen83) · 21 October, 2:04 pm:**

~~~text
Well, Mujoco doesn't have as good of a physics engine as HYFYDY, for example, and is limited to rigid tendons. Python also is harder to use because it is code without a GUI. Simscape multibody allows for importing parts from solidworks. Some sort of block-programming language like Simulink is a natural choice for building our style of connectome. Mujoco-SNS Toolbox allows for real-time control and computer in the loop commands.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I now explain the complementary needs: scripted neural/mechanical coupling, graphical neural construction, and SolidWorks/MATLAB integration. The capability table separates demonstrated interfaces from unfinished combinations, and the text states the relevant simulator limitations.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-132"></a>

### M-132 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 12:53 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): answer

**Earlier locator:** SNS-07, source line 10 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 12:53 pm:**

~~~text
This is all super awkward. You need to identify the gap in knowledge and how you fill that gap in knowledge. Much of what you wrote here should be saved for discussion, where you talk about the implications of the work.
~~~

**Ben (bbolen83) · 21 October, 2:00 pm:**

~~~text
I thought that was for research hypothesis not for "here are some tools we made".
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I put the knowledge gaps before the new tools and kept the main results to representative verification. Broader sweeps and supported-walker results are in Appendix D; future integration belongs in the Discussion chapter.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-133"></a>

### M-133 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 12:59 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): verification-first approach

**Earlier locator:** SNS-09, source line 20 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 12:59 pm:**

~~~text
I don't understand the connection.
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I rewrote the transition to connect the identified modeling gaps to the specific neural/mechanical interfaces developed.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-134"></a>

### M-134 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 1:00 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): \begin{figure}[htbp]

**Earlier locator:** SNS-11, source line 26 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 1:00 pm:**

~~~text
Overlapping text in figure
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I replaced the pipeline schematic with the corrected vector figure, with separated columns and labels.
~~~

**Work class:** Figure/layout. **Next check:** Check the actual rendered figure and the reader-facing purpose; no model experiment required.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-135"></a>

### M-135 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 1:02 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): build order

**Earlier locator:** SNS-10, source line 24 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 1:02 pm:**

~~~text
huh?
~~~

**Ben (bbolen83) · 23 October, 6:27 am:**

~~~text
I removed the unexplained inventory and retained the interface’s purpose, its demonstrated check, and its limitation.
~~~

**Work class:** Clarification; scope to verify. **Next check:** Read the actual marked passage and surrounding argument. Prefer a local definition, specific referent, or transition; escalate only if the claim lacks evidence.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-136"></a>

### M-136 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 1:08 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): \subsection{New Tools}\label{sec:sns_new_tools}

**Earlier locator:** SNS-14, source line 39 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 1:08 pm:**

~~~text
Wait, that is all you are saying about the other tools? It is hard for the reader to understand the limitations and gaps without a deeper description of them.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I expanded the explanation of the environments’ capabilities and limitations and added a table showing which interfaces are demonstrated and which combinations remain incomplete.
~~~

**Work class:** Substantial explanation. **Next check:** Explain the missing capability/knowledge, why it mattered, the work performed, and its demonstrated outcome in reading order.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.

<a id="m-137"></a>

### M-137 · 60-sns.tex · Comment thread

**Author/date:** ajh26 · 21 October, 1:11 pm  
**Section:** Original current section unavailable: old anchor collapsed at end of replaced file.  
**Frozen source line:** [60-sns.tex:117](<D:/Github/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/matrix_source/60-sns.tex:117>) (end-of-file anchor only, not original location).  
**Anchor confidence:** Collapsed at end of replaced file; current original location unavailable.

**Highlighted text:** Historical exact highlight (22 Oct. capture): Simulink blocks

**Earlier locator:** SNS-15, source line 41 in the 22 October capture; not the current source line.

**Advisor (ajh26) · 21 October, 1:11 pm:**

~~~text
Did you compare and discuss against the libraries myself and Fletcher made here?
~~~

**Ben (bbolen83) · 21 October, 2:00 pm:**

~~~text
No. You've never encouraged us to use Simulink.
~~~

**Ben (bbolen83) · 23 October, 6:21 am:**

~~~text
I removed the Hunt/Young library comparison because those citations do not establish a shared neuromechanical Simulink library. The section now describes the implemented blocks and matched-circuit verification directly. A comparison with the specific libraries you mean still requires checking those libraries.
~~~

**Work class:** Evidence/comparison. **Next check:** Identify the actual earlier libraries, compare functions and limitations; do not substitute an unrelated citation or performance claim.

**Self-score:** Unassessed: no new revision was made or tested in this pass. **Actual new revision/evidence:** no source edit in this capture; scores of 0 concern prior agent clarity work and the quoted reader feedback. **Acceptance:** not asserted.


---

# Earlier snapshots and accountability record (preserved verbatim)

# Final review capture — 10/23/2026, 12:18:55 PM Pacific

Ben requested cancellation of the half-hour monitor after the advisor left for the day. **The automation was deleted successfully.** This final read-only pass found **no new, changed, or missing Overview entries** compared with the latest saved snapshot.

Final open-review inventory: **63 advisor comment threads, 63 pending advisor tracked edits, and 11 pending Ben tracked edits** (137 entries). All expandable text was opened before comparison. The full current capture is saved as Notes/advisor_review_20261023/entries_final_20261023.json; the dated comment/edit additions, suggestions, and repeated-failure audit remain below. This confirms capture of the open review inventory, not advisor approval, acceptance of tracked edits, or a re-audit of resolved comments.

No Overleaf source, comments, replies, resolutions, or tracked-change decisions were changed. No source anchors were followed or PDF scrolling performed. The temporary review tab was closed. Only local report/evidence files were written; no commits or pushes. Scheduled monitoring is stopped.

---

# Monitoring update — 10/23/2026, 11:43:33 AM Pacific

**Eight new advisor threads and twelve new pending advisor edits, all in 60-sns.tex.** Current inventory: **63 open advisor threads, 63 pending advisor edits, and 11 pending Ben edits** (137 Overview entries). Prior snapshots remain intact below.

## What changed and what it means

The advisor is rewriting the five-system introduction into explicit first-person actions: what Ben adapted, developed, connected, and tested. That is concrete guidance on how to repair the earlier abstract “missing combinations” framing. These remain pending edits, not accepted changes or independently verified claims that all five systems are complete. One replacement ends “and tested it...”; it requires an actual verified test description before finalizing.

New comments include “This all still feels a bit awkward.” (11:26), “I don't really know what this means” and “huh?” (11:28), and “What BPA interface? Non seq” (11:32). **Reinforce the major deductions for unclear SNS explanations and disconnected prose.** These repeat concerns already raised on 21 October and earlier today. The word “still” explicitly signals continuing dissatisfaction. No further intervention by an agent on those exact paragraphs was verified between comments.

At 11:31 the advisor explains the audit problem in detail: name actual defects found and corrected, or explain the checks without implying that defects existed. His examples of wrong muscle origins and left-handed joint orientations are expressly hypothetical (“I dunno”), **not findings to copy into the dissertation**. This adds a specific deduction against calling a checking procedure an “audit” without explaining its concrete outcome.

## Suggested local dispositions

- **Five-system introduction:** follow the advisor's concrete action/result structure, but verify each system's completed scope and test. Do not turn proposed or partial integration into a completed system because the pending wording says “I built”.
- **“audit” passage:** describe the specific properties checked, then the observed result. If a defect was found, identify the real defect, correction, and confirming check from campaign evidence. If checks passed unchanged, say that. Do not invent the advisor's illustrative problems.
- **“What BPA interface? Non seq”:** define the actual inputs, outputs, and mechanical/neural connection at first mention, and place the explanation where that interface is introduced. Exact marked text was not recovered in this passive check, so no replacement is asserted ready to paste.
- **“need a better word” and “*”:** exact highlighted targets remain unverified; retain both comments for contextual review rather than guessing the requested word or interpreting the asterisk.
- **“This all still feels a bit awkward”:** reopen the paragraph-level revision. A passing build and an added definition do not close this concern. Check the whole explanation in reading order.

## New review items

### CHECK-1141-01 · 60-sns.tex · edit · position 3650

~~~text
ajh26
23 October, 11:13 am
Changed: Within this project, five model combinations were missing: a
to I developed several software tools to aid in the development of a human-mimetic robot driven by a neural model of control. To do this, I built the following five systems. 1) I adapted the
~~~

### CHECK-1141-02 · 60-sns.tex · edit · position 3859

~~~text
ajh26
23 October, 11:15 am
Changed: coupled to a CPG in MuJoCo--SNS-Toolbox; its
to in OpenSIM to an equivalent model in the physics simulator MuJoCo. I then connected it to the SNS-Toolbox and tested it with a CPG controller. 2) I developed a simulated
~~~

### CHECK-1141-03 · 60-sns.tex · edit · position 4058

~~~text
ajh26
23 October, 11:17 am
Added: to the Gait2392 Model, also
~~~

### CHECK-1141-04 · 60-sns.tex · edit · position 4102

~~~text
ajh26
23 October, 11:17 am
Changed: ; the
to controller. 3) I
~~~

### CHECK-1141-05 · 60-sns.tex · edit · position 4129

~~~text
ajh26
23 October, 11:18 am
Added: the
~~~

### CHECK-1141-06 · 60-sns.tex · edit · position 4152

~~~text
ajh26
23 October, 11:18 am
Changed: controller
to neural model
~~~

### CHECK-1141-07 · 60-sns.tex · edit · position 4191

~~~text
ajh26
23 October, 11:18 am
Changed: ;
to and tested it... 4) I built a
~~~

### CHECK-1141-08 · 60-sns.tex · comment · position 4232

~~~text
ajh26
23 October, 11:19 am
huh?
~~~

### CHECK-1141-09 · 60-sns.tex · edit · position 4298

~~~text
ajh26
23 October, 11:19 am
Changed: and
to 5) I coupled
~~~

### CHECK-1141-10 · 60-sns.tex · edit · position 4340

~~~text
ajh26
23 October, 11:19 am
Deleted:  coupled
~~~

### CHECK-1141-11 · 60-sns.tex · edit · position 4387

~~~text
ajh26
23 October, 11:19 am
Changed: These missing combinations prevented a
to These tools enable
~~~

### CHECK-1141-12 · 60-sns.tex · edit · position 4427

~~~text
ajh26
23 October, 11:20 am
Added: s
~~~

### CHECK-1141-13 · 60-sns.tex · edit · position 4431

~~~text
ajh26
23 October, 11:20 am
Deleted:  the same
~~~

### CHECK-1141-14 · 60-sns.tex · comment · position 4646

~~~text
ajh26
23 October, 11:21 am
need a better word
~~~

### CHECK-1141-15 · 60-sns.tex · comment · position 4695

~~~text
ajh26
23 October, 11:21 am
*
~~~

### CHECK-1141-16 · 60-sns.tex · comment · position 4749

~~~text
ajh26
23 October, 11:26 am
This all still feels a bit awkward.
~~~

### CHECK-1141-17 · 60-sns.tex · comment · position 7520

~~~text
ajh26
23 October, 11:28 am
I don't really know what this means
~~~

### CHECK-1141-18 · 60-sns.tex · comment · position 7551

~~~text
ajh26
23 October, 11:28 am
huh?
~~~

### CHECK-1141-19 · 60-sns.tex · comment · position 7777

~~~text
ajh26
23 October, 11:31 am
I feel like this all seems to imply that there were problems. Perhaps if you are more specific here, that would clarify. Example You ran MyoConverter and half hte origin points were wrong and needed to be corrected. The joint orientations were left hand instead of right hand. I dunno. Either there were problems, and it is a good thing you checked and fixed them. Or there weren't and it feels weird that you are calling out you created an 'audit' system.
~~~

### CHECK-1141-20 · 60-sns.tex · comment · position 7957

~~~text
ajh26
23 October, 11:32 am
What BPA interface? Non seq
~~~


## Verification and scope

All 137 Overview entries were captured with expandable text open. Evidence: Notes/advisor_review_20261023/entries_20261023_1141.json. Compared against entries_latest.json; twenty new items, with no previously captured text changes detected. A temporary separate tab was used and closed; no comment anchors were followed and no PDF scrolling was performed. No Overleaf source, replies, resolutions, or accept/reject decisions changed. Only local report/evidence files were written; no commit, push, or LaTeX-source changes.

---

# Monitoring update — 10/23/2026, 11:13:13 AM Pacific

One existing SNS thread gained an advisor follow-up; no new threads or pending edits were detected. Current inventory remains 55 open advisor threads, 51 pending advisor edits, and 11 pending Ben edits (117 entries).

## Repeated clarity failure: SNS paragraph

The 23 October, 10:52 a.m. “huh?” thread in 60-sns.tex now also contains, at 11:11 a.m.: **“This paragraph confuses me.”** This is a second explicit criticism within the same thread. It reinforces the earlier accountability deduction for SNS clarity. No intervening revision to that exact paragraph was verified, so this is repeated feedback, not proof that a further attempted fix failed.

Suggested local disposition: reopen this paragraph as NOT CLEAR TO THE REVIEWER. Recover its exact anchor before proposing replacement prose. Explain its single purpose and concrete result, or propose removing it if it repeats material elsewhere. Another abstract transition or expanded tool inventory would not answer the criticism. No paragraph was edited or deleted in this check.

~~~text
ajh26
23 October, 10:52 am
huh?
ajh26
23 October, 11:11 am
This paragraph confuses me.
~~~

Evidence: Notes/advisor_review_20261023/entries_20261023_1111.json. All Overview text was expanded in a temporary separate tab; no comment anchors were followed, no PDF scrolling was performed, and no Overleaf content or review decisions changed. Previous report retained verbatim below.

---

# Advisor recheck and repeated-failure audit — 10/23/2026, 11:08:20 AM Pacific

**12 new advisor threads and 10 new pending advisor edits** since the 10:11 snapshot. Six new threads are in Chapter 5 and six in SNS. Current inventory: **55 open advisor threads, 51 pending advisor edits, and 11 pending Ben edits** (117 Overview entries). Four newly observed Ben edits are recorded separately, not attributed to the advisor. All prior report text is preserved below.

## Accountability: repeated criticisms that our revisions did not adequately fix

These are deductions against the agent writing/review process, including my reporting. They are not deductions against Ben. A retained old thread alone does not prove a failed second attempt; the cases below distinguish a later recurrence from several comments made in one reading. The evidence does not reliably identify which individual agent authored every sentence, so it would be misleading to assign every defect specifically to ZCode. The prior Revision_C/COMMENT_DISPOSITIONS.md and this report's completion claims are part of the record being challenged, not evidence of successful communication.

| Deduction | Earlier criticism and claimed fix | Repeated criticism | Judgment and correction needed |
|---|---|---|---|
| **Major: database overview still fails the reader** | 21 Oct., 11:39: “You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized?” Revision C DB-22 says the opening explains the resource, five record types, and organization (COMMENT_DISPOSITIONS.md:117). Ben's 23 Oct., 6:21 reply makes the same claim. | 23 Oct., 9:51: “This is unclear”; 9:52: “organization is still unclear”; 10:26: “I don't really understand this. Also, be clear about what a paper cites, vs what cites the paper”. | **Failed revision at the reader-understanding level.** Listing tables and mentioning a screenshot did not explain the relationships. Give one concrete paper-to-pathway-to-review example, then distinguish outgoing references from incoming citations. Recheck against the actual schema. Do not mark this addressed merely because a definition was inserted. |
| **Major: unexplained terminology keeps returning** | Chapter 5 already had “what are notes?” and several “huh?” threads. DB-10–12 claim definitions and removal of jargon (disposition lines 105–107); DB-21 claims notes are defined (line 116). | 23 Oct., 9:52: “Unclear why there would be an 'archived' record.” At 10:22 two further “huh?” threads, one followed by “sensory feedback pathway?”; at 10:29 “animal?” | **Repeated clarity failure, with exact new highlighted phrases not recovered in this pass.** Define each category at its first use, state what it contains, and provide an example. Archive is a disposition and needs a reason; animal may mean species, subject, or model and must be specified from context. Do not substitute another abstract label for the one questioned. |
| **Major: SNS still lacks a clear account of results and lessons** | 21 Oct., 11:48 asks for the goal, knowledge gap, and contribution. At 12:53 the advisor again asks how the gap is filled and says implications belong in Discussion. SNS-01/SNS-07 claim the causal opening and verification account resolve this (disposition lines 122/128). Ben's 23 Oct., 6:21 replies describe the changes. | 23 Oct., 10:51: “Huh? Just delete all this?”; 10:52: “Results? Lessons?” and another “huh?”. | **The rewrite has not reliably communicated the purpose and outcome.** The deletion is a proposal to disposition locally, not an instruction executed here. Replace process narration with what was attempted, what occurred, what the attempt established, and what remained unachieved. State the known AnimatLab outcome plainly: stepping-like movement in air; hand tuning did not achieve walking. Tools and workflow are not a substitute for results. |
| **Major: disconnected prose is a recurring editing defect** | Background includes 22 Oct., 11:49 “non seq”, 11:56 “I feel like this is jumping all over the place”, and 12:13 “non seq”. BG-19/BG-43/BG-50 claim repaired purposes and transitions (disposition lines 64/87/94). | 23 Oct., 8:56 “Non seq”, plus new SNS “Huh? Just delete all this?” at 10:51. | **Recurring defect across sections; exact same-sentence recurrence is not established.** We revised locally without consistently checking the paragraph's place in the argument. Every paragraph needs a reason to follow its predecessor. Do not add another transition sentence around content that does not belong. |
| **Moderate: contribution attribution and pronouns need an explicit check** | Ben's 10:10 reply to Introduction's “double check”: “I didn't redesign the linkages”. The historical report had already emphasized separating completed work and attribution. | SNS has two distinct “We?” threads at 10:49 and 10:50, and “for each leg?” at 10:49. | **At least two direct pronoun criticisms in this pass; cross-pass recurrence is not proven for these exact pronouns.** Name contributors or use “I” for Ben's individual work and “we” only for documented collaborative work. Distinguish actuator routing from linkage design. Verify bilateral scope before supplying a count. |
| **Major: completion reporting was too confident** | Prior report says every advisor ID has a local disposition; the disposition table repeatedly says a concept is explained. Earlier records use four-gate GO judgments. | The later comments above again identify unclear organization, language, missing results, and attribution. | **Deduction against our review standard, including my summaries.** A successful build, implemented wording, or populated disposition row cannot establish that the advisor's concern is resolved. Future status should be “revised locally; reader acceptance unverified” unless the evidence supports more. In particular, DB-22 and the SNS opening claims need reopening as editorial concerns. No old record has been silently rewritten. |

The strongest documented second-hit cases are database organization and SNS purpose/results. The repeated “huh?” and “non seq” pattern is also clear, but not every short comment can be assigned to the same sentence or a particular agent without exact anchors and authorship history.

## Suggested local fixes for the newly observed comments

These are proposals, not changes made to Overleaf. Exact new highlights were not followed to avoid moving the PDF. Review positions and full wording are in the ledger below.

- **Citation direction (Chapter 5, 10:26):** possible wording: “A paper's references identify earlier works cited by that paper. Its citing papers are later works that cite it.” Then explain which of those relationships the application actually stores or displays. Verify implementation before describing automatic links or completeness.
- **“expand” (Chapter 5, 10:21):** recover the highlighted source before drafting an expansion. Do not assume this means add more general database history.
- **Two “huh?” threads (Chapter 5, 10:22):** one explicitly asks “sensory feedback pathway?” Define a pathway by the actual stored fields and a source-linked example; avoid implying the database fully specifies a physiological circuit. Exact intended wording remains uncertain.
- **“animal?” (Chapter 5, 10:29):** replace the ambiguous term with species, experimental subject, or musculoskeletal model as appropriate to the marked passage. Context is required before choosing one.
- **Future work (Chapter 5, 10:34):** the advisor explicitly requests connecting feedback pathways with papers, and reviews with both papers and pathways. Suggested direction: “Further work should link each feedback pathway to the papers supporting it and connect each review to the papers and pathways it discusses.” Verify what links already exist so completed relationships are not incorrectly described as missing. Discuss remaining linking/verification work concretely rather than promising exhaustive coverage.
- **SNS “for each leg?” (10:49):** specify whether the described controller structure is unilateral, duplicated bilaterally, or interconnected across legs. No neuron count or bilateral claim is inferred from the question alone.
- **SNS two “We?” comments:** identify the actual people responsible before choosing the pronoun. Do not automatically replace all instances with “I”.
- **SNS “Results? Lessons?”:** pair the AnimatLab attempt with its known outcome and resulting need. Possible wording: “The combined model produced stepping-like movement in air, but hand tuning did not achieve walking. This limitation motivated tools for inspecting neural activity and systematically testing controller parameters.” Check that the cited attempt specifically motivated the named tools before using the causal sentence.

## New or changed review items since 10:11

### RECHECK-01 · 30-torque.tex · pending Ben edit · position 49483

~~~text
bbolen83
23 October, 10:56 am
Changed: \begin{sidewaystable}[p] \centering\small\setlength{\tabcolsep}{4pt} \caption[Pre-optimization pinned-flexor goodness of fit]{Pre-optimization pinned-flexor goodness of fit. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. Scoring-point restrictions are stated in the text.} \label{tab:FlxPinPreGof} \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}lccccccccc@{}} \toprule & \multicolumn{3}{c}{\textbf{Rigid body}} & \multicolumn{3}{c}{\textbf{Improved char.\ alone}} & \multicolumn{3}{c}{\textbf{Back-calculated (hybrid)}} \\ \cmidrule(lr){2-4} \cmidrule(lr){5-7} \cmidrule(l){8-10} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{48.5}{\cm}$ (training) & 5.779 & 1.1761 & 9.798 & 4.871 & 0.8354 & 9.865 & 1.635 & 0.0991 & 3.842 \\ $l_{0}=\qty{41.0}{\cm}$ (training) & 6.633 & 3.2313 & 11.409 & 6.074 & 2.7090 & 9.262 & 2.187 & 0.1184 & 4.243 \\ $l_{0}=\qty{45.7}{\cm}$ (validation) & 5.315 & 0.9408 & 9.991 & 4.465 & 0.6640 & 9.536 & 1.246 & 0.0409 & 2.275 \\ $l_{0}=\qty{40.6}{\cm}$, tendon (validation) & 3.945 & 0.6750 & 7.527 & 3.226 & 0.4512 & 7.376 & 2.085 & 0.0755 & 3.073 \\ \midrule Mean (training) & 6.206 & 2.2037 & 10.604 & 5.472 & 1.7722 & 9.563 & 1.911 & 0.1088 & 4.043 \\ Mean (validation) & 4.630 & 0.8079 & 8.759 & 3.845 & 0.5576 & 8.456 & 1.666 & 0.0582 & 2.674 \\ \bottomrule \end{tabular} \end{sidewaystable}
to \begin{table}[H] \caption{Pre-optimization pinned-knee flexor torque-prediction goodness of fit per test, comparing the three pre-optimization torque predictions of Figure~\ref{fig:FlxPin_preopt}: the Hunt-equation rigid-body prediction, the unoptimized prediction of the improved BPA characterization, and the back-calculated (hybrid) torque. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:FlxPinPreGof} \begingroup \renewcommand{\arraystretch}{1}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc ccccccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Hunt equation}} & \multicolumn{3}{c}{\textbf{Improved char.\ alone}} & \multicolumn{3}{c}{\textbf{Back-calculated (hybrid)}} \\ \cmidrule(lr){3-5} \cmidrule(lr){6-8} \cmidrule(l){9-11} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule 48.5 & -- & 5.779 & 1.1761 & 9.798 & 4.871 & 0.8354 & 9.865 & 1.635 & 0.0991 & 3.842 \\ 41.0 & -- & 6.633 & 3.2313 & 11.409 & 6.074 & 2.7090 & 9.262 & 2.187 & 0.1184 & 4.243 \\ 45.7 & -- & 5.315 & 0.9408 & 9.991 & 4.465 & 0.6640 & 9.536 & 1.246 & 0.0409 & 2.275 \\ 40.6 & 35 & 3.945 & 0.6750 & 7.527 & 3.226 & 0.4512 & 7.376 & 2.085 & 0.0755 & 3.073 \\ \bottomrule \end{tabular} \endgroup \end{table}
~~~

### RECHECK-02 · 30-torque.tex · pending Ben edit · position 53187

~~~text
bbolen83
23 October, 10:59 am
Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-flexor goodness of fit]{Optimized pinned-flexor goodness of fit at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.} \label{tab:FlxPinGof} \begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}} \toprule & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){2-4} \cmidrule(l){5-7} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{48.5}{\cm}$ (training) & 7.617 & 2.1278 & 25.531 & 2.165 & 0.1718 & 5.467 \\ $l_{0}=\qty{41.0}{\cm}$ (training) & 9.694 & 2.1149 & 21.076 & 1.084 & 0.0277 & 1.766 \\ $l_{0}=\qty{45.7}{\cm}$ (validation) & 8.803 & 1.9380 & 26.264 & 1.368 & 0.0468 & 3.614 \\ $l_{0}=\qty{40.6}{\cm}$, tendon (validation) & 6.373 & 0.5935 & 16.859 & 1.611 & 0.0427 & 5.891 \\ \midrule Mean (training) & 8.656 & 2.1213 & 23.303 & 1.624 & 0.0998 & 3.616 \\ Mean (validation) & 7.588 & 1.2657 & 21.562 & 1.490 & 0.0447 & 4.752 \\ \bottomrule \end{tabularx} \end{table}
to \begin{table}[H] \caption{Pinned-knee flexor torque-prediction goodness of fit: the improved BPA characterization alone (the zero-correction prediction of the text) and the corrected model at the adopted correction terms (Table~\ref{tab:AdoptedXi}); training and validation tests follow the adopted cross-validation fold. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:FlxPinGof} \begingroup \renewcommand{\arraystretch}{1.4}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc cccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){3-5} \cmidrule(l){6-8} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule \multicolumn{8}{@{}l}{\textbf{Training}} \\ 48.5 & -- & 7.617 & 2.1278 & 25.531 & 2.165 & 0.1718 & 5.467 \\ 41.0 & -- & 9.694 & 2.1149 & 21.076 & 1.084 & 0.0277 & 1.766 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 8.656 & 2.1213 & 23.303 & 1.624 & 0.0998 & 3.617 \\ \midrule \multicolumn{8}{@{}l}{\textbf{Validation}} \\ 45.7 & -- & 8.803 & 1.9380 & 26.264 & 1.368 & 0.0468 & 3.614 \\ 40.6 & 35 & 6.373 & 0.5935 & 16.859 & 1.611 & 0.0427 & 5.891 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 7.588 & 1.2657 & 21.562 & 1.490 & 0.0447 & 4.752 \\ \bottomrule \end{tabular} \endgroup \end{table}
~~~

### RECHECK-03 · 30-torque.tex · pending Ben edit · position 59478

~~~text
bbolen83
23 October, 11:01 am
Changed: \begin{table}[p] \centering\small \caption[Optimized pinned-extensor goodness of fit]{Optimized pinned-extensor goodness of fit at the adopted correction terms. RMSE and maximum absolute residual (Max.) are in N\,m; FVU is dimensionless. The baseline uses the improved BPA characterization without transmission corrections.} \label{tab:ExtPinGof} \begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}*{6}{>{\centering\arraybackslash}X}@{}} \toprule & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){2-4} \cmidrule(l){5-7} \textbf{Test} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule $l_{0}=\qty{43.6}{\cm}$ (training) & 2.437 & 0.9618 & 4.284 & 0.978 & 0.1634 & 1.818 \\ $l_{0}=\qty{43.2}{\cm}$, tendon (training) & 2.878 & 4.7304 & 4.817 & 1.525 & 1.3270 & 2.481 \\ $l_{0}=\qty{45.7}{\cm}$ (training) & 2.178 & 1.0586 & 4.974 & 0.490 & 0.0535 & 1.206 \\ $l_{0}=\qty{40.0}{\cm}$ (validation) & 3.092 & 0.5638 & 4.842 & 1.244 & 0.1303 & 2.196 \\ $l_{0}=\qty{40.0}{\cm}$, tendon (validation) & 3.622 & 6.3919 & 7.845 & 1.350 & 0.8878 & 1.890 \\ $l_{0}=\qty{46.5}{\cm}$ (validation) & 0.949 & 0.3152 & 2.272 & 0.699 & 0.1714 & 1.559 \\ \midrule Mean (training) & 2.498 & 2.2502 & 4.692 & 0.997 & 0.5146 & 1.835 \\ Mean (validation) & 2.554 & 2.4236 & 4.986 & 1.098 & 0.3965 & 1.882 \\ \bottomrule \end{tabularx} \end{table}
to \begin{table}[H] \caption{Pinned-knee extensor torque-prediction goodness of fit, at the adopted correction terms with the effective-stiffness pair locked to the flexor-identified values (Table~\ref{tab:AdoptedXi}): the improved BPA characterization alone (the zero-correction prediction of the text) and the corrected model. RMSE and maximum absolute residual (Max.) are in \unit{\newton\meter}, and FVU denotes the fraction of variance unexplained.} \label{tab:ExtPinGof} \begingroup \renewcommand{\arraystretch}{1}% \fontsize{10pt}{12pt}\selectfont \begin{tabular}{@{}cc cccccc@{}} \toprule & & \multicolumn{3}{c}{\textbf{Improved BPA char.\ alone}} & \multicolumn{3}{c}{\textbf{Corrected model}} \\ \cmidrule(lr){3-5} \cmidrule(l){6-8} $\boldsymbol{l_{0}}$ (\unit{\cm}) & $\boldsymbol{l_{\text{ten}}}$ (\unit{\mm}) & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} & \textbf{RMSE} & \textbf{FVU} & \textbf{Max.} \\ \midrule \multicolumn{8}{@{}l}{\textbf{Training}} \\ 43.6 & -- & 2.437 & 0.9618 & 4.284 & 0.978 & 0.1634 & 1.818 \\ 43.2 & 22 & 2.878 & 4.7304 & 4.817 & 1.525 & 1.3270 & 2.481 \\ 45.7 & -- & 2.178 & 1.0586 & 4.974 & 0.490 & 0.0535 & 1.206 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 2.498 & 2.2503 & 4.692 & 0.997 & 0.5146 & 1.835 \\ \midrule \multicolumn{8}{@{}l}{\textbf{Validation}} \\ 40.0 & -- & 3.092 & 0.5638 & 4.842 & 1.244 & 0.1303 & 2.196 \\ 40.0 & 40 & 3.622 & 6.3919 & 7.845 & 1.350 & 0.8878 & 1.890 \\ 46.5 & -- & 0.949 & 0.3152 & 2.272 & 0.699 & 0.1714 & 1.559 \\ \multicolumn{2}{@{}l}{\hspace*{1em}Mean} & 2.554 & 2.4236 & 4.986 & 1.098 & 0.3965 & 1.882 \\ \bottomrule \end{tabular} \endgroup \end{table}
~~~

### RECHECK-04 · 30-torque.tex · pending Ben edit · position 71969

~~~text
bbolen83
23 October, 10:53 am
Deleted: +
~~~

### RECHECK-05 · 50-database.tex · open thread · position 2911

~~~text
ajh26
23 October, 10:26 am
I don't really understand this. Also, be clear about what a paper cites, vs what cites the paper
~~~

### RECHECK-06 · 50-database.tex · pending advisor edit · position 3328

~~~text
ajh26
23 October, 10:19 am
Changed: Zotero libraries
to libraries
~~~

### RECHECK-07 · 50-database.tex · pending advisor edit · position 3338

~~~text
ajh26
23 October, 10:19 am
Changed: of the AARL and myself, and papers already curated in Airtable
to used to provide background and understanding for the previous publications on locomotion and postural balance control in the AARL lab
~~~

### RECHECK-08 · 50-database.tex · open thread · position 3537

~~~text
ajh26
23 October, 10:21 am
expand
~~~

### RECHECK-09 · 50-database.tex · pending advisor edit · position 3554

~~~text
ajh26
23 October, 10:21 am
Deleted:  and the laboratory's PDF collection
~~~

### RECHECK-10 · 50-database.tex · pending advisor edit · position 4066

~~~text
ajh26
23 October, 10:24 am
Added:  and placed in an archive
~~~

### RECHECK-11 · 50-database.tex · open thread · position 4093

~~~text
ajh26
23 October, 10:22 am
huh?
~~~

### RECHECK-12 · 50-database.tex · open thread · position 4181

~~~text
ajh26
23 October, 10:22 am
huh?
ajh26
23 October, 10:22 am
sensory feedback pathway?
~~~

### RECHECK-13 · 50-database.tex · pending advisor edit · position 4287

~~~text
ajh26
23 October, 10:23 am
Deleted: ; missing text was flagged rather than supplied from recollection
~~~

### RECHECK-14 · 50-database.tex · pending advisor edit · position 4289

~~~text
ajh26
23 October, 10:24 am
Deleted: The retained and archived sets distinguish the working collection from the broader audit.
~~~

### RECHECK-15 · 50-database.tex · pending advisor edit · position 5649

~~~text
ajh26
23 October, 10:29 am
Added: of the database
~~~

### RECHECK-16 · 50-database.tex · open thread · position 5682

~~~text
ajh26
23 October, 10:29 am
animal?
~~~

### RECHECK-17 · 50-database.tex · pending advisor edit · position 6569

~~~text
ajh26
23 October, 10:31 am
Added:  \section{Discussion}
~~~

### RECHECK-18 · 50-database.tex · pending advisor edit · position 7192

~~~text
ajh26
23 October, 10:33 am
Added:  \subsection{Further improvement and work}
~~~

### RECHECK-19 · 50-database.tex · open thread · position 7193

~~~text
ajh26
23 October, 10:34 am
Expand what all still needs to be done. Some thoughts. Connecting feedback pathways and papers. Connecting reviews with papers and feedback pathways.
~~~

### RECHECK-20 · 60-sns.tex · pending advisor edit · position 789

~~~text
ajh26
23 October, 10:48 am
Added: For a course project in ME 610: Neuromechanical Modeling,
~~~

### RECHECK-21 · 60-sns.tex · open thread · position 1364

~~~text
ajh26
23 October, 10:49 am
for each leg?
~~~

### RECHECK-22 · 60-sns.tex · open thread · position 1371

~~~text
ajh26
23 October, 10:49 am
We?
~~~

### RECHECK-23 · 60-sns.tex · open thread · position 1558

~~~text
ajh26
23 October, 10:50 am
We?
~~~

### RECHECK-24 · 60-sns.tex · open thread · position 1882

~~~text
ajh26
23 October, 10:51 am
Huh? Just delete all this?
~~~

### RECHECK-25 · 60-sns.tex · open thread · position 2333

~~~text
ajh26
23 October, 10:52 am
Results? Lessons?
~~~

### RECHECK-26 · 60-sns.tex · open thread · position 2455

~~~text
ajh26
23 October, 10:52 am
huh?
~~~

## Verification and scope

Captured 117 rendered Overview entries; all expandable text was opened. Compared full visible text with entries_latest.json from the previous check. Evidence: Notes/advisor_review_20261023/entries_manual_20261023_1110.json (filename is a capture identifier; this heading carries the actual completion timestamp). No source anchors were followed and no PDF scrolling was performed. A separate temporary review tab was used after the user's tab did not switch reliably. No Overleaf source, replies, resolutions, or tracked-change decisions were changed. Only this local report and evidence were written; no LaTeX build was required. No commits or pushes.

---

# Monitoring update — 23 October 2026, 10:15 a.m. Pacific

**New advisor activity:** three new Chapter 5 comments and fourteen new pending Chapter 5 edits since the 9:40 a.m. snapshot. Current Overview: **43 open advisor threads, 41 pending advisor edits, and 7 pending Ben edits** (91 entries). The older inventories below remain historical snapshots.

| File | Open advisor threads | Advisor edits | Ben edits |
|---|---:|---:|---:|
| 10-introduction.tex | 6 | 7 | 0 |
| 15-background.tex | 16 | 20 | 2 |
| 30-torque.tex | 0 | 0 | 3 |
| 50-database.tex | 12 | 14 | 0 |
| 60-sns.tex | 9 | 0 | 2 |

The three new advisor comments in 50-database.tex are “This is unclear” (9:51 a.m., position 1637), “Unclear why there would be an 'archived' record.” (9:52 a.m., position 1718), and “organization is still unclear” (9:52 a.m., position 1816). The full changes are transcribed below. Positions are review-card offsets, not source lines. No comment anchors were followed during this check, to avoid editor/PDF scrolling; exact marked passages were not recaptured.

Suggested response strategy: explain the resource’s contents and table relationships before its implementation details. Define what an archived record is, why it is retained, and how it differs from a retained paper, using the actual inclusion/exclusion rules. If accurate for the current app, possible wording is: “Records excluded from the retained collection remain in the archive so that their exclusion decisions can be inspected.” Verify that purpose against the curation record before using it. Do not invent archive behavior or treat the advisor’s proposed table naming as verified implementation.

Ben also replied to Introduction’s “double check” at 10:10 a.m.: “I didn't redesign the linkages”. This clarification narrows NEW-01 above: distinguish actuator attachment/routing redesign from mechanical-linkage design, and preserve Steele’s linkage attribution. Five newly observed pending Ben edits concern the torque figure and SNS verification; they are not advisor edits.

Evidence: Notes/advisor_review_20261023/entries_20261023_1011.json. This capture used the existing alternate Overleaf tab’s Overview, expanded text only, and did not follow source anchors, scroll the PDF, edit source, post replies, resolve threads, or accept/reject changes. The full prior report is preserved below.

## New or changed items since the prior capture

### CHECK-1011-01 · 10-introduction.tex · comment · position 8404

~~~text
ajh26
23 October, 9:10 am
double check
bbolen83
23 October, 10:10 am
I didn't redesign the linkages
~~~

### CHECK-1011-02 · 30-torque.tex · edit · position 64747

~~~text
bbolen83
23 October, 9:48 am
Changed: addresses.
to ...addresses. Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of Section~\ref{sec:redesign_physical}: across its fifteen tested measurements from \ang{-119} to \ang{5}, the measured torque magnitude met or exceeded the human biceps-femoris short-head maximum at every point, with margins from \qty{0.42}{\newton\meter} at \ang{-105} to \qty{17.0}{\newton\meter} at \ang{-14}.
~~~

### CHECK-1011-03 · 30-torque.tex · edit · position 65172

~~~text
bbolen83
23 October, 9:50 am
Changed: \begin{center} \includegraphics[width=1\textwidth]{FlxGrp.eps} \end{center}
to \begin{center} \includegraphics[width=0.72\textwidth]{FlxGrp.eps} \par\smallskip \begin{minipage}[t]{0.30\textwidth}\vspace{0pt}\centering {\textbf{(C)}}\par\smallskip \includegraphics[width=\linewidth]{Flx20mm_meetsHuman.pdf} \end{minipage} \end{center}
~~~

### CHECK-1011-04 · 30-torque.tex · edit · position 66396

~~~text
bbolen83
23 October, 9:53 am
Added:  \textbf{(C)} The redesigned two-BPA flexor route
~~~

### CHECK-1011-05 · 50-database.tex · edit · position 576

~~~text
ajh26
23 October, 9:48 am
Changed: The knowledge gap was the lack of a
to When I began this work, there was no
~~~

### CHECK-1011-06 · 50-database.tex · edit · position 632

~~~text
ajh26
23 October, 9:48 am
Deleted: , available to this project,
~~~

### CHECK-1011-07 · 50-database.tex · edit · position 699

~~~text
ajh26
23 October, 9:49 am
Changed: support
to be used to design
~~~

### CHECK-1011-08 · 50-database.tex · edit · position 1101

~~~text
ajh26
23 October, 9:49 am
Added:  the
~~~

### CHECK-1011-09 · 50-database.tex · comment · position 1637

~~~text
ajh26
23 October, 9:51 am
This is unclear
~~~

### CHECK-1011-10 · 50-database.tex · comment · position 1718

~~~text
ajh26
23 October, 9:52 am
Unclear why there would be an 'archived' record.
~~~

### CHECK-1011-11 · 50-database.tex · comment · position 1816

~~~text
ajh26
23 October, 9:52 am
organization is still unclear
~~~

### CHECK-1011-12 · 50-database.tex · edit · position 1898

~~~text
ajh26
23 October, 9:54 am
Changed: in
to is contained in
~~~

### CHECK-1011-13 · 50-database.tex · edit · position 1922

~~~text
ajh26
23 October, 9:54 am
Added: , an online database program. It
~~~

### CHECK-1011-14 · 50-database.tex · edit · position 1971

~~~text
ajh26
23 October, 9:57 am
Added: tables for
~~~

### CHECK-1011-15 · 50-database.tex · edit · position 1998

~~~text
ajh26
23 October, 9:56 am
Added:  Pathways
~~~

### CHECK-1011-16 · 50-database.tex · edit · position 2017

~~~text
ajh26
23 October, 9:56 am
Added: and
~~~

### CHECK-1011-17 · 50-database.tex · edit · position 2034

~~~text
ajh26
23 October, 9:57 am
Deleted: , and Recorders tables
~~~

### CHECK-1011-18 · 50-database.tex · edit · position 2438

~~~text
ajh26
23 October, 9:57 am
Changed: The other tables organize review literature and recording methods
to Review papers contain a short summary
~~~

### CHECK-1011-19 · 50-database.tex · edit · position 2648

~~~text
ajh26
23 October, 10:00 am
Deleted:  I curate the linked Airtable base and rebuild the read-only application from its exported records.
~~~

### CHECK-1011-20 · 50-database.tex · edit · position 3160

~~~text
ajh26
23 October, 10:01 am
Added: This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication.
~~~

### CHECK-1011-21 · 50-database.tex · edit · position 3708

~~~text
ajh26
23 October, 10:01 am
Deleted: This was a curated research collection, rather than a registered systematic review or an exhaustive search of every relevant publication.
~~~

### CHECK-1011-22 · 60-sns.tex · edit · position 6785

~~~text
bbolen83
23 October, 9:43 am
Changed: & Mechanism import achieved & Standalone reflex demonstrations do not establish control of imported models
to ...& Mechanism import achieved & The imported lower-body mechanism's ankle joints remain rigid, and standalone reflex demonstrations do not establish control of imported models
~~~

### CHECK-1011-23 · 60-sns.tex · edit · position 12681

~~~text
bbolen83
23 October, 9:44 am
Changed: The MuJoCo--Simulink bridge was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} closed-loop run. Commands and measurements passed through the interface. However, under tonic drive alone the Simulink realization reached a fixed point while the NumPy reference oscillated. The network-level discrepancy remains a limit on transferring tuned controller parameters between implementations.
to The MuJoCo--Simulink bridge, built on the official MathWorks blockset that embeds the MuJoCo model inside Simulink as the mechanical model, was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} open-loop run, without afferent return; the closed-loop reflex experiment exercised commands and sensory measurements through the same interface. However, under tonic drive alone the Simulink realization reached a fixed point while the NumPy reference oscillated. The network-level discrepancy remains a limit on transferring tuned controller parameters between implementations.
~~~


---

# Current advisor review — 23 October 2026

Read-only snapshot of [Bolen-Dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb), captured 10/23/2026, 9:40:23 AM Pacific. **No Overleaf source, comments, replies, resolutions, or tracked-change decisions were changed.** Suggestions are local proposals for Ben.

## Current inventory, superseding historical counts below

**40 open advisor threads (42 advisor comment messages), 27 pending advisor edits, and 2 pending Ben edits.** All 69 rendered Overview entries were expanded and captured. Threads with Ben replies remain open; a reply does not establish advisor acceptance or independently verified implementation. The resolved archive was not re-audited.

| File | Open advisor threads | Advisor edits | Ben edits |
|---|---:|---:|---:|
| 10-introduction.tex | 6 | 7 | 0 |
| 15-background.tex | 16 | 20 | 2 |
| 50-database.tex | 9 | 0 | 0 |
| 60-sns.tex | 9 | 0 | 0 |

Six threads are new this morning, one in Introduction and five in Background. All 27 pending advisor edits are dated this morning. Older inventories and revision logs are preserved verbatim under Historical report, not presented as current status. Missing older threads are not classified as resolved without consulting the archive.

Evidence: Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023/ contains the full expanded overview, DOM-derived entries and thread identifiers, source windows for the six new comments, and the pre-pass report. Positions below are live review-card character offsets, not source line numbers. Many older anchors now share a position after chapter replacement. Clicking new cards exposed contextual source windows but no exact text selection; suggested targets below are explicitly interpretations.

## Suggested edits for this morning’s comments

### NEW-01 · Introduction · 9:10 am · “double check” · position 8404

The clicked source window starts with the master’s-project paragraph, placement optimization, Steele’s knee linkage, and Muscle Mutt/DoggyDeux, followed by Research Objectives. The exact highlighted phrase was unavailable. Check which contributions belong to the master’s project versus this dissertation, and what the cited quadruped actually demonstrated.

Suggested narrower wording, after checking the cited work: “My master’s work compared the biomimetic robot leg with the Gait2392 model and developed MATLAB tools for evaluating muscle lengths, moment arms, and isometric joint torques. This dissertation extends that framework through experimental actuator characterization, corrected torque predictions, and actuator-route redesign.” Retain citations and separate attribution of collaborators’ hardware. The underlying publications were not independently checked in this pass.

### NEW-02 · Background · 8:39 am · position 27746

Advisor: “Chapter 5 attempts to collect and organize all known studies on sensory feedback.”

The source window contains Sensory Feedback in Locomotion. Likely intent: connect this overview to the database chapter. “All known” describes an aim, not demonstrated exhaustive coverage.

Suggested addition: “Chapter 5 presents the Sensory Afferent Database, which collects and organizes studies of sensory feedback relevant to locomotor control. Its coverage and the limits of the collection are described there.” Use the existing chapter label in LaTeX.

### NEW-03 · Background · 8:42 am · “huh?” · position 33161

The clicked window includes the Neuromechanical Simulation paragraph ending “the robustness benefit of each depends on which muscle carries it.” That wording obscures both the strategies and what “carries” means. This is a likely target, not a recovered exact highlight.

Suggested replacement for that clause: “They found that the contribution of feed-forward activation and reflex feedback to walking robustness varied among muscles.” Preserve the existing Haeufle citation, and check its result before using more specific wording. This is a clarity proposal, not independent literature verification.

### NEW-04 and NEW-05 · Background · 8:56 am · two “huh?” threads · positions 39194 and 39226

Both clicks expose the ethical-case opening: a robot as a replacement method, every simulated intervention avoiding a living subject, duplicated models as statistically identical subjects, and repeatability “no breeding colony can offer.” The two threads remain distinct despite identical displayed source windows. Exact marked ranges were unavailable. Ben’s pending repeatability wording is recorded separately in the ledger.

Likely problem: reproducible initial conditions do not establish independent biological replication, reduced animal use, or biological validation. The Three Rs connection needs conditional wording.

Suggested replacement for the replacement/repeatability sentences: “Robotic and computational models allow proposed neural interventions to be tested under controlled conditions before animal experiments are considered. The same model and initial conditions can be reused across interventions, helping isolate the effect of a changed pathway. These tests may inform the design of biological experiments, but their biological relevance depends on validation against experimental observations.” Retain the existing Three Rs citations with their definition. Do not claim an observed reduction in animal use without evidence.

### NEW-06 · Background · 8:56 am · “Non seq” · position 39294

The clicked window also contains the ethical-case paragraph and transition to the decerebrate-cat figure. Nearby, “The failure modes of these experiments damage hardware, not tissue” shifts abruptly from repeatability to hardware damage and does not describe purely computational interventions. Exact highlighting was unavailable.

Suggested transition, if this is the intended target: “These models permit interventions on synthetic neural circuits without performing the corresponding intervention on a living animal.” Follow with a factual account of the figure’s scientific relevance. Any removal or relocation of the figure remains a proposal.

## Pending-edit issues before acceptance

- Introduction, 9:07 am: the proposed replacement contains “into a complete into a complete” and “rigurous.” Possible cleanup: “they had not been combined into a complete human neuromechanical model for rigorous testing of different hypotheses.” Verify the novelty claim and distinguish the human model from the BPA robot model.
- Background, 8:47 am: the displayed resulting source contains a standalone “SNS Toolbox” citation without a sentence. Possible wording: “SNS-Toolbox simulates the neural controller and exchanges signals with the mechanical simulation.” Check the interface attribution against Chapter 6.
- Background, 8:49–8:50 am: the current passage claims thousands of parameter combinations overnight, hundreds of gait cycles per condition, unmatched biological throughput, and statistically powered experiments. These need evidence. Possible qualitative wording: “The coupled model supports repeated parameter, stimulation, and pathway-removal tests under controlled simulation conditions.” Throughput alone does not establish statistical power. No benchmark numbers were changed.
- Four Background additions at 8:47 am contain only newlines; their exact content survives in the captured DOM.

## Complete current open review ledger

Each entry retains visible author, date, and full text. ajh26 is the advisor; bbolen83 is Ben. Suggestions above are not posted replies. LIVE identifiers belong to this snapshot; older INTRO/BG/DB/SNS identifiers remain in the preserved history.

### LIVE-001 · 10-introduction.tex · open thread · position 1897

~~~text
ajh26
14 October, 10:59 am
This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too.
bbolen83
23 October, 6:20 am
I revised the opening to state the sensory-evidence and model-integration gaps directly and explain how the dissertation addresses them. The counts are now in to Chapter 5.
~~~

### LIVE-002 · 10-introduction.tex · pending advisor edit · position 3304

~~~text
ajh26
23 October, 9:02 am
Changed: repair the discrepancy
to fill in the gaps
~~~

### LIVE-003 · 10-introduction.tex · pending advisor edit · position 3735

~~~text
ajh26
23 October, 9:03 am
Changed: rational
to accurate
~~~

### LIVE-004 · 10-introduction.tex · pending advisor edit · position 4873

~~~text
ajh26
23 October, 9:04 am
Changed: N
to Attempts to implement n
~~~

### LIVE-005 · 10-introduction.tex · pending advisor edit · position 4910

~~~text
ajh26
23 October, 9:04 am
Added: models
~~~

### LIVE-006 · 10-introduction.tex · pending advisor edit · position 5282

~~~text
ajh26
23 October, 9:06 am
Added: , and was difficult to sift through for new researchers to the field
~~~

### LIVE-007 · 10-introduction.tex · pending advisor edit · position 5379

~~~text
ajh26
23 October, 9:06 am
Changed: tested
to validated
~~~

### LIVE-008 · 10-introduction.tex · pending advisor edit · position 5683

~~~text
ajh26
23 October, 9:07 am
Changed: their availability did not establish controller behavior after model conversion, actuator substitution, or the addition of contact feedback
to they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses
~~~

### LIVE-009 · 10-introduction.tex · open thread · position 8404

~~~text
ajh26
23 October, 9:10 am
double check
~~~

### LIVE-010 · 10-introduction.tex · open thread · position 12469

~~~text
ajh26
22 October, 9:40 am
expand
bbolen83
23 October, 6:27 am
I expanded the explanation of the gap and its connection to controller design, rather than retaining the undefined label.
~~~

### LIVE-011 · 10-introduction.tex · open thread · position 12469

~~~text
ajh26
22 October, 9:41 am
You should have a good idea of how many now that you have gathered them.
bbolen83
23 October, 6:21 am
Chapter 5 now reports the audited corpus size and retained-paper coverage, with explicit limits on what those counts establish.
~~~

### LIVE-012 · 10-introduction.tex · open thread · position 12469

~~~text
ajh26
22 October, 9:44 am
use mirror language here and for your chapter order. Your chapters should be in the order you present the gaps/you should present the gaps in the order your chapters are.
bbolen83
23 October, 6:21 am
I aligned the gap statements, contribution order, and chapter overview so they follow the same sequence.
~~~

### LIVE-013 · 10-introduction.tex · open thread · position 12469

~~~text
ajh26
22 October, 9:49 am
It is not clear how this work connects to the Biomimetic Humanoid Project and how your work builds on it.
bbolen83
23 October, 6:21 am
I rewrote this transition to explain how the work extends the lab’s biomimetic humanoid project through actuator characterization, corrected joint-torque predictions, and neural/mechanical model interfaces.
~~~

### LIVE-014 · 15-background.tex · open thread · position 6563

~~~text
ajh26
22 October, 10:42 am
Motor center is a little off
bbolen83
23 October, 6:21 am
I checked the replacement motor schematic: the rotor and stator now share a center.
~~~

### LIVE-015 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 10:45 am
This figure is confusing and probably not needed,
bbolen83
23 October, 6:21 am
I removed the confusing graphical abstract and retained the figures that support specific comparisons in the text.
~~~

### LIVE-016 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 10:50 am
not on your collected data though, right? Be more specific.
bbolen83
23 October, 6:21 am
These are OpenSim model-derived moment arms, not measurements from my robot tests. I made that distinction explicit.
~~~

### LIVE-017 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 10:50 am
huh?
ajh26
22 October, 10:50 am
definitions are stressed?
bbolen83
23 October, 6:27 am
I replaced the unclear terminology with a signed fixed-axis moment-arm definition and a separate explanation of OpenSim’s coordinate-length sensitivities.
~~~

### LIVE-018 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 10:52 am
huh?
bbolen83
23 October, 6:27 am
I removed the repository digression and rewrote this explanation around the mechanics and the Gait2392 model.
~~~

### LIVE-019 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 10:53 am
I'm now a little lost...which method does 2392 use?
bbolen83
23 October, 6:21 am
I clarified the Gait2392 knee kinematics and distinguished OpenSim’s coordinate-length sensitivity calculation from the fixed-axis geometric moment-arm expression.
~~~

### LIVE-020 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 10:53 am
I assume the 23 corresponds to the DoF. What about the 92?
bbolen83
23 October, 6:21 am
I now define Gait2392 as a model with 23 degrees of freedom and 92 musculotendon actuators.
~~~

### LIVE-021 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 11:49 am
non seq
bbolen83
23 October, 6:27 am
I shortened the geometry explanation and stated what the accompanying figures compare.
~~~

### LIVE-022 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 11:56 am
I feel like this is jumping all over the place
bbolen83
23 October, 6:27 am
I reorganized this passage around the prior model, the human/robot comparison, and the subsequent placement optimization, removing the file-saving chronology.
~~~

### LIVE-023 · 15-background.tex · open thread · position 26009

~~~text
ajh26
22 October, 11:57 am
what comparison was made? You don't say
bbolen83
23 October, 6:21 am
I now state the human-versus-robot comparison and the role of the subsequent attachment-placement optimization explicitly.
~~~

### LIVE-024 · 15-background.tex · open thread · position 27746

~~~text
ajh26
23 October, 8:39 am
Chapter 5 attempts to collect and organize all known studies on sensory feedback.
~~~

### LIVE-025 · 15-background.tex · open thread · position 30753

~~~text
ajh26
22 October, 12:13 pm
non seq
bbolen83
23 October, 6:21 am
I removed the misplaced dissertation-controller sentence and kept this section focused on computational neuroscience generally. CPG-specific material is in the CPG discussion.
~~~

### LIVE-026 · 15-background.tex · open thread · position 33161

~~~text
ajh26
23 October, 8:42 am
huh?
~~~

### LIVE-027 · 15-background.tex · pending advisor edit · position 35491

~~~text
ajh26
23 October, 8:43 am
Changed: ing
to s
~~~

### LIVE-028 · 15-background.tex · pending advisor edit · position 35597

~~~text
ajh26
23 October, 8:44 am
Added: elaborate
~~~

### LIVE-029 · 15-background.tex · pending advisor edit · position 35631

~~~text
ajh26
23 October, 8:44 am
Changed: too elaborate to integrate by hand, joined to
to and complex
~~~

### LIVE-030 · 15-background.tex · pending advisor edit · position 35657

~~~text
ajh26
23 October, 8:44 am
Changed: too large to wire by hand.
to .
~~~

### LIVE-031 · 15-background.tex · pending advisor edit · position 35694

~~~text
ajh26
23 October, 8:45 am
Deleted: , and the ones adjacent to it in the same ecosystem,
~~~

### LIVE-032 · 15-background.tex · pending advisor edit · position 36193

~~~text
ajh26
23 October, 8:46 am
Deleted: , and near the ends of a joint's range the routed paths of sparsely wrapped muscles can bunch up or pass through bone, an interference documented within Gait2392-derived routings themselves
~~~

### LIVE-033 · 15-background.tex · pending advisor edit · position 36427

~~~text
ajh26
23 October, 8:47 am
Added: 
~~~

### LIVE-034 · 15-background.tex · pending advisor edit · position 36616

~~~text
ajh26
23 October, 8:47 am
Added: 
~~~

### LIVE-035 · 15-background.tex · pending advisor edit · position 36852

~~~text
ajh26
23 October, 8:47 am
Added: 
~~~

### LIVE-036 · 15-background.tex · pending advisor edit · position 37054

~~~text
ajh26
23 October, 8:47 am
Added: 
~~~

### LIVE-037 · 15-background.tex · pending advisor edit · position 37279

~~~text
ajh26
23 October, 8:47 am
Changed: Combined with the SNS-Toolbox's MuJoCo backend \citep{nourse_sns_2023},
to SNS Toolbox \citep{nourse_sns_2023}
~~~

### LIVE-038 · 15-background.tex · pending advisor edit · position 37340

~~~text
ajh26
23 October, 8:48 am
Changed: convert the laboratory's
to develop a
~~~

### LIVE-039 · 15-background.tex · pending advisor edit · position 37350

~~~text
ajh26
23 October, 8:49 am
Changed: modified
to robotic
~~~

### LIVE-040 · 15-background.tex · pending advisor edit · position 37375

~~~text
ajh26
23 October, 8:49 am
Added: , convert it
~~~

### LIVE-041 · 15-background.tex · pending advisor edit · position 37397

~~~text
ajh26
23 October, 8:49 am
Added: ,
~~~

### LIVE-042 · 15-background.tex · pending advisor edit · position 37504

~~~text
ajh26
23 October, 8:49 am
Changed: Simulation speed is the practical enabler
to This new neuromechanical model can simulate very quickly
~~~

### LIVE-043 · 15-background.tex · pending advisor edit · position 37689

~~~text
ajh26
23 October, 8:50 am
Changed: ,
to .
~~~

### LIVE-044 · 15-background.tex · pending advisor edit · position 37691

~~~text
ajh26
23 October, 8:50 am
Changed: t
to T
~~~

### LIVE-045 · 15-background.tex · pending advisor edit · position 37754

~~~text
ajh26
23 October, 8:50 am
Changed: ,
to .
~~~

### LIVE-046 · 15-background.tex · pending advisor edit · position 37756

~~~text
ajh26
23 October, 8:50 am
Changed: and that
to This new model
~~~

### LIVE-047 · 15-background.tex · pending Ben edit · position 39166

~~~text
bbolen83
23 October, 9:27 am
Changed: form of reduction
to n experimental repeatability
~~~

### LIVE-048 · 15-background.tex · open thread · position 39194

~~~text
ajh26
23 October, 8:56 am
huh?
~~~

### LIVE-049 · 15-background.tex · pending Ben edit · position 39225

~~~text
bbolen83
23 October, 9:28 am
Changed: Refinement follows as well: t
to T
~~~

### LIVE-050 · 15-background.tex · open thread · position 39226

~~~text
ajh26
23 October, 8:56 am
huh?
~~~

### LIVE-051 · 15-background.tex · open thread · position 39294

~~~text
ajh26
23 October, 8:56 am
Non seq
~~~

### LIVE-052 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:14 am
There are still lots of unknowns
bbolen83
23 October, 6:27 am
I removed the suggestion that the physiology is fully specified. The chapter now states the uneven coverage and the interpretation still needed before using a pathway in a controller.
~~~

### LIVE-053 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:16 am
After this sentence, define the gap in knowledge and what you did to fill it.
bbolen83
23 October, 6:21 am
I put the knowledge gap immediately after introducing SADb: the project lacked a structured, source-linked account of afferent measurements and pathways for selecting feedback rules. The next sentences explain how SADb addresses it.
~~~

### LIVE-054 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:20 am
I thought it was presented 2x, once without the AI help, and the other with.
bbolen83
23 October, 6:21 am
Both presentations are now identified and cited separately: the database in 2023 and AI-assisted triage in 2024.
~~~

### LIVE-055 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:25 am
I thought you migrated it to a different system
bbolen83
23 October, 6:21 am
I clarified the two roles: Airtable remains the editable linked working base, and the offline browser application is its read-only distribution. The app now carries all 1,393 pre-prune records; pivot and citation views remain paper-only.
~~~

### LIVE-056 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:33 am
You don't have to follow this exactly, but you should generally put this in the Chochrane form. https://www.cochrane.org/sites/default/files/uploads/PDFs/guide_to_the_contents_of_a_cochrane_methodology_protocol_and_review.pdf
bbolen83
23 October, 6:21 am
I reorganized Curation Methods around sources, identity matching, inclusion criteria, source-text checking, and my decision authority. I also state that this was a curated research collection rather than an exhaustive or registered systematic review.
~~~

### LIVE-057 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:37 am
what are notes?
bbolen83
23 October, 6:21 am
I now define curation notes as short, source-grounded summaries of the mechanism, experimental conditions, and relevance to feedback design.
~~~

### LIVE-058 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:39 am
You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized?
bbolen83
23 October, 6:21 am
The chapter now opens with the finished resource, its purpose, record types, organization, and a screenshot. Implementation history and supporting counts no longer precede that overview.
~~~

### LIVE-059 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:40 am
approved by who?
bbolen83
23 October, 6:21 am
I made my responsibility explicit: I decide inclusion, classification, and deletion after checking the source.
~~~

### LIVE-060 · 50-database.tex · open thread · position 7769

~~~text
ajh26
21 October, 11:40 am
how was this decided?
ajh26
21 October, 11:41 am
how were they identified?
bbolen83
23 October, 6:21 am
I added the relevance criteria and described DOI/title matching and source-text verification. Missing source text is flagged rather than reconstructed from recollection.
~~~

### LIVE-061 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 11:48 am
So many details and I don't understand what the goal is. What is the gap in knowledge? What did you do to fill it?
bbolen83
21 October, 2:09 pm
Gaps: Gait2392_simbody with CPG in Mujoco-SNS Toolbox
- Gait2392_robotbody (BPAs and BPA routing) with CPG
- Block-programming on a fast physics simulator (Animatlab no longer supported)
- Robot solidworks assembly with BPAs from Matlab and a neural controller
- Li walker with Deng neural controller and foot contact feedback.
bbolen83
23 October, 6:21 am
I rebuilt the opening around the Li/Deng AnimatLab walker and my interconnected left/right RG extension. It then identifies the five missing model/controller combinations in my earlier reply and explains the tools and verification needed to investigate them. Detailed supporting results are in Appendix D.
~~~

### LIVE-062 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 12:45 pm
huh?
bbolen83
23 October, 6:27 am
I removed the build-order shorthand. The opening now connects the bilateral AnimatLab construction to the missing interfaces and the verification each new tool provides.
~~~

### LIVE-063 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 12:49 pm
why would we need several?
bbolen83
21 October, 2:04 pm
Well, Mujoco doesn't have as good of a physics engine as HYFYDY, for example, and is limited to rigid tendons. Python also is harder to use because it is code without a GUI. Simscape multibody allows for importing parts from solidworks. Some sort of block-programming language like Simulink is a natural choice for building our style of connectome. Mujoco-SNS Toolbox allows for real-time control and computer in the loop commands.
bbolen83
23 October, 6:21 am
I now explain the complementary needs: scripted neural/mechanical coupling, graphical neural construction, and SolidWorks/MATLAB integration. The capability table separates demonstrated interfaces from unfinished combinations, and the text states the relevant simulator limitations.
~~~

### LIVE-064 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 12:53 pm
This is all super awkward. You need to identify the gap in knowledge and how you fill that gap in knowledge. Much of what you wrote here should be saved for discussion, where you talk about the implications of the work.
bbolen83
21 October, 2:00 pm
I thought that was for research hypothesis not for "here are some tools we made".
bbolen83
23 October, 6:21 am
I put the knowledge gaps before the new tools and kept the main results to representative verification. Broader sweeps and supported-walker results are in Appendix D; future integration belongs in the Discussion chapter.
~~~

### LIVE-065 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 12:59 pm
I don't understand the connection.
bbolen83
23 October, 6:27 am
I rewrote the transition to connect the identified modeling gaps to the specific neural/mechanical interfaces developed.
~~~

### LIVE-066 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 1:00 pm
Overlapping text in figure
bbolen83
23 October, 6:21 am
I replaced the pipeline schematic with the corrected vector figure, with separated columns and labels.
~~~

### LIVE-067 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 1:02 pm
huh?
bbolen83
23 October, 6:27 am
I removed the unexplained inventory and retained the interface’s purpose, its demonstrated check, and its limitation.
~~~

### LIVE-068 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 1:08 pm
Wait, that is all you are saying about the other tools? It is hard for the reader to understand the limitations and gaps without a deeper description of them.
bbolen83
23 October, 6:21 am
I expanded the explanation of the environments’ capabilities and limitations and added a table showing which interfaces are demonstrated and which combinations remain incomplete.
~~~

### LIVE-069 · 60-sns.tex · open thread · position 14367

~~~text
ajh26
21 October, 1:11 pm
Did you compare and discuss against the libraries myself and Fletcher made here?
bbolen83
21 October, 2:00 pm
No. You've never encouraged us to use Simulink.
bbolen83
23 October, 6:21 am
I removed the Hunt/Young library comparison because those citations do not establish a shared neuromechanical Simulink library. The section now describes the implemented blocks and matched-circuit verification directly. A comparison with the specific libraries you mean still requires checking those libraries.
~~~

## Verification and scope

All 69 rendered Overview entries are included: 40 threads, 27 advisor edits, and 2 Ben edits. No “show more” controls remained after expansion. Six new-comment source windows were captured. Only this report and its local evidence folder were written in this pass. No thesis source or bibliography was edited, so a LaTeX build was not applicable. No commits, pushes, uploads, replies, accept/reject actions, or resolutions were performed. Existing changes from other work were preserved.

---

# Historical report and earlier revision records (preserved verbatim)

The following inventory, source line numbers, completion statements, and uploads describe earlier passes, not this read-only review.

# Advisor review and revision handoff

Updated from the live **Bolen-Dissertation** Overleaf Review → Overview, final capture **22 October 2026, 12:22:09 Pacific**. This was a review of advisor contents and their anchors, not another full-document quality review. The advisor is actively reviewing; this is a snapshot, not a claim that no later comments exist.

## Current inventory

**73 open advisor threads (76 advisor messages) and 51 pending advisor tracked edits.** Ben's replies are explicitly separated below. Resolved comments and Ben-only appendix/main entries are excluded from the advisor totals. The previous 36-thread/8-edit inventory is superseded and preserved in Notes/advisor_review_20261022/CHATGPT_REPORT_before_this_pass.md. Existing SNS-01..15, DB-01..26, INTRO-01, BG-01, and BG-02 identifiers are retained.

| File | Advisor threads | Pending advisor edits |
|---|---:|---:|
| 10-introduction.tex | 14 | 20 |
| 15-background.tex | 24 | 25 |
| 50-database.tex | 24 | 2 |
| 60-sns.tex | 11 | 4 |

Line numbers and context refer to the fresh pre-upload source snapshot in tmp/advisor_review_20261022/overleaf_before/chapters, except the later live background snapshot saved as Notes/advisor_review_20261022/background_final.tex. They do not refer to the differing local ProofFinal chapters or PDF line numbers. A long paragraph occupies one source line. Positions are the rendered review panel's character offsets. Exact highlighted text was captured from the editor where available; a context window is explicitly labeled when no highlight was returned. Thread IDs are retained in the JSON evidence. Source offsets can drift after concurrent edits. Tracked edits remain pending; recording them does not accept them.

## Ben's current controlling instructions

These supersede conflicting older campaign plans and the earlier ZCode prompt:

1. Reduce the dissertation, especially SADb and SNS. Quality and a clear argument take precedence over more text. Do not bury a cut by adding an equally long appendix.
2. SNS does **not** require Methods/Results/Discussion headings. Lead with Ben's PhD-time combination of the Deng CPG and Li walker in AnimatLab: stepping-like movement in air, but no walking achieved by hand tuning. Use that limitation to motivate the knowledge/capability gaps in Ben's replies below. Follow with brief results of new software tools and concrete verification showing what works and what gap is closed. Future work belongs mainly in Chapter 7.
3. Preserve the advisor's active review. Supply narrowly targeted changes wherever possible; prefer a few paragraphs or one section to replacement of entire chapters. Log moves/deletions, retain labels and citations, and do not overwrite Overleaf from stale local chapters.
4. Bring Chapter 4 into line with verified recent torque-paper work from GitHub, distinguishing new measurements, uncorrected estimates, completed reductions, and predicted results.
5. Preserve Chapter 3's isometric-force figures and MATLAB figures. Reports and prompts stay local; **never upload CHATGPT_REPORT or handoff/administrative files to Overleaf.**

## Advisor comments and tracked changes

## 10-introduction.tex

### INTRO-01 · line 8 · comment · 14 October, 10:59 am

**Highlighted source:**

```latex
Artificial
```

**Advisor — 14 October, 10:59 am:**

> This paragraph does a great job identifying one of the gaps of knowledge you worked to fill, but there are a few others you can worked on too.


### INTRO-02 · line 12 · tracked edit · 22 October, 9:38 am

**Context near the review anchor (not an exact highlight):**

```latex
rque data, so even the nominal torque envelope of the robot leg rested on unverified calculations.

There are additional gaps concerning the development of more detailed neural control models of locomotion and balance. No simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kin
```

**Pending advisor edit:**

```text
Changed: The remaining
to There are additional
```


### INTRO-03 · line 12 · tracked edit · 22 October, 9:38 am

**Context near the review anchor (not an exact highlight):**

```latex
rque envelope of the robot leg rested on unverified calculations.

There are additional gaps concerning the development of more detailed neural control models of locomotion and balance. No simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kind and a biologically grounded neu
```

**Pending advisor edit:**

```text
Added: ing
```


### INTRO-04 · line 12 · tracked edit · 22 October, 9:36 am

**Context near the review anchor (not an exact highlight):**

```latex
elope of the robot leg rested on unverified calculations.

There are additional gaps concerning the development of more detailed neural control models of locomotion and balance. No simulation pipeline existed for closing the loop between a musculoskeletal robot model of this kind and a biologically grounded neural cont
```

**Pending advisor edit:**

```text
Changed: controller
to development of more detailed neural control models of locomotion and balance
```


### INTRO-05 · line 12 · tracked edit · 22 October, 9:36 am

**Context near the review anchor (not an exact highlight):**

```latex
 existed for closing the loop between a musculoskeletal robot model of this kind and a biologically grounded neural controller. The neuromechanical modeling literature provides the ingredients, including two-level spinal central pattern generator architectures \citep{rybak_modelling_2006, mccrea_organization_2008}, neu
```

**Pending advisor edit:**

```text
Changed: realistic
to grounded
```


### INTRO-06 · line 12 · tracked edit · 22 October, 9:37 am

**Context near the review anchor (not an exact highlight):**

```latex
}, neuromechanical simulation frameworks \citep{markin_afferent_2010}, and synthetic nervous system design tools \citep{nourse_sns_2023}, but it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing. 

Wiring realistic feedback into such a controller exposes a fourth gap, this on
```

**Pending advisor edit:**

```text
Changed: tooling
to design tools
```


### INTRO-07 · line 12 · comment · 22 October, 9:40 am

**Highlighted source:**

```latex
but it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing.
```

**Advisor — 22 October, 9:40 am:**

> expand


### INTRO-08 · line 12 · comment · 22 October, 9:37 am

**Highlighted source:**

```latex
their
```

**Advisor — 22 October, 9:37 am:**

> undefined


### INTRO-09 · line 12 · tracked edit · 22 October, 9:41 am

**Context near the review anchor (not an exact highlight):**

```latex
t it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing. 

Wiring realistic feedback into such a controller exposes a fourth gap, this one in the literature rather than in methods: the physiological properties of the group Ia, group Ib, and group II afferents that should infor
```

**Pending advisor edit:**

```text
Added: 
```


### INTRO-10 · line 14 · comment · 22 October, 9:42 am

**Context near the review anchor (not an exact highlight):**

```latex
it does not provide their assembly for a BPA-driven bipedal robot with human-like muscle routing. 

Wiring realistic feedback into such a controller exposes a fourth gap, this one in the literature rather than in methods: the physiological properties of the group Ia, group Ib, and group II afferents that should inform 
```

**Advisor — 22 October, 9:42 am:**

> Awkward


### INTRO-11 · line 14 · comment · 22 October, 9:41 am

**Highlighted source:**

```latex
hundreds if not thousands
```

**Advisor — 22 October, 9:41 am:**

> You should have a good idea of how many now that you have gathered them.


### INTRO-12 · line 16 · tracked edit · 22 October, 9:43 am

**Context near the review anchor (not an exact highlight):**

```latex
er or an animal researcher could draw them.

Four gaps thus stood between the existing robot design methods and a robot that can walk under neural control: no actuator force model valid at arbitrary resting length, no quantified correction for the torque lost through brackets, tendons, and wrapping, no closed-loop simu
```

**Pending advisor edit:**

```text
Added: methods
```


### INTRO-13 · line 16 · comment · 22 October, 9:44 am

**Highlighted source:**

```latex
no closed-loop simulation pipeline joining the robot's musculoskeletal model to a biologically realistic neural controller, and no consolidated resource for the afferent data that such a controller requires.
```

**Advisor — 22 October, 9:44 am:**

> use mirror language here and for your chapter order. Your chapters should be in the order you present the gaps/you should present the gaps in the order your chapters are.


### INTRO-14 · line 22 · comment · 22 October, 9:49 am

**Highlighted source:**

```latex
In parallel, a biomimetic four-bar knee linkage with a migrating instantaneous center of rotation was designed and validated by Steele \citep{steele_development_2017, steele_biomimetic_2018, steele_biomimetic_2018-1}, and a canine-inspired quadruped platform, Muscle Mutt, formerly called DoggyDeux, was developed with antagonistic BPA pairs and synthetic nervous system control \citep{scharzenberger_design_2019}. That platform extended an earlier AARL demonstration of central pattern generator control of a pneumatically actuated robot \citep{hunt_development_2017}.
```

**Advisor — 22 October, 9:49 am:**

> It is not clear how this work connects to the Biomimetic Humanoid Project and how your work builds on it.


### INTRO-15 · line 26 · comment · 22 October, 9:49 am

**Highlighted source:**

```latex
first two gaps
```

**Advisor — 22 October, 9:49 am:**

> undef

**Advisor — 22 October, 9:50 am:**

> You should name them specifically


### INTRO-16 · line 26 · comment · 22 October, 9:49 am

**Highlighted source:**

```latex
fourth
```

**Advisor — 22 October, 9:49 am:**

> undef


### INTRO-17 · line 26 · comment · 22 October, 9:50 am

**Highlighted source:**

```latex
third
```

**Advisor — 22 October, 9:50 am:**

> undef


### INTRO-18 · line 41 · tracked edit · 22 October, 9:54 am

**Context near the review anchor (not an exact highlight):**

```latex
are adapted from the author's publication of this work \citep{bolen_2026}.

Chapter~\ref{ch:torque} improves the actuator model at the joint. It describes the robot leg model, the pinned-joint and biomimetic knee test stands, and the equipment used to collect isometric torque data; formulates the hybrid torque calculat
```

**Pending advisor edit:**

```text
Changed: carries
to improves
```


### INTRO-19 · line 41 · tracked edit · 22 October, 9:54 am

**Context near the review anchor (not an exact highlight):**

```latex
's publication of this work \citep{bolen_2026}.

Chapter~\ref{ch:torque} improves the actuator model at the joint. It describes the robot leg model, the pinned-joint and biomimetic knee test stands, and the equipment used to collect isometric torque data; formulates the hybrid torque calculation method and the improved
```

**Pending advisor edit:**

```text
Deleted: s
```


### INTRO-20 · line 41 · tracked edit · 22 October, 9:54 am

**Context near the review anchor (not an exact highlight):**

```latex
s publication of this work \citep{bolen_2026}.

Chapter~\ref{ch:torque} improves the actuator model at the joint. It describes the robot leg model, the pinned-joint and biomimetic knee test stands, and the equipment used to collect isometric torque data; formulates the hybrid torque calculation method and the improved 
```

**Pending advisor edit:**

```text
Changed: to
to at
```


### INTRO-21 · line 41 · comment · 22 October, 9:55 am

**Highlighted source:**

```latex
Gait2392
```

**Advisor — 22 October, 9:55 am:**

> Too specific


### INTRO-22 · line 45 · tracked edit · 22 October, 9:57 am

**Context near the review anchor (not an exact highlight):**

```latex
on living subjects. It opens with the preliminary AnimatLab bipedal walker results, whose hand-tuned, two-layer CPG produced suspended stepping, with neither weight-bearing walking nor any role for afferent feedback in weight-bearing motion demonstrated; that outcome set the verification-first approach the developed to
```

**Pending advisor edit:**

```text
Added: ,
```


### INTRO-23 · line 45 · comment · 22 October, 9:58 am

**Highlighted source:**

```latex
nor any role for afferent feedback in weight-bearing motion demonstrated
```

**Advisor — 22 October, 9:58 am:**

> huh?


### INTRO-24 · line 45 · tracked edit · 22 October, 9:59 am

**Context near the review anchor (not an exact highlight):**

```latex
eedback in weight-bearing motion demonstrated; that outcome set the verification-first approach the developed tools use. The tools themselves are the OpenSim-to-MuJoCo conversion of the robot's musculoskeletal model with its synthetic BPA knee, the Simulink SNS library built to prototype reflex topologies, the spiking-
```

**Pending advisor edit:**

```text
Added: developed
```


### INTRO-25 · line 45 · tracked edit · 22 October, 9:59 am

**Context near the review anchor (not an exact highlight):**

```latex
t-bearing motion demonstrated; that outcome set the verification-first approach the developed tools use. The tools themselves are the OpenSim-to-MuJoCo conversion of the robot's musculoskeletal model with its synthetic BPA knee, the Simulink SNS library built to prototype reflex topologies, the spiking-neuron mirrors, 
```

**Pending advisor edit:**

```text
Changed: follow
to use
```


### INTRO-26 · line 45 · comment · 22 October, 10:00 am

**Highlighted source:**

```latex
the spiking-neuron mirrors
```

**Advisor — 22 October, 10:00 am:**

> huh?


### INTRO-27 · line 45 · tracked edit · 22 October, 10:02 am

**Context near the review anchor (not an exact highlight):**

```latex
, and the cross-platform verification protocol that ties the implementations together. Each tool is validated before the complete model is assembled. The stimulus and deletion experiments complete the validation process, and the chapter closes by discussing a five-stage training curriculum, the spiking-neuron mirrors, 
```

**Pending advisor edit:**

```text
Changed: demonstrated on the fundamental exercises that establish what it does and does not show: cross-engine reproduction of the contact-driven rhythm core, headless re-runs of the archived AnimatLab reference models, the BPA-actuated walker in MuJoCo, a gait-validation sweep against a locomotion reference library, and a comparison of the production walker variants
to validated before the complete model is assembled
```


### INTRO-28 · line 45 · tracked edit · 22 October, 10:02 am

**Context near the review anchor (not an exact highlight):**

```latex
alidated before the complete model is assembled. The stimulus and deletion experiments complete the validation process, and the chapter closes by discussing a five-stage training curriculum, the spiking-neuron mirrors, and the new tools the program leaves other researchers.

Chapter~\ref{ch:overall} discusses the work 
```

**Pending advisor edit:**

```text
Changed: demonstrations
to validation process
```


### INTRO-29 · line 45 · tracked edit · 22 October, 10:03 am

**Context near the review anchor (not an exact highlight):**

```latex
ulus and deletion experiments complete the validation process, and the chapter closes by discussing a five-stage training curriculum, the spiking-neuron mirrors, and the new tools the program leaves other researchers.

Chapter~\ref{ch:overall} discusses the work as a whole. It shows how the four chapters work together 
```

**Pending advisor edit:**

```text
Changed: the
to a
```


### INTRO-30 · line 45 · comment · 22 October, 10:03 am

**Highlighted source:**

```latex
spiking-neuron mirrors
```

**Advisor — 22 October, 10:03 am:**

> again, huh?


### INTRO-31 · line 47 · tracked edit · 22 October, 10:04 am

**Context near the review anchor (not an exact highlight):**

```latex
s the program leaves other researchers.

Chapter~\ref{ch:overall} discusses the work as a whole. It shows how the four chapters work together to develop improved biomimetic humanoid robots. Its also scopes the future work these results enable for other researchers, from BPA characterization and routing expansion throug
```

**Pending advisor edit:**

```text
Changed: opens with the outcome, that improved first-pass robot design is now achievable, and attributes it: the design factors that remove the rigid-body idealization, the improved BPA characterization, the tiered attachment and routing optimization, the sensory afferent database, and the new neuromechanical simulation tools
to shows how the four chapters work together to develop improved biomimetic humanoid robots
```


### INTRO-32 · line 47 · tracked edit · 22 October, 10:08 am

**Context near the review anchor (not an exact highlight):**

```latex
e. It shows how the four chapters work together to develop improved biomimetic humanoid robots. Its also scopes the future work these results enable for other researchers, from BPA characterization and routing expansion through database and simulation tools to real-time neural control of the physical robot. Chapter~\re
```

**Pending advisor edit:**

```text
Changed: remaining subsections
to also
```


### INTRO-33 · line 47 · tracked edit · 22 October, 10:08 am

**Context near the review anchor (not an exact highlight):**

```latex
s how the four chapters work together to develop improved biomimetic humanoid robots. Its also scopes the future work these results enable for other researchers, from BPA characterization and routing expansion through database and simulation tools to real-time neural control of the physical robot. Chapter~\ref{ch:concl
```

**Pending advisor edit:**

```text
Added: s
```


### INTRO-34 · line 47 · tracked edit · 22 October, 10:08 am

**Context near the review anchor (not an exact highlight):**

```latex
 other researchers, from BPA characterization and routing expansion through database and simulation tools to real-time neural control of the physical robot. Chapter~\ref{ch:conclusion} consolidates the principal contributions of the dissertation. Figure~\ref{fig:roadmap} summarizes the flow.

\begin{figure}[htbp]
\cent
```

**Pending advisor edit:**

```text
Changed: tooling
to tools
```


## 15-background.tex

### BG-03 · line 36 · comment · 22 October, 10:42 am

**Highlighted source:**

```latex
\begin{figure}[htbp]
```

**Advisor — 22 October, 10:42 am:**

> Motor center is a little off


### BG-01 · line 48 · tracked edit · 14 October, 1:31 pm

**Context near the review anchor (not an exact highlight):**

```latex
th relationship. Predicting their force as a function of length and pressure is therefore the first improvement of the torque tools developed in this dissertation, and the state of that prediction is the subject of the next subsection.

Figure~\ref{fig:festo_bpa_cad} shows the commercial actuator assembly used in the l
```

**Pending advisor edit:**

```text
Changed: requirement
to improvement
```


### BG-04 · line 100 · comment · 22 October, 10:45 am

**Highlighted source:**

```latex
\begin{figure}[htbp]
```

**Advisor — 22 October, 10:45 am:**

> This figure is confusing and probably not needed,


### BG-05 · line 111 · comment · 22 October, 10:53 am

**Highlighted source:**

```latex
Gait2392
```

**Advisor — 22 October, 10:53 am:**

> I assume the 23 corresponds to the DoF. What about the 92?


### BG-06 · line 113 · comment · 22 October, 10:49 am

**Highlighted source:**

```latex
the full-muscle rat hindlimb model of \citet{young_analyzing_2019} makes that adjacent-joint dependence a first-class result of its attachment-point technique.
```

**Advisor — 22 October, 10:49 am:**

> awkward


### BG-07 · line 113 · comment · 22 October, 10:49 am

**Highlighted source:**

```latex
first-class
```

**Advisor — 22 October, 10:49 am:**

> primary?


### BG-31 · line 113 · tracked edit · 22 October, 10:49 am

**Context near the review anchor (not an exact highlight):**

```latex
-point technique. OpenSim evaluates moment arms as muscle-length sensitivities of this kind, so the moment arms used throughout this dissertation are properties of the Gait2392 geometry rather than measured inputs \citep{delp_opensim_2007, seth_opensim_2018}.

The human knee is the joint where these definitions are str
```

**Pending advisor edit:**

```text
Added: moment
```


### BG-09 · line 113 · comment · 22 October, 10:50 am

**Highlighted source:**

```latex
used throughout this dissertation
```

**Advisor — 22 October, 10:50 am:**

> not on your collected data though, right? Be more specific.


### BG-10 · line 115 · comment · 22 October, 10:50 am

**Highlighted source:**

```latex
stressed
```

**Advisor — 22 October, 10:50 am:**

> huh?

**Advisor — 22 October, 10:50 am:**

> definitions are stressed?


### BG-11 · line 115 · comment · 22 October, 10:51 am

**Highlighted source:**

```latex
 Its instantaneous center of rotation migrates posteriorly with flexion, a behavior captured by four-bar linkage descriptions of the joint \citep{romero_four-bar_2020, steele_biomimetic_2018}, and \citet{romero_four-bar_2020} treat locating that moving center as a kinematic problem in its own right, solved with natural coordinates so that the answer stays cheap to evaluate when an optimization must recompute it many times
```

**Advisor — 22 October, 10:51 am:**

> run on. awkward


### BG-12 · line 115 · tracked edit · 22 October, 10:52 am

**Context near the review anchor (not an exact highlight):**

```latex
femoral joint give different moment arms for identical muscle paths, so the modeling choice must be determined carefully. The ISB joint coordinate system recommendations were written against exactly this definitional freedom, to keep reports of joint motion comparable across laboratories \citep{wu_isb_2002}. Throughout
```

**Pending advisor edit:**

```text
Changed: matched to the question asked of it
to determined carefully
```


### BG-13 · line 115 · comment · 22 October, 10:52 am

**Highlighted source:**

```latex
written against
```

**Advisor — 22 October, 10:52 am:**

> huh?


### BG-14 · line 115 · comment · 22 October, 10:53 am

**Highlighted source:**

```latex
judged
```

**Advisor — 22 October, 10:53 am:**

> I'm now a little lost...which method does 2392 use?


### BG-15 · line 119 · comment · 22 October, 10:55 am

**Highlighted source:**

```latex
It is a one-degree-of-freedom four-bar linkage whose crossed links reproduce the posterior migration of the human knee's instantaneous center of rotation (ICR) during flexion, so the rolling-and-sliding behavior of the human joint emerges from the link geometry itself rather than from cams or gear profiles; the design was calculated to reproduce an average sliding-to-rolling ratio of 0.079 over a \ang{151} range of motion \citep{steele_biomimetic_2018-1}.
```

**Advisor — 22 October, 10:55 am:**

> run on(s)


### BG-16 · line 119 · comment · 22 October, 10:56 am

**Highlighted source:**

```latex
which is exactly the behavior a torque prediction tool must capture
```

**Advisor — 22 October, 10:56 am:**

> non sequitor


### BG-17 · line 131 · tracked edit · 22 October, 10:58 am

**Context near the review anchor (not an exact highlight):**

```latex
mination_2019}. Hand placement of muscle origins and insertions was a painstaking and iterative task. This resulted in the creation of a placement optimization algorithm \citep{morrow_optimization_2020}, discussed in Chapter~\ref{ch:torque}. Figure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skelet
```

**Pending advisor edit:**

```text
Changed: , and the
to . This resulted in the creation of a
```


### BG-18 · line 131 · tracked edit · 22 October, 10:58 am

**Context near the review anchor (not an exact highlight):**

```latex
ions was a painstaking and iterative task. This resulted in the creation of a placement optimization algorithm \citep{morrow_optimization_2020}, discussed in Chapter~\ref{ch:torque}. Figure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed fo
```

**Pending advisor edit:**

```text
Changed: that replaced it is dissertation work
to algorithm
```


### BG-19 · line 131 · comment · 22 October, 11:49 am

**Highlighted source:**

```latex
igure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed for it
```

**Advisor — 22 October, 11:49 am:**

> non seq


### BG-20 · line 131 · tracked edit · 22 October, 10:59 am

**Context near the review anchor (not an exact highlight):**

```latex
_optimization_2020}, discussed in Chapter~\ref{ch:torque}. Figure~\ref{fig:aarl_opensim_comparison} shows the laboratory's bipedal skeleton beside the OpenSim muscle-routing reference developed for it.

\begin{figure}[htbp]
\centering
\begin{minipage}[c]{0.67\textwidth}
\centering
\textbf{(A)}\par\smallskip
\includegra
```

**Pending advisor edit:**

```text
Changed: places
to shows
```


### BG-21 · line 152 · tracked edit · 22 October, 11:50 am

**Context near the review anchor (not an exact highlight):**

```latex
n_determination_2019, morrow_optimization_2020}.}
\label{fig:aarl_opensim_comparison}
\end{figure}

An overview of the mechanics of that master's study are described below. Moment arms were computed geometrically from attachment points in the body frame, as the cross product of a joint-axis unit vector with the vector 
```

**Pending advisor edit:**

```text
Changed: T
to An overview of t
```


### BG-22 · line 152 · tracked edit · 22 October, 11:50 am

**Context near the review anchor (not an exact highlight):**

```latex
{fig:aarl_opensim_comparison}
\end{figure}

An overview of the mechanics of that master's study are described below. Moment arms were computed geometrically from attachment points in the body frame, as the cross product of a joint-axis unit vector with the vector from the joint center to the muscle's line of action, as
```

**Pending advisor edit:**

```text
Changed: worth recording because the dissertation's torque work reuses them
to described below
```


### BG-23 · line 152 · tracked edit · 22 October, 11:51 am

**Context near the review anchor (not an exact highlight):**

```latex
t of a joint-axis unit vector with the vector from the joint center to the muscle's line of action, as described in Section~\ref{sec:bg_biomechanics}. PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfe
```

**Pending advisor edit:**

```text
Changed: the same construction
to as
```


### BG-24 · line 152 · comment · 22 October, 11:56 am

**Highlighted source:**

```latex
PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfering with one another, and the result was saved as the gait2392\_robotbody model that the later chapters of this dissertation continue to refine \citep{bolen_bipedal_2019}. 
```

**Advisor — 22 October, 11:56 am:**

> Delete? 


### BG-25 · line 152 · comment · 22 October, 11:51 am

**Highlighted source:**

```latex
routings
```

**Advisor — 22 October, 11:51 am:**

> origin and insertion points?


### BG-26 · line 152 · comment · 22 October, 11:51 am

**Highlighted source:**

```latex
 moved to respect Festo BPA contraction limits
```

**Advisor — 22 October, 11:51 am:**

> huh?


### BG-27 · line 152 · comment · 22 October, 11:52 am

**Highlighted source:**

```latex
and the result was saved
```

**Advisor — 22 October, 11:52 am:**

> I don't need to be told that you saved a file.


### BG-28 · line 152 · comment · 22 October, 11:56 am

**Highlighted source:**

```latex
PAM routings were then placed by hand in OpenSim beside the Gait2392 paths they imitate, moved to respect Festo BPA contraction limits and to keep actuators from interfering with one another, and the result was saved as the gait2392\_robotbody model that the later chapters of this dissertation continue to refine \citep{bolen_bipedal_2019}. 
```

**Advisor — 22 October, 11:56 am:**

> I feel like this is jumping all over the place


### BG-29 · line 152 · comment · 22 October, 11:59 am

**Highlighted source:**

```latex
followed
```

**Advisor — 22 October, 11:59 am:**

> is calculated?


### BG-30 · line 152 · comment · 22 October, 11:57 am

**Highlighted source:**

```latex
The biceps femoris long head shows what the comparison exposed:
```

**Advisor — 22 October, 11:57 am:**

> what comparison was made? You don't say


### BG-32 · line 182 · tracked edit · 22 October, 12:05 pm

**Context near the review anchor (not an exact highlight):**

```latex
tor controller developed in Chapter~\ref{ch:sns} is built directly on this two-level organization.


\section{Sensory Feedback in Locomotion}\label{sec:bg_sensory}

Sensory feedback is not a passive correction layer but an integral part of locomotor pattern generation. Reviews of afferent control of the CPG \citep{van_
```

**Pending advisor edit:**

```text
Deleted: \begin{sidewaysfigure}[p] \centering \includegraphics[width=0.97\textheight,height=0.82\textwidth,keepaspectratio]{circuit_literature.pdf} % ALT TEXT: A two-panel synthesis shows a bilateral spinal locomotor network % and a detailed knee-control motif with rhythm-generator, pattern-formation, % motoneuron, commissural, and proprioceptive pathways. \caption[Locomotor circuit architecture synthesized from the literature]{Locomotor circuit architecture synthesized from the literature. The bilateral organization combines rhythm-generation, pattern-formation, motoneuron, commissural, and proprioceptive pathways used across two-level CPG and neuromechanical models \citep{rybak_modelling_2006, rybak_organization_2015, shevtsova_modeling_2016, deng_neuromechanical_2019}. The detailed knee-control motif makes the reciprocal antagonist and sensory pathways explicit.} \label{fig:locomotor_circuit_synthesis} \end{sidewaysfigure}
```


### BG-33 · line 185 · tracked edit · 22 October, 12:04 pm

**Context near the review anchor (not an exact highlight):**

```latex
stems to locomotion are the subject of ongoing study \citep{akay_relative_2021, akay_sensory_2020}. Tonic central and sensory stimuli have been shown to facilitate involuntary air-stepping in humans \citep{selionov_tonic_2009}, demonstrating that simple graded stimulation of motoneuronal and rhythmic circuits can elici
```

**Pending advisor edit:**

```text
Changed: Of particular relevance to the tests planned in Chapter~\ref{ch:futurework}, t
to T
```


### BG-34 · line 185 · tracked edit · 22 October, 12:07 pm

**Context near the review anchor (not an exact highlight):**

```latex
ve and contact feedback channels designed in Chapter~\ref{ch:sns} follow this afferent organization.

\begin{figure}[htbp]
\centering
\includegraphics[width=0.55\textwidth]{rybak_cpg_circuit.pdf}
% ALT TEXT: A three-tier schematic shows the rhythm-generator half-centers
% above a pattern-formation network of interneuro
```

**Pending advisor edit:**

```text
Deleted:  Three published circuits, reproduced in Figures~\ref{fig:deng_connectome}, \ref{fig:afferent_connectomes}, and~\ref{fig:shinohara_connectome}, place these afferent channels in their circuit context. The rat hindlimb controller of \citet{deng_neuromechanical_2019} is organized exactly as Section~\ref{sec:bg_neural} describes: a rhythm-generator half-center and a pattern-formation layer sit above a sensory--motor network in which extensor and flexor muscle receptors feed back through interneurons and Renshaw cells onto the same motoneurons the CPG drives. The human-scale model of \citet{dirusso_investigating_2023} carries the idea to a full lower limb: each agonist--antagonist muscle pair carries Ia, II, and Ib reflex populations with Renshaw cells, and pattern-formation nodes from the locomotor controller terminate on the same motoneuron pools. The cat interlimb model of \citet{shinohara_mechanisms_2025} adds the left--right half of the problem, coupling the two rhythm generators through V3 and commissural interneuron populations while force, length, and velocity feedback returns from every muscle. The schematics differ in animal, scale, and emphasis, but they wire the same afferent classes of this section into the circuit at the level each hypothesis concerns, whether that is the rhythm generator, the pattern-formation layer, or the motor pool. \begin{figure}[htbp] \centering \includegraphics[width=1.0\textwidth]{deng_cpg_connectome.pdf} % ALT TEXT: A four-panel circuit schematic: single-limb control hierarchy, % rhythm-generator configuration, pattern-formation network, and a % sensory-motor network with extensor and flexor receptors, interneurons, % Renshaw cells, motoneurons, and muscles. \caption[Deng two-layer CPG and sensory--motor network]{Synthetic nervous system of the rat hindlimb model: \textbf{(A)} single limb control hierarchy, \textbf{(B)} rhythm-generator configuration, \textbf{(C)} pattern-formation network, and \textbf{(D)} sensory--motor network. Reproduced from \citet[Fig.~4]{deng_neuromechanical_2019} under CC BY 4.0.} \label{fig:deng_connectome} \end{figure} \begin{figure}[htbp] \centering \includegraphics[width=0.62\textwidth]{dirusso_spinal_network.pdf} % ALT TEXT: A spinal network between a muscle and its antagonist with pattern % nodes, interneurons, sensory neurons, motoneurons, and two muscles. \caption[Di Russo spinal network with afferent feedback]{Spinal network between a muscle and its antagonist in the human locomotion model of \citet{dirusso_investigating_2023}: Ia, II, and Ib reflex populations with Renshaw cells, and pattern input from the locomotor controller onto the same motoneuron pools. Reproduced from \citet[Fig.~5]{dirusso_investigating_2023} under CC BY 4.0.} \label{fig:afferent_connectomes} \end{figure} \begin{figure}[htbp] \centering \includegraphics[width=0.62\textwidth]{shinohara_connectome.pdf} % ALT TEXT: A bilateral neuromusculoskeletal cat model: rhythm-generator and % pattern-formation layers, commissural and V3 interneurons, motoneuron % pools, and a cat hindlimb on a treadmill with labeled muscles. \caption[Shinohara bilateral neuromusculoskeletal model]{Neuromusculoskeletal cat model of \citet{shinohara_mechanisms_2025}: bilateral rhythm-generator and pattern-formation layers coupled by V3 and commissural interneurons, motoneuron pools for seven muscles per hindlimb, and force, length, and velocity afferent feedback. Reproduced from \citet[Fig.~1]{shinohara_mechanisms_2025} under CC BY 4.0.} \label{fig:shinohara_connectome} \end{figure} The two-level organization these circuits share descends from the modeling program of the Rybak group, and its foundational schematics are reproduced in Figures~\ref{fig:rybak_cpg_circuit} and~\ref{fig:shevtsova_spinal_circuit}. In the concept model of \citet{rybak_modelling_2006}, afferent feedback ascends to both the rhythm generator and the pattern-formation network, so a sensory perturbation can reset the rhythm at the generator or delete motoneuron activity at the pattern level without necessarily doing both, precisely the deletion phenomenology that motivated the two-level split. The extended model of \citet{shevtsova_modeling_2016} carries the same organization to the two-joint muscles, embedding dedicated pattern-formation populations for the bifunctional posterior biceps--semitendinosus and rectus femoris pools, and that biarticular convergence problem is the one a robot muscle routed across two joints presents to the controllers developed in this dissertation.
```


### BG-35 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex


Computational neuroscience expresses hypotheses about locomotor circuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and exte
```

**Pending advisor edit:**

```text
Changed: locomotor
to computational
```


### BG-36 · line 210 · tracked edit · 22 October, 12:08 pm

**Context near the review anchor (not an exact highlight):**

```latex
otheses about locomotor circuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally
```

**Pending advisor edit:**

```text
Deleted:  differential
```


### BG-37 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
tor circuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central 
```

**Pending advisor edit:**

```text
Deleted: every
```


### BG-38 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
cuitry as executable models. A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populat
```

**Pending advisor edit:**

```text
Added: s
```


### BG-39 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
A computational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \cite
```

**Pending advisor edit:**

```text
Deleted: s
```


### BG-40 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
omputational model is a system of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{b
```

**Pending advisor edit:**

```text
Deleted:  an
```


### BG-41 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
m of coupled equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{brown_intrinsic_1911}. Modern 
```

**Pending advisor edit:**

```text
Added: s
```


### BG-42 · line 210 · tracked edit · 22 October, 12:09 pm

**Context near the review anchor (not an exact highlight):**

```latex
 equations in which elements of the schematic correspond to identified neural populations or pathways. The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{brown_intrinsic_1911}. Modern instantiatio
```

**Pending advisor edit:**

```text
Changed: ; the wiring diagram is not an illustration of the model, it is the model
to s
```


### BG-43 · line 210 · comment · 22 October, 12:10 pm

**Highlighted source:**

```latex
The founding hypothesis is Brown's half-center: alternating flexor and extensor activity arises from two reciprocally inhibitory central populations with opposing bias \citep{brown_intrinsic_1911}.
```

**Advisor — 22 October, 12:10 pm:**

> non seq - this is specifically about CPGs


### BG-44 · line 210 · comment · 22 October, 12:10 pm

**Highlighted source:**

```latex
Modern instantiations of the half-center endow each neuron with persistent sodium channels, whose slow inactivation terminates each burst so that reciprocal inhibition switches the populations without any phasic input \citep{deng_neuromechanical_2019}. Rhythmogenesis, in this tradition, is a property of the intrinsic membrane dynamics, not of a clock elsewhere in the network.
```

**Advisor — 22 October, 12:10 pm:**

> More stuff specific to CPGs, not Comp Neuro


### BG-45 · line 212 · comment · 22 October, 12:11 pm

**Highlighted source:**

```latex
The two-level organization described in Section~\ref{sec:bg_neural} is the field's current working architecture, and its models account for specific electrophysiological observations: the fictive-locomotion deletions reproduced by \citet{rybak_modelling_2006}, the genetically identified interneuron classes mapped onto the model by \citet{rybak_organization_2015}, and the organization of circuits controlling two-joint muscles reviewed by \citet{shevtsova_modeling_2016}. Coupling such circuits to musculoskeletal models shows how afferent feedback shapes CPG output \citep{markin_afferent_2010} and how biomechanics and sensory feedback regularize the behavior of different CPG architectures \citep{deng_biomechanical_2022}. Physiologically plausible models of human locomotion use the same machinery to dissect the relative contributions of reflexes and CPGs, with locomotor primitives, Ia, II, and Ib pathways, and Renshaw cells represented explicitly \citep{dirusso_investigating_2023}, and recent models address interlimb coordination under abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}.
```

**Advisor — 22 October, 12:11 pm:**

> Seems like it should be in the CPG section


### BG-46 · line 214 · tracked edit · 22 October, 12:12 pm

**Context near the review anchor (not an exact highlight):**

```latex
 abrupt perturbations such as the sudden loss of ground support \citep{shinohara_mechanisms_2025}.

One way of interpreting computational neuroscience diagrams is from a controls tradition in which the block diagram is the controller itself. Each block is an identifiable neural population with conductance-based dynamic
```

**Pending advisor edit:**

```text
Changed: For a mechanical engineer, the productive way to read this literature is as
to One way of interpreting computational neuroscience diagrams is from
```


### BG-47 · line 214 · tracked edit · 22 October, 12:12 pm

**Context near the review anchor (not an exact highlight):**

```latex
uroscience diagrams is from a controls tradition in which the block diagram is the controller itself. Each block is an identifiable neural population with conductance-based dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sn
```

**Pending advisor edit:**

```text
Changed: :
to .
```


### BG-48 · line 214 · tracked edit · 22 October, 12:12 pm

**Context near the review anchor (not an exact highlight):**

```latex
oscience diagrams is from a controls tradition in which the block diagram is the controller itself. Each block is an identifiable neural population with conductance-based dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sns}
```

**Pending advisor edit:**

```text
Changed: e
to E
```


### BG-49 · line 214 · tracked edit · 22 October, 12:13 pm

**Context near the review anchor (not an exact highlight):**

```latex
dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sns} comes from this tradition: its two-layer organization follows the rat hindlimb model of \citet{deng_neuromechanical_2019}, and its verification tests are drawn from the d
```

**Pending advisor edit:**

```text
Deleted:  The connection schematics published with these models are therefore worth reading as carefully as their equations.
```


### BG-50 · line 214 · comment · 22 October, 12:13 pm

**Context near the review anchor (not an exact highlight):**

```latex
dynamics, and each connection is a synaptic conductance with a physiologically interpretable weight. The controller developed in Chapter~\ref{ch:sns} comes from this tradition: its two-layer organization follows the rat hindlimb model of \citet{deng_neuromechanical_2019}, and its verification tests are drawn from the d
```

**Advisor — 22 October, 12:13 pm:**

> non seq


### BG-02 · line 260 · tracked edit · 14 October, 2:05 pm

**Context near the review anchor (not an exact highlight):**

```latex
cientific case are the same case: the robot is the instrument that makes these experiments possible.

```

**Pending advisor edit:**

```text
Deleted: these tools were
```


## 50-database.tex

### DB-01 · line 6 · comment · 21 October, 11:13 am

**Highlighted source:**

```latex
and design needs parameters, not only principles
```

**Advisor — 21 October, 11:13 am:**

> huh?


### DB-02 · line 6 · comment · 21 October, 11:13 am

**Highlighted source:**

```latex
For each sensory channel, the designer has to know what the receptor encodes, over what ranges of length, velocity, and force it operates and with what receptor dynamics,
```

**Advisor — 21 October, 11:13 am:**

> not sure what this has to do with the database


### DB-03 · line 6 · comment · 21 October, 11:14 am

**Highlighted source:**

```latex
whether the target is ipsilateral or contralateral
```

**Advisor — 21 October, 11:14 am:**

> there are more options


### DB-04 · line 6 · comment · 21 October, 11:14 am

**Highlighted source:**

```latex
Every element of that specification was established experimentally
```

**Advisor — 21 October, 11:14 am:**

> There are still lots of unknowns


### DB-05 · line 6 · comment · 21 October, 11:16 am

**Highlighted source:**

```latex
literature
```

**Advisor — 21 October, 11:16 am:**

> After this sentence, define the gap in knowledge and what you did to fill it.


### DB-06 · line 6 · tracked edit · 21 October, 11:14 am

**Context near the review anchor (not an exact highlight):**

```latex
se, rat, and human preparations, and the results are distributed over an enormous primary literature. A controller designer is not the only reader stranded by that scatter. An animal researcher planning the next experiment wants to know which preparations established a pathway property, and how the property compares ac
```

**Pending advisor edit:**

```text
Deleted:  rather than collected anywhere
```


### DB-07 · line 6 · comment · 21 October, 11:15 am

**Highlighted source:**

```latex
A controller designer is not the only reader stranded by that scatter. An animal researcher planning the next experiment wants to know which preparations established a pathway property, and how the property compares across species, before committing a preparation to it. A roboticist wants to know which receptor encoding a sensor should imitate. Those readers ask different questions of the same literature, and the literature answers none of them without being gathered, structured, and searchable first.
```

**Advisor — 21 October, 11:15 am:**

> This should probably be the beginning of the paragraph


### DB-08 · line 8 · comment · 21 October, 11:19 am

**Highlighted source:**

```latex
Reviews
```

**Advisor — 21 October, 11:19 am:**

> The tone of this paragraph is as an opening paragraph, and yet it is the 2nd. Probably doesn't need to be here. Discussion of the reviews as a starting point probably goes in your 'methods' section of this chapter


### DB-09 · line 10 · comment · 21 October, 11:20 am

**Highlighted source:**

```latex
Neuroscience 2023 \citep{bolen_sensory_2023}
```

**Advisor — 21 October, 11:20 am:**

> I thought it was presented 2x, once without the AI help, and the other with.


### DB-10 · line 10 · comment · 21 October, 11:20 am

**Highlighted source:**

```latex
distilled mechanism notes
```

**Advisor — 21 October, 11:20 am:**

> huh?


### DB-11 · line 10 · comment · 21 October, 11:21 am

**Highlighted source:**

```latex
first-class record
```

**Advisor — 21 October, 11:21 am:**

> huh?


### DB-12 · line 10 · comment · 21 October, 11:21 am

**Highlighted source:**

```latex
pathway records
```

**Advisor — 21 October, 11:21 am:**

> huh?


### DB-13 · line 10 · comment · 21 October, 11:25 am

**Highlighted source:**

```latex
The database lives in a five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders)
```

**Advisor — 21 October, 11:25 am:**

> I thought you migrated it to a different system


### DB-14 · line 10 · comment · 21 October, 11:25 am

**Highlighted source:**

```latex
canonical
```

**Advisor — 21 October, 11:25 am:**

> database?


### DB-15 · line 12 · comment · 21 October, 11:26 am

**Highlighted source:**

```latex
This chapter and the controller chapter are therefore built as one program.
```

**Advisor — 21 October, 11:26 am:**

> huh?


### DB-16 · line 12 · comment · 21 October, 11:28 am

**Highlighted source:**

```latex
The database informs the sensory channels of Chapter~\ref{ch:sns}: the afferent classes the controller monitors, the signal encoding each channel carries, and the connection rules its circuit implements all trace to records held here, and Section~\ref{sec:sadb_use} documents that mapping.
```

**Advisor — 21 October, 11:28 am:**

> This is a discussion point for the end of the chapter.


### DB-17 · line 13 · tracked edit · 21 October, 11:29 am

**Context near the review anchor (not an exact highlight):**

```latex
 follow describe how the corpus was assembled, audited, pruned, and packaged into its present form.

Methods

The audit and curation effort reported here began from a three-way reconciliation, because the database's three homes did not agree. The author's personal Zotero collection held 538 items, the AARL group librar
```

**Pending advisor edit:**

```text
Added:  Methods
```


### DB-18 · line 14 · comment · 21 October, 11:33 am

**Highlighted source:**

```latex
Methods
```

**Advisor — 21 October, 11:33 am:**

> You don't have to follow this exactly, but you should generally put this in the Chochrane form. https://www.cochrane.org/sites/default/files/uploads/PDFs/guide_to_the_contents_of_a_cochrane_methodology_protocol_and_review.pdf


### DB-19 · line 16 · comment · 21 October, 11:29 am

**Highlighted source:**

```latex
audit
```

**Advisor — 21 October, 11:29 am:**

> ?


### DB-20 · line 16 · comment · 21 October, 11:37 am

**Highlighted source:**

```latex
The author's personal Zotero collection held 538 items, the AARL group library's Biology subtree held 911 items across 27 collections, and the Airtable base held 99 papers. Matching record by record on DOI and normalized title showed that the personal collection was essentially a subset of the group library (all 483 personal DOIs present) and that the Airtable base was a small curated subset rather than a mirror: 95 of its 99 papers resolved to records in the group library's Biology subtree, while 443 personal items had never been entered at all. The gap was imported in verified batches, with every written record re-read and checked after writing, bringing the Papers table to 532 records. The corpus then kept growing while curation caught up: citation mining in the style of Research Rabbit contributed a Crossref-verified digest of recent papers, newly saved references flowed in from the laboratory's Zotero libraries, and a bulk pass imported 372 papers auto-classified from the laboratory's staged collection of PDFs.
```

**Advisor — 21 October, 11:37 am:**

> This is super confusing. We don't need all these details. The airtable base is the database...it can't 'start' with anything if you 'start from the beginning'. Also, starting from the beginning does not mean you tell a chronological story of everything that happened for the database curation. Start with the finished product, then tell us what we need to know about how it could be re-created. 


### DB-21 · line 18 · comment · 21 October, 11:37 am

**Highlighted source:**

```latex
note
```

**Advisor — 21 October, 11:37 am:**

> what are notes?


### DB-22 · line 18 · comment · 21 October, 11:39 am

**Highlighted source:**

```latex
A note is written from the paper's own text, from the abstract at minimum and from the full text when available, never from the curator's recollection of the paper;
```

**Advisor — 21 October, 11:39 am:**

> You are giving details before giving us the big picture. What does that database look like? What is in it? How is it organized?


### DB-23 · line 22 · comment · 21 October, 11:40 am

**Highlighted source:**

```latex
approved
```

**Advisor — 21 October, 11:40 am:**

> approved by who?


### DB-24 · line 22 · comment · 21 October, 11:40 am

**Highlighted source:**

```latex
left as off-topic or orphaned
```

**Advisor — 21 October, 11:40 am:**

> how was this decided?

**Advisor — 21 October, 11:41 am:**

> how were they identified?


### DB-25 · line 24 · comment · 21 October, 11:41 am

**Highlighted source:**

```latex
the laboratory
```

**Advisor — 21 October, 11:41 am:**

> don't you mean you? Take credit in your dissertation.


### DB-26 · line 24 · comment · 21 October, 11:42 am

**Highlighted source:**

```latex
curation tool, not a reader
```

**Advisor — 21 October, 11:42 am:**

> huh?


## 60-sns.tex

### SNS-01 · line 6 · comment · 21 October, 11:48 am

**Highlighted source:**

```latex
Two definitions cover most of the neuroscience vocabulary used below. First, the locomotor controller is a two-level central pattern generator: a rhythm generator sets the timing of the step cycle, and a pattern-formation layer shapes that timing into phase-appropriate activation of antagonist muscle groups \citep{rybak_modelling_2006}. Each leg's rhythm generator is a pair of mutually inhibiting antagonist populations, one extensor and one flexor, so that only one is strongly active at a time. Second, in the non-spiking (graded) formulation used for most of this work, each neuron is electrically a nonlinear RC circuit: synaptic current charges a membrane capacitance through a leak conductance, and the resulting membrane voltage, tens of millivolts, is both the neuron's state and its output \citep{nourse_sns_2023}. Graded chemical synapses couple these units: a depolarized presynaptic cell opens a conductance in the postsynaptic membrane whose sign makes the connection excitatory or inhibitory. The sensory channels that close the loop are the muscle afferents, Ia signaling muscle length and velocity, Ib signaling muscle force, and group II signaling length, together with heel and toe contact sensors.
```

**Advisor — 21 October, 11:48 am:**

> So many details and I don't understand what the goal is. What is the gap in knowledge? What did you do to fill it?

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:09 pm:**

> Gaps: Gait2392_simbody with CPG in Mujoco-SNS Toolbox
> - Gait2392_robotbody (BPAs and BPA routing) with CPG
> - Block-programming on a fast physics simulator (Animatlab no longer supported)
> - Robot solidworks assembly with BPAs from Matlab and a neural controller
> - Li walker with Deng neural controller and foot contact feedback.


### SNS-02 · line 8 · comment · 21 October, 12:45 pm

**Highlighted source:**

```latex
the reflex-based and learned alternatives \citep{geyer_musclereflex_2010, song_deep_2021} are treated there.
```

**Advisor — 21 October, 12:45 pm:**

> huh?


### SNS-03 · line 10 · comment · 21 October, 12:47 pm

**Highlighted source:**

```latex
tooling
```

**Advisor — 21 October, 12:47 pm:**

> tooling tends to be a manufacturing term. Software Tool Development is more like what you are looking for here.


### SNS-04 · line 10 · tracked edit · 21 October, 12:48 pm

**Context near the review anchor (not an exact highlight):**

```latex
hapters~\ref{ch:methods}--\ref{ch:torque}.

This chapter concerns tooling: its object is a verified virtual neuromechanical laboratory for the BPA-actuated robot. In this laboratory a controller can be built in any of several simulation engines, checked against its own re-implementations in the other engines and agains
```

**Pending advisor edit:**

```text
Added: virtual
```


### SNS-05 · line 10 · tracked edit · 21 October, 12:48 pm

**Context near the review anchor (not an exact highlight):**

```latex
tooling: its object is a verified virtual neuromechanical laboratory for the BPA-actuated robot. In this laboratory a controller can be built in any of several simulation engines, checked against its own re-implementations in the other engines and against the published behavior of the models it descends from, and then 
```

**Pending advisor edit:**

```text
Changed: that
to this
```


### SNS-06 · line 10 · comment · 21 October, 12:49 pm

**Highlighted source:**

```latex
built in any of several simulation engines
```

**Advisor — 21 October, 12:49 pm:**

> why would we need several?

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:04 pm:**

> Well, Mujoco doesn't have as good of a physics engine as HYFYDY, for example, and is limited to rigid tendons. Python also is harder to use because it is code without a GUI. Simscape multibody allows for importing parts from solidworks. Some sort of block-programming language like Simulink is a natural choice for building our style of connectome. Mujoco-SNS Toolbox allows for real-time control and computer in the loop commands.


### SNS-07 · line 10 · comment · 21 October, 12:53 pm

**Highlighted source:**

```latex
answer
```

**Advisor — 21 October, 12:53 pm:**

> This is all super awkward. You need to identify the gap in knowledge and how you fill that gap in knowledge. Much of what you wrote here should be saved for discussion, where you talk about the implications of the work.

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:00 pm:**

> I thought that was for research hypothesis not for "here are some tools we made".


### SNS-08 · line 14 · tracked edit · 21 October, 3:34 pm

**Context near the review anchor (not an exact highlight):**

```latex
e closing with the new tools and the program they enable for Chapter~\ref{ch:futurework}.

\section{The AnimatLab Bipedal Walker}\label{sec:prior_walker}

As preliminary work of this dissertation, built before the simulation toolchain of this chapter existed, the author and Morrow constructed a neuromechanical model of
```

**Pending advisor edit:**

```text
Deleted: Preliminary Results:
```


### SNS-09 · line 20 · comment · 21 October, 12:59 pm

**Highlighted source:**

```latex
verification-first approach
```

**Advisor — 21 October, 12:59 pm:**

> I don't understand the connection.


### SNS-10 · line 24 · comment · 21 October, 1:02 pm

**Highlighted source:**

```latex
build order
```

**Advisor — 21 October, 1:02 pm:**

> huh?


### SNS-11 · line 26 · comment · 21 October, 1:00 pm

**Highlighted source:**

```latex
\begin{figure}[htbp]
```

**Advisor — 21 October, 1:00 pm:**

> Overlapping text in figure


### SNS-12 · line 35 · tracked edit · 21 October, 1:04 pm

**Context near the review anchor (not an exact highlight):**

```latex
peline}
\end{figure}

\subsection{Existing Tools and Data Resources}\label{sec:sns_existing_tools}

OpenSim supplies the musculoskeletal modeling side, the robot's source model, and the moment-arm and torque standards the conversions are audited against \citep{delp_opensim_2007, seth_opensim_2018}. MuJoCo supplies the 
```

**Pending advisor edit:**

```text
Deleted: The laboratory stands on tools this program did not build.
```


### SNS-13 · line 37 · comment · 21 October, 1:06 pm

**Highlighted source:**

```latex
The SADb of Chapter~\ref{ch:database} is the data resource among the existing tools \citep{bolen_sensory_2023}. Its segments are the five-table Airtable base (Papers, Feedback, Models, Review Papers, and Recorders), the laboratory's Zotero libraries that feed it, citation mining in the Research Rabbit style that expanded it, and the laboratory's own packaged explorer application that distributes it. The SADb framed this chapter's benchmarks and validation measurements: its pathway records fix what each afferent class encodes and what each connection rule does, so the sensory channels, the reflex topologies, and the acceptance phenomena against which the controllers are validated (the deafferentation and stimulation results cited in Section~\ref{sec:walker_variant_results}) are looked up rather than re-derived from the primary literature.
```

**Advisor — 21 October, 1:06 pm:**

> This whole paragraph seems unnecessary. 


### SNS-14 · line 39 · comment · 21 October, 1:08 pm

**Highlighted source:**

```latex
\subsection{New Tools}\label{sec:sns_new_tools}
```

**Advisor — 21 October, 1:08 pm:**

> Wait, that is all you are saying about the other tools? It is hard for the reader to understand the limitations and gaps without a deeper description of them.


### SNS-15 · line 41 · comment · 21 October, 1:11 pm

**Highlighted source:**

```latex
Simulink blocks
```

**Advisor — 21 October, 1:11 pm:**

> Did you compare and discuss against the libraries myself and Fletcher made here?

**Ben's reply (author context, not an advisor instruction or independent verification) — 21 October, 2:00 pm:**

> No. You've never encouraged us to use Simulink.


## What the advisor's new notes require

The introduction has newly annotated undefined references to numbered gaps, unexplained spiking-neuron mirrors, an overly specific roadmap, and a weak connection to the Biomimetic Humanoid Project. Name the gaps directly and present them in chapter order. Explain Ben's work as a continuation of the humanoid project, rather than inserting adjacent laboratory projects without a clear connection.

The final background batch also identifies CPG-specific material misplaced in the broader computational-neuroscience section (BG-31 and BG-43–45). Move or condense it rather than repeating it in both sections. The new background notes target the motor drawing, the graphical abstract's usefulness, Gait2392's name, the moment-arm discussion, joint-center modeling, run-on sentences, and the chronology of the master's model. State what is compared and how it supports the present work. Replace abstractions such as “first-class,” “definitions are stressed,” and “written against” with concrete descriptions. Do not simply add explanatory paragraphs after every “huh?”; replace or remove the unclear statement and its unnecessary detail.

SADb should start with the finished resource: what it contains, its organization, and who can use it. Then give only the methods needed to reproduce the curation, adapting the advisor's Cochrane suggestion without pretending it was a registered systematic review. Verify the two conference presentations, AI-assisted workflow, decision authority, and actual inclusion criteria from evidence. Ben's role should be explicit. Drop the import chronology and operational record counts unless needed to interpret the final corpus.

For SNS, the old suggestion to compress AnimatLab to one generic paragraph is superseded by Ben's requested causal narrative. The account must identify what Ben combined, what the air-stepping result shows, what failed under hand tuning, and why the new tools were needed. A compact capability/evidence table can replace repeated inventories. Verification should identify the expected behavior, the test or comparison, the observed outcome, and its limit. Do not expand this into a catalogue of all experiments.

## GitHub torque update checked in this pass

Fetch succeeded. The active branch is KneeTestSetup_BenBo_stw. During the pass HEAD and origin both reached **65b05aae161d6a4e87f561a38e57dbae90ec6f95** (22 October, 12:00 Pacific), committed by another session/user; this agent did not commit or pull over local changes. Re-fetch at the start of ZCode work because other machines remain active.

- **fd1c13f4**, 22 October 02:22: updates Knee_Torque_Test/ProofFinal/Main.tex:540–542 with the six-point extensor cam-route sweep, explicitly pending reaction-point geometry and adjoint torque reduction. Read as manuscript status, not independently verified final torque values. Chapter 4 should not continue to say no full angle sweep exists if the test record confirms this sweep.
- The **same manuscript still says “full sweep ... pending” at line 556 and “first ... points ... near full extension” at line 560**. Resolve this internal inconsistency against the actual test sheets before adapting it. Do not copy both claims.
- **7862e4e2**, 22 October 04:31: adds ExtTest_3_01..03 records, updates Results_table_20mm.xlsx, Knee_Extensor_20mm.m, and Ext20mm_compare_20261022.png, and adds the new BiPam_pulley2_20261022.m plus sweep33/escalation records. These are newer than the six-point paragraph; inspect their provenance before selecting the current validation status.
- The 21 October commits also update torque and margin surfaces and the transmission doctrine. PULLEY_SWEEP_DIRECTIVE_20261021.md is the current scope: nBPA and tackle parts each span 1–3, re-solve actuator sizes, account for moving-via load and Bowden compliance, and audit OpenSim path artifacts. The existence of CSVs/mats is not proof that every configuration is feasible or implemented in the simulation model.

The dissertation's 30-torque.tex also differs between live Overleaf and local ProofFinal. Neither was overwritten in this pass. Integrate in a small evidence-backed patch against fresh Overleaf text, not by replacing the chapter with the torque manuscript.

## Figure synchronization in this pass

Uploaded the 11 revised PDF assets to their existing figs/10-introduction, figs/15-background, figs/60-sns, and figs/70-overall folders. Kept editable PPTX, EPS, alt text, and generators local. Chapter 3 assets were not uploaded or edited.

Applied seven exact-fragment changes to figure references, sizing, captions, and the sentence introducing the functional neuron view. No whole chapter was uploaded. Subsequent layout and technical checks led to nine further exact-fragment changes to the two combined figures: shorter captions and 0.94-textwidth panels retain figure text above 10 pt, and “fix axial displacement” precisely identifies the transverse constraint. The separate patches are recorded in the JSON evidence. Advisor tracked edits and comment status were not accepted, rejected, or resolved.

No report, prompt, or administrative file was uploaded to Overleaf; the fresh pre-upload ZIP contains no chatgpt_report file. The source snapshot revealed six divergent chapters (introduction, background, torque, SADb, SNS, and Discussion). Local drafts and live advisor edits must be reconciled, not treated as interchangeable. A fresh source ZIP is evidence at its download time, not a permanent canonical version.

## ZCode execution contract

Use Notes/new_session_prompt_20261022.md. It replaces the stale broad task list with a focused reduction/revision pass, explicit prose and caption limits, the new SNS structure, small patch batches, and evidence requirements. The old prompt and report are preserved locally in Notes/advisor_review_20261022. No chapter reduction has been claimed completed in this pass: ZCode is assigned that heavy revision work.

## Final synchronization checks

- **Live Overleaf:** 295-page current preview; Errors 0, Warnings 4. Both new combined-figure size warnings were eliminated. The four remaining warning cards concern pre-existing floats and a database table outside this figure-sync scope. The advisor's live source was not replaced with the different local draft.
- **Local verification:** ProofFinal/build_advisor_figures_20261022_final/main.pdf, 298 pages; four MiKTeX passes exit 0, zero undefined references/citations, and zero missing files. This local draft differs from the live review copy, so its page count is not a shortening result and must not be presented as the live Overleaf page count.
- Final live SNS and Discussion sources match the fresh pre-upload source plus only the intended figure patches (UTF-8 BOM ignored for comparison). The earlier copied introduction matches its intended single caption patch. No general prose rewrite was applied on Overleaf.
- The final rendered pages for Figures 6.6 and 7.2 pass independent visual and technical review. Minimum editable figure text is 10.152 pt. Their lowest text lies at 672.51 and 671.01 pt on a 792 pt page, respectively. No clipping remains.
- Chapter 3's 19 protected files remain byte-identical. Reports, prompts, native editing sources, and evidence stayed local. No advisor comments were resolved, no tracked edits were dispositioned, and no commits or pushes were made by this agent.
- GitHub was fetched again at completion; HEAD and origin still match 65b05aae. The latest torque-manuscript source commit remains fd1c13f4. No remote change was merged over the dirty working tree.
- Overleaf source/PDF export attempts later in the pass were blocked by Chrome. Final verification therefore uses the visible live PDF, compile UI, and exact source-copy comparisons, plus the independent local build. A completed earlier PDF download was inspectable, but it is not claimed as the final live PDF.

### Per-directive disposition

| Ben's request | Completed action / evidence |
|---|---|
| Review advisor contents, not the whole dissertation | Review Overview capture, exact anchors, pending edits, and relevant reply threads; no new full quality review. |
| Update local CHATGPT_REPORT | Current 73-thread/51-edit inventory above; previous report preserved locally. |
| Insert changed figures | Eleven PDF assets uploaded to the existing chapter figure folders; exact-fragment integration patches and Overleaf proof image saved in Notes/advisor_review_20261022/. |
| Prepare ZCode for reduction and new SNS structure | Replaced Notes/new_session_prompt_20261022.md; explicit AnimatLab → gap → tools/results → verification sequence, limited future-work handoff, and no mandatory Methods/Results/Discussion headings. |
| Avoid more bloated prose and captions | Prompt requires measured cuts, brief captions, comment-specific repairs, and independent review; both combined-figure captions shortened in this pass. |
| Check latest torque work | Fetch and commit/source inspection recorded above, including contradictory sweep status and newer test/sweep records. No unverified numerical transplant. |
| Keep active advisor review intact | No chapter-file uploads, wholesale replacements, comment resolution, or tracked-edit acceptance/rejection; future prose changes remain targeted local proposals. |
| Never save reports on Overleaf | No report/prompt uploaded; only figure PDFs and targeted figure-related source fragments changed there. |
| Protect isometric-force/MATLAB figures | All 19 protected Chapter 3 files retain their hashes. |

Start the new ZCode session using the entire prompt file. It is a handoff for work still to do, not a claim that SADb/SNS reduction or torque revision has already been completed.


## Final Revision C sweep — 2026-10-23

The expanded live Overleaf Review Overview was checked against this ledger. The inventory remains **73 advisor comment threads and 51 pending advisor tracked edits**: Introduction 34 items, Background 49, SADb 26, SNS 15. Ben’s replies and the four Ben-only appendix/main items remain separate. No new open advisor items were found. No Overleaf source, comment resolution, or tracked-edit disposition was changed.

Every advisor ID now has a local disposition and source-line evidence in `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/COMMENT_DISPOSITIONS.md`. The original expanded entries are retained in `review/final_sweep_20261023/recorded_comments.json`. These are editorial dispositions, not claims of advisor approval.

This pass makes the SNS knowledge gap match all five combinations in Ben’s reply, keeps the interconnected bilateral RG construction in the causal opening, adds an early SADb evidence-gap statement, and puts the bubble chart in the SADb main chapter. The current app payload was independently read: it now includes all 1,393 records, so the old paper-only app limitation is superseded. The separate monthly scraper remains under development.

The completed physical-test plots from the latest commits are inserted; no spreadsheet sanity-check torque values are used. Four full fit tables move to Appendix C, with a compact mean-error summary in Chapter 4. The cat caption is shortened, and Figures 2.8 and 2.9 share one page. The active bipulley run is not reported as completed. Original sources and the protected force chapter/figures remain unchanged by this pass.


## Revision C uploaded; Chapter 4 tables restored — 2026-10-23

Uploaded Revision C into the existing Bolen-Dissertation project at Ben's explicit authorization. Replaced individual chapter files and main/support files, preserving folders and review threads. Updated the roadmap, Background figures, torque figures including final physical-test exports, SADb screenshots, SNS figures, Discussion schematics, and Appendix B circuit figures. Existing unchanged force-chapter and reference figures remain. No reports or review artifacts uploaded.

Chapter 4 now contains separate complete pre-optimization, optimized pinned-flexor, and optimized pinned-extensor GoF tables (4.1–4.3), plus a two-row flexor/extensor correction table (4.4), with four parameter columns chi0–chi3 as Ben confirmed. The earlier cross-configuration fit remains in Appendix C. No GoF results were recalculated or changed. The extensor-identified wrapping coefficient is explicitly distinguished from the pinned-flexor fit.

Final Overleaf compilation: Errors 0, Warnings 0; inherited typesetting warnings remain. Overleaf PDF has 259 physical pages; local MiKTeX verification PDF has 257. Both show Tables 4.2 and 4.3 on printed page 96. The Overleaf table page was visually inspected. No source text or comments were added to the protected isometric-force chapter.

Advisor review: 40 straightforward comments resolved; 34 substantive threads remain open, each with a new response describing the revision. This does not claim advisor acceptance. Accepted all 13 remaining grouped tracked-edit markers after replacement; Overview shows no pending tracked edits. The upload regrouped the original 51 advisor tracked edits and Ben's 3 tracked deletions, so grouped-marker counts differ from the original inventory. review_actions.json retains thread text, submitted replies, resolutions, and accepted groups.

The source-ZIP download event could not be verified. A labeled pre-upload history checkpoint was created instead. User authorization covered uploads and routine advisor replies/resolutions. No commit or push was performed.

Evidence: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/review/overleaf_upload_20261023/`.


## Ben's figure/prose correction sweep — 2026-10-23

Working copy: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C`. Originals in ProofFinal remain unchanged. Ben explicitly authorizes Chapter 3 layout changes in this pass. ZCode-owned items are requests to verify, not facts to invent. Progress is updated in this checklist and `review/user_sweep_20261023/checklist.json`.

- [x] U01 [DONE] Remove the quoted stiffness-pair nonuniqueness/external-overshoot-selection paragraph; do not substitute a defensive Pareto justification. Evidence/status: 30-torque.tex: quoted sentence removed; matching Appendix C selection narrative removed. Before snapshot preserves it.
- [x] U02 [DONE] Replace Fig. 2.3 CAD BPA with supplied inflated/deflated photo; preserve originals, give both camera photos descriptive filenames. Evidence/status: Copied original JPEGs without resampling to figs/15-background/festo_bpas_inflated_deflated.jpg and figs/30-torque/cam_extensor_hardware.jpg; original phone files untouched. Photo wired into Fig.2.3; camera EXIF rotation handled in TeX. Visual page check pending.
- [x] U03 [DONE] Fig. 2.4: uppercase panel letters. Evidence/status: steeleknee.pdf: uppercase A/B native PDF overlays; caption uppercase too.
- [x] U04 [DONE] Fig. 2.5: panel letters at upper left. Evidence/status: 15-background.tex: three top-aligned minipages, letters at upper left.
- [x] U05 [DONE] Reorganize Figs. 2.5–2.6: skeleton, robot PAM routes, human Gait2392 paths together; remove repeated models and awkward letters. Evidence/status: Fig.2.5 now skeleton / robot PAM routes / human Gait2392 paths; split original masters_leg_model.pdf into robot_pam_routes.pdf and human_gait2392_paths.pdf. Fig.2.6 replaced duplicate renderings with moment_arm_methods.pdf.
- [x] U06 [DONE] Resize fig:Fmax_fit; fit Figs. 3.3 and 3.4 within page text area without altering MATLAB plot data. Evidence/status: 20-muscles.tex: 3 figures use page floats, centered, bounded panel width/height. Source EPS/MATLAB plots unchanged. Full page verification pending. Full 5-pass build passed; figures3.2–3.4 rendered and visually checked within page bounds.
- [x] U07 [DONE] Center Table 3.3 and increase text size; preserve every value. Evidence/status: Table model_compare unpivoted into 12 readable rows (model per row), centered, small rather than tiny, no negative adjustwidth. Numeric conservation/page check pending. All 36 numeric results checked against original; centered readable table visually verified.
- [x] U08 [DONE] Add diagram connecting Eq. 4.2 dot-product moment-arm construction to cross-product torque method; reuse appropriate Fig. 2.6 geometry where useful. Evidence/status: figs/15-background/moment_arm_methods.pdf: cross-product torque and dot-product moment-arm projection; Fig.2.6 label preserved and cross-referenced near Eq.4.2. Fixed prior self-orthogonal definitions of p_proj versus p_f in nearby prose; equation unchanged. Diagram visually checked.
- [x] U09 [DONE] Fig. 4.4 panel letters bold. Evidence/status: bktFrame.pdf: native text A–D replaced with bold labels; CAD/image content preserved.
- [x] U10 [DONE] Locate and show bending schematic with Eq. 4.24. Evidence/status: No existing BPA-wrapping schematic was located in searched dissertation/paper assets (bracket-bending assets are a different effect). Added native vector bpa_wrapping_loss.pdf as panel B of existing Fig.4.5, referenced at Eq.4.24, preserving later figure numbering. Diagram states fitted loss is not arc-minus-chord length.
- [x] U11 [DONE] Shorten captions of Tables 4.1–4.3; move essential scoring/fold context to prose, preserve values. Evidence/status: Tables4.1–4.3 captions shortened; point-set restrictions and counts moved before Table4.1; fold membership already in rows; excluded flexor test retained in prose. Numbers unchanged.
- [x] U12 [DONE] Table 4.4 chi0 in meters, with correct unit conversion. Evidence/status: Table4.4 header meters; +3.94mm -> +0.00394m, -6.38mm -> -0.00638m. Other parameters unchanged.
- [x] U13 [DONE] Explain anatomical origin/insertion search bounds, biomimetic routing/RoM tradeoff, and subsequent force/part-failure tradeoff. Evidence/status: Added paragraph before route redesign: anatomical bounds, relaxation for RoM, subsequent lower-force hardware tradeoff; cam paragraph simplified. No new bounds or optimizer outcomes invented.
- [x] U14 [DONE 2026-10-23 ZCode] Fig. 4.11 panel C showing torque magnitude across RoM: ZCode-owned; do not duplicate concurrent plot work. Evidence/status: DONE. Panel (C) (Flx20mm_meetsHuman — measured torque vs corrected 620 kPa prediction vs human BFSH maximum, 15 holds −119..+5°) added to fig:FlxGrp in 30-torque.tex; margins +0.42 N·m (@−105°) to +17.0 N·m (@−14°) stated in the panel-B prose paragraph; caption numbers match the torque paper's fig:Flx20mm record and were re-verified by re-running gen_Flx20mm_3panel_20261023.m (15/15 meet-or-exceed). Asset copied byte-identically into figs/30-torque/. A/B panels and caption preserved; label survives once; build fits (float-too-large resolved at 0.72/0.30 widths). Report: Notes/advisor_revision_20261022/Revision_C/review/u14u21_20261023/REPORT.md.
- [x] U15 [DONE] Add supplied cam-knee hardware photograph to Fig. 4.12. Evidence/status: cam_extensor_hardware.jpg copied byte-identically and included below cam route panels in Fig.4.12. Needs full page visual check.
- [x] U16 [DONE] SADb gap: structured database rather than structured account. Evidence/status: 50-database.tex: structured account -> structured database.
- [x] U17 [DONE] SADb source libraries: Zotero libraries of AARL and Ben, plus papers already in Airtable. Evidence/status: 50-database.tex: Zotero libraries of the AARL and myself, and papers already curated in Airtable.
- [x] U18 [DONE] Read SfN24_Bolen.pdf; replace vague AI claims with evidence-based account of limited help and categorization difficulties. Evidence/status: Read SfN24_Bolen.pdf, extracted text and visually inspected poster. Paragraph names ResearchRabbit, Litmaps, txyz, SciSpace, MirrorThink; categorization: none good, txyz adequate, others insufficient. Limited benefits and substantial correction stated. Existing bolen_using_2024 citation retained.
- [x] U19 [DONE] SNS 6.1: initial Connor/Ben controller had ONE RG controlling both sides; air stepping at NIH2021; Ben added afferent and heel/toe feedback; distinguish later bilateral RG expansion. Evidence/status: 60-sns.tex: initial Connor/Ben model ONE RG both sides; air stepping at NIH2021; distinguishes Ben feedback/contact additions and later bilateral interconnected RG model.
- [x] U20 [DONE] State AnimatLab unsupported, based on old physics engine, and slow in this workflow. Evidence/status: 60-sns.tex: AnimatLab unsupported, older physics engine, slow in this workflow.
- [x] U21 [DONE 2026-10-23 ZCode] ZCode verification: Table 6.1 split verification/limits; identify actual authored models versus existing MuJoCo/SNS tools, external MuJoCo–Simulink bridge attribution and direction, controller ports/demos/curriculum scope, resolved welded joint, actual muscle integration status. Do not guess answers. Evidence/status: DONE (verification pass; three text corrections). (1) Bridge attributed (official MathWorks mujoco-simulink-blockset, MuJoCo 2.3.7, AARL keyframe patch + a/b/c tests — BRIDGE_REPORT.md) and direction stated (MuJoCo model embeds inside Simulink; ctrl in, sensors out) in 60-sns.tex §6.4.3. (2) "92 actuators 5 s closed-loop run" corrected to open-loop without afferent return (E2 e2_cpg_open_loop.m); closed-loop credit moved to the E1 reflex experiment. (3) Table 6.1 row-4 limit now states the imported lower-body mechanism's ankle joints remain rigid (digest_2026-09-30: AH_001 ankle DOF open; knee DOF present in the pivoted BA_003 import). Verified-no-edit-needed: authored-vs-existing attributions, W2L rhythm-core numbers (1.000 s, 10/10, r≥0.99996 — CROSS_PLATFORM_REPORT §1.3), Appendix D claims incl. anchor −161.57, ports/demos/curriculum scope, muscle-integration account. Welded-joint claim stays dead (zero hits). Supervisor four-gate GO. Report: Notes/advisor_revision_20261022/Revision_C/review/u14u21_20261023/REPORT.md.
- [x] U22 [DONE] Remove unsupported Hunt/Young neuromechanical Simulink-library comparison; flag prior advisor reply for correction. Evidence/status: Unsupported comparison removed from manuscript. Earlier Overleaf reply edited in place: it no longer claims a Hunt/Young library comparison and leaves comparison with the actual libraries pending verification.
- [x] U23 [DONE] Define W2L at first use, grounded in actual model naming. Evidence/status: Defined Walker Two-Layer (W2L) at first mention; grounded in spinal/w2l_cpg/README.md / Walker_2_Layer_CPG_BilateralRG model naming.
- [x] U24 [DONE] Gap includes tested neuromechanical models with Gait2392-scale actuator counts and degrees of freedom; scope novelty honestly. Evidence/status: 60-sns.tex gap now explicitly includes many degrees of freedom and 92 musculotendon actuators in Gait2392. Does not claim no such models exist anywhere.
- [x] U25 [DONE] Fresh full build, render changed pages, conservation/line-evidence audit, clean upload bundle, sync authorized finished changes to Overleaf. No report uploads. Evidence/status: Fresh 5-pass local build passed. All 36 Table3.3 results and all Table4.1–4.3 data bodies conserved; Chapter3 source plots byte-identical. Changed pages rendered and visually inspected. Uploaded six changed chapters and eight changed/new figure assets individually to Overleaf. Live compilation: Errors0, Warnings0; existing typesetting warnings remain. Live photo and Fig3.3 page inspected; plot and caption fit. Prior misleading reply corrected in place. Reports not uploaded. ZCode-owned U14/U21 remain outstanding.


### Handoff checkpoint: user correction sweep

All edits are in Revision_C, NOT ProofFinal. Nothing from this new sweep has been uploaded to Overleaf yet. Do not run older one-shot revision scripts or packaging scripts: they would overwrite current work or assert the old table/protection state. Current sources are UTF-8; read/write explicitly with UTF-8. Bibliography remains untouched. Source MATLAB/experimental plots remain untouched; Chapter3 layout and table reformatting were explicitly authorized by Ben in this request.

Before snapshot: `Revision_C/review/user_sweep_20261023/before/` (all chapters and figs). Checklist: same folder `checklist.json`; evidence includes poster text/render, asset renders, and original copies. New scripts in root tmp: `revision_c_user_edit.py` (ONE-SHOT, already applied), `revision_c_user_assets.py` (vector PDF authoring, do not rerun on modified assets without restoring snapshot), `revision_c_user_layout2.py` (already applied), `revision_c_user_progress.py` (this checkpoint). Generating vector diagram code is in the assets script. Existing full build helper: `tmp/revision_c_build.py`.

ZCode priorities: Figure4.11 panel C, Table6.1 factual verification and attribution/direction questions, final layout QA, and correcting the prior Overleaf Hunt/Young reply. Do not invent a benchmark or claim Hunt/Young shared a neuromechanical Simulink library. Ben's original technical questions remain in U21 above. Fig4.5 now has wrapping schematic panel B to avoid renumbering the downstream torque figures.


### User correction sweep: final uploaded status

The pre-upload checkpoint above is superseded by this status. All six revised chapter sources and eight figure assets are now on Overleaf. Compile panel shows Errors 0 and Warnings 0; inherited typesetting warnings remain. The live BPA photograph and resized Figure 3.3 were visually checked. The earlier Hunt/Young advisor reply was corrected in place. No new comments were resolved in this sweep; substantive threads remain open.

Completed: U01–U13, U15–U20, U22–U25. ZCode-owned: U14 (Figure4.11 panel C) and factual verification under U21 (Table6.1 attribution, MuJoCo–Simulink bridge direction, actual imported mechanism/actuator integration, and scope of demos/curriculum). U21 layout and the clearly false welded-knee statement were corrected; do not infer that the remaining technical questions have been answered. No performance benchmark was invented.

Latest local review PDF: `Revision_C/Revision_C_user_sweep_20261023.pdf`. Clean upload folder: `Revision_C/OVERLEAF_UPLOAD`. Latest ZIP: `Revision_C/Revision_C_user_sweep_Overleaf.zip`. Evidence, before snapshots, complete diff, checklist, numeric conservation, and screenshots: `Revision_C/review/user_sweep_20261023/`. Original phone photos and MATLAB plot files are preserved. Photo copies have descriptive names. No commit/push or ProofFinal edits.


## U14 + U21 ZCode pass — 2026-10-23

Both ZCode-owned items are done locally against the fresh Overleaf ZIP (downloaded 08:25, extracted to tmp/overleaf_live_20261023). U14: panel (C) added to Fig. 4.11 (measured torque magnitude vs human BFSH across the 15-hold RoM, paper numbers of record re-verified by re-running the generator). U21: bridge attribution/direction stated, the 92-actuator run corrected to open-loop, and the Table 6.1 imports row now carries the ankle-rigid limit; all other table claims verified against campaign records with no correction needed. Supervisor four-gate GO. Local verification build 256 pp, 0 errors, 0 undefined refs; renders of printed pages 101/122/129 checked. NOTHING uploaded; no comments resolved; no commits or pushes. Full record: Notes/advisor_revision_20261022/Revision_C/review/u14u21_20261023/REPORT.md.


<!-- SESSION_1_SIMULATION_RECOVERY_20261027_START -->
## Session 1 simulation recovery update

Reporting window: 26–27 October 2026. Updated 27 October 2026 on **easteregg2**.
Repository root: **`D:\GitHub\Bipedal_Robot`**. Python:
**`D:\Anaconda\envs\myo\python.exe`**. Paths below are relative to that repository
root unless given in full. This block records Session 1's simulation work;
Session 2's entries and the existing advisor/dissertation record are preserved.

### Verified W2L recovery and speed changes

- Corrected gain registration at the evaluator boundary, where the previous
  fresh-process verification discarded per-target settings. Unknown settings
  now raise an error rather than disappearing silently.
- Recovered trial 458 with all effective settings delivered. Its 30 s trajectory
  has 28 left and 21 right contact events, median phase 0.496063, forward
  displacement 2.246231 m, and fitted speed 0.091786 m/s. Four independent
  replays match exactly. All 12 fresh-process perturbations of the saved numeric
  settings by ±1%, excluding the fixed horizon, pass the supported-stepping gate.
- Audited 56 effective conductance classes and three sensor dimensions. Two
  saved settings have no corresponding edges. These are conductance classes,
  not 56 independently tuned anatomical synapses.
- The support ladder reaches nominal `wsup=0`, but the PD harness still carries
  about 16.76% of body weight. No harness-free walking claim is supported.
- Replaced the unsafe sparse-equivalence claim with an exact backend that
  preserves reduction order. Full 30 s traces match the reference exactly;
  isolated neural-step speed improves 1.425×. Independent processes use one
  BLAS thread each. The old FastTrial physics loop remains unsuitable.
- Re-ran perturbation panels with effective gains. Commissural deletion fails
  bilateral stepping; removing all feedback fails the contact-phase gate while
  cycling persists. These are supported W2L adaptation experiments, not an
  original Deng replication.

Evidence: `Code\MuJoCo_SNS\spinal\reports_20261027_codex\REPORT.md`,
`verification.json`, `full_basin_results.json`, `speed_parity_exact_fast.json`,
`perturbation_results.json`, and `w2l_verified_config.json` in that folder.

### Di Russo replication

Found the original author repository, `DiRussoAndrea/Spinal_controller`, and
custom `simgait/SCONE` fork. The previously cited neumovelab repository is a
comparative reimplementation. The original model is the reduced planar
`gait9dof18musc_Thelen_20170320.osim`, with nine active coordinates and 18 Thelen
muscles. It initializes locally in OpenSim 4.4.1. The optimized parameter file
contains 256 parameters, now inventoried.

Corrected multiple equation, sensor, and initialization errors in the existing
Di Russo-inspired adaptation. Its compiled neural implementation matches its
corrected Python equations, but that does not establish equivalence to the
author's controller. A 960-trial duration/speed search selected a somersault;
the candidate is rejected as walking. Further optimization of that surrogate
was stopped. Published walking remains unreplicated, and the original custom
SegmentalController runtime is not available locally.

Original sources are staged at `tmp\codex_dirusso_original\` and
`tmp\codex_simgait_scone\`. Evidence and rejected-result plots are in
`Code\MuJoCo_SNS\spinal\reports_20261027_codex\`.

### Actual robot geometry and corrected conventions

- Refreshed the open `10_AH_001.SLDASM` after Ben adjusted the left knee.
  Preserved both four-bar knees and both spherical hips. Current joint anchors
  meet within numerical precision. Source CAD was not saved by Session 1.
- Corrected Session 1's earlier frame mistake. The current robot pipeline uses
  **x anterior, y superior, and z right**, including gravity, ground plane,
  geometry, inertia, and muscle geometry. A pelvis-part datum transform removes
  the assembly's global repositioning while retaining the knee adjustment.
- Ben superseded the lumbar hinge with a **ball-and-socket joint**. The spine
  remains a rigid 0.25 m bar. Payload center is 100 mm anterior and 75 mm above
  the joint. The current model splits the balancing payload equally between
  torso and a rigid pelvis/ilium load.
- Loaded `Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design,
  Newest (7-17-20)\Foot Ankle.SLDASM`: all 22 components resolved with no open
  errors or warnings. Exported the geometry, mirrored it for the opposite side,
  and welded all foot internal joints and the ankle attachment to each shank.
- Total mass remains 75 kg: 5.356664 kg CAD lower body with feet, 0.12 kg spine,
  and 69.523336 kg balancing payload. Each payload portion is 34.761668 kg.
  Neutral torso extension demand is 34.1012 N·m. Printed parts use solid Onyx;
  hardware densities and the 20 mm square spine cross-section remain engineering
  assumptions rather than measured hardware properties.

The current model is `Code\MuJoCo_SNS\robot_torso\ah001_neural.xml`.
Earlier z-up routing files, including `opensim_routes.json` and
`ah001_with_opensim_routes.xml`, are superseded and must not be used as the
current robot model. Legacy converted human/W2L models were not globally
rotated; their verified baselines remain separate.

### Robotbody routes and neural control

The source is now **`Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim`**,
using its right-side routes and explicit left mirrors. The existing
`Code\Matlab\Mesh_Optimization\Routing27_v2\actuator_map_27.json` selects
**27 actuators per side, 54 total**. No new hip-muscle consolidation was performed.
The full source inventory retains other route definitions for auditing; they
are not additional robot actuators.

Moving path points use exact piecewise cubic functions exported from OpenSim,
including derivative and endpoint checks. Conditional membership updates at
runtime. Complete path Jacobians include moving-guide terms in force and
velocity calculations. Checks cover 492 source function evaluations, 25 native
conditional checks, and 75 constrained poses. Function errors are below
1.4 × 10⁻¹⁷ m, and full path-derivative errors are below 1.4 × 10⁻⁹ m.

The updated audit reports **zero reversed neutral sagittal actions**. Mirrored
attachment corrections for adductor magnus 2, gracilis, and semimembranosus are
recorded explicitly; maximum point shifts are 15.6, 21.0, and 18.1 mm. The
previous TFL reversal clears with the refreshed CAD pose. Conditional switching
of unselected semitendinosus routes has small length jumps after retargeting;
those routes are not actuators in the 27-per-side model and remain flagged.

Constructed a **262-neuron, 230-input SNS controller** driving the 54 actuator
controls, with Ia/II/Ib and contact inputs, postural force allocation, and
RG-driven stepping references. Reference tensions have a 30 ms activation
filter and the existing 27-actuator map's MIF caps, including 2500 N for each
erector. These are simulated reference capacities, not a calibrated BPA force
model. The latest test mode uses explicitly compensated MN target currents;
it is a hybrid controller, not a faithful Deng circuit replication.

**Current control status: a 20-second supported crouched-stance test is verified;
locomotion remains unverified.** Initial tests collapsed or lost posture. After
correcting force allocation, motor-neuron tracking, and the explicit support rig,
the final five seconds of the successful stance run had knee excursions of
1.35° right and 0.95° left, maximum pelvis tilt of 3.92°, and maximum lumbar
rotation of 0.58°. Mean foot loads were 17.63 N right and 21.61 N left. Actual
vertical rig support averaged **94.68% of body weight** in that interval, with
angular stabilization also active. This establishes supported postural control,
not independent balance or harness-free standing. Disabling all muscle forces
under the same support settings caused lumbar failure after 0.202 seconds.

The reference and exact accelerated neural backends produced **bit-identical
complete 20-second stance trajectories**, including positions, velocities,
activation, force, contact, route, and neural logs. Evidence is in
`Code\MuJoCo_SNS\robot_torso\controller_assessment.json` and
`runs\stand_20s_reference.*` / `runs\stand_20s_exact.*`. A replay image and
animation are `robot_torso\neural_stance.png` and `neural_stance.gif`.

The original stepping drive latched rather than oscillated. A separate rhythm
probe found bilateral oscillation at drive 3.0, and the corresponding supported
robot run completed 20 seconds without a numerical or posture failure. Its
whole-run forward displacement was 0.446 m, and mean actual vertical support
was 93.74% of body weight. These totals include settling and do **not** establish
alternating gait or locomotion; contact, phase, and displacement after stepping
onset still require assessment. Reducing vertical rig stiffness from 800 to
300 N/m failed the pelvis-height gate after 4.352 seconds, before stepping
onset. These runs are `runs\step_drive3.*` and
`runs\step_drive3_lower_rig.*`. Task 3's `s3k`, `w2lvar`, and `syn6` campaign
has not been launched.

Current evidence: `Code\MuJoCo_SNS\robot_torso\robotbody_routes.json`,
`robotbody_source_inventory.json`, `route_sign_corrections.json`,
`dynamic_route_verification.json`, `run_neural_robot.py`, and `runs\`.
Fresh CAD/foot inventories are in
`Code\MuJoCo_SNS\spinal\reports_20261027_robot_control\`.

### Scope and next steps

Continue stepping diagnostics and inspect contact, phase, and support forces
before classifying locomotion. Supported stance and complete-trajectory backend
parity are now verified within the assistance conditions stated above. Updated
ZCode prompts name easteregg2 and provide full paths, environment, and runnable
commands in
`Code\MuJoCo_SNS\spinal\reports_20261027_robot_control\ZCODE_PROMPTS_SESSION1.md`.
No
commits, pushes, outgoing ZCode messages, or dissertation source edits were
made in this simulation-recovery work. This Session 1 block is additive; it does
not replace Session 2's work or resolve advisor comment threads.
<!-- SESSION_1_SIMULATION_RECOVERY_20261027_END -->

