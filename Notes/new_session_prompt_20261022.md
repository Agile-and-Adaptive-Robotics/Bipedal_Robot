# Focused dissertation reduction and advisor revision

Paste the prompt below into a new ZCode session. The prior prompt is preserved in `Notes/advisor_review_20261022/new_session_prompt_20261022_before_this_pass.md`.

---

Work in D:\Github\Bipedal_Robot. Shorten and improve Ben Bolen's dissertation, chiefly SADb and SNS, address the advisor's latest comments, and reconcile Chapter 4 with verified recent torque-paper work. The advisor is actively reviewing Overleaf. Do the heavy local revision work and deliver targeted patches. Do not upload prose, reports, prompts, or administrative files to Overleaf. Do not commit or push.

## Establish the baseline

Check hostname first, then read the dissertation skill, AGENTS.md, required campaign context, and **CHATGPT_REPORT.md**. That report contains 73 open advisor threads, 76 advisor messages, and 51 tracked edits captured from live Overleaf on 22 October, with anchors and Ben's separate replies. The advisor added comments during capture; counts are a snapshot. Current instructions below supersede conflicting old campaign plans.

Fetch GitHub, inspect status, and preserve all existing changes. Other machines are active. At handoff HEAD and origin/KneeTestSetup_BenBo_stw matched 65b05aae. Use GitHub Desktop's bundled git on EB475WS4. Do not pull/reset/restore/stash another session's work to tidy the tree.

Get a fresh Overleaf ZIP before drafting patches. Introduction, background, torque, SADb, SNS, and Discussion differ from local ProofFinal. Compare them and preserve local-only revisions separately. Stage proposals against a dated copy of fresh live source. Never replace a live chapter with a stale local copy. Source ZIPs do not contain the review-system comments or pending-edit history. Evidence is in Notes/advisor_review_20261022/ and tmp/advisor_review_20261022/.

## SNS: Ben's new controlling structure

Do **not** impose Methods, Results, and Discussion headings. Tell this argument in order:

1. **The AnimatLab combination.** Ben combined the Deng CPG with the Li walker during his PhD. Explain the combination, show stepping-like motion in air, and state that hand tuning did not produce walking. This belongs here, not Background. Retain enough detail to motivate the next step, rather than reducing it to a generic paragraph.
2. **What that attempt exposed.** Use the gaps Ben identified in his comment replies: human and BPA-routed robot models coupled to neural control; the Li walker and Deng controller with contact feedback; accessible graphical circuit construction on maintained simulation tools; and CAD-derived mechanics, BPA models, and neural controllers connected together. Explain the actual limitations of existing tools. Replies establish Ben's intent, not proof that every capability was completed.
3. **New software tools and brief results.** Explain what Ben built, what was inherited, and what the combination adds. A compact capability/evidence table can replace repeated inventories. Distinguish completed components, partial imports, provisional parameters, and proposed integration.
4. **Verification.** Select a few convincing tests. State expected behavior, test/comparison, observed result, and limit. Show that the tools work and address the gap without cataloguing every experiment. Include the advisor's requested comparison with Hunt/Fletcher libraries using actual sources, not an invented comparison.
5. **Short handoff.** Summarize delivered capabilities and remaining limitations. Detailed future work belongs in Chapter 7.

Retain “no harness-free walking demonstrated” unless new verified evidence changes it. Removing tonic drive is not removing physical support. Draft circuits, rule translations, imports, GUI existence, and compilation are not verified locomotion. Do not start a new simulation or optimization campaign to fill the chapter.

## Shorten the document

First create a small cut map: section/paragraph, keep/replace/remove/move, purpose, comments addressed, and estimated words saved. Then execute bounded local batches; do not ask permission for every routine cut already authorized by Ben.

Planning targets: roughly 40–50% less SNS narrative and 25–35% less SADb narrative against fresh live text. These are working goals, not quotas. Preserve the argument and report a justified different result rather than padding or cutting blindly. Measure narrative words, caption words, chapter pages, and total PDF pages before/after. Moving text into an appendix does not reduce dissertation length. Keep only essential reproducibility material there; operational history belongs in existing project records.

SADb starts with the **finished product**: contents, organization, table/record relationships, and useful queries. Follow with only the methods needed to recreate it, meaningful results, and limitations. Adapt the advisor's Cochrane suggestion without claiming a registered systematic review. Drop import chronology, batch bookkeeping, saved-file stories, repeated segment lists, and repeated corpus counts. Verify the two conference presentations, AI-assisted workflow, decision authority, and actual inclusion criteria. Define “notes,” give Ben credit for his work, and distinguish retained from audited corpora. Never invent missing decisions or reconcile disagreeing numbers by guessing.

## Address the “huh?” passages

Address the report's current “huh?”, undefined-term, awkward, non-sequitur, and run-on comments, including introduction and background. Do not conduct an unrelated stylistic rewrite of all chapters.

- Replace the unclear claim, or remove it if unnecessary. Do not append a long explanation to a confusing paragraph.
- Name gaps and models directly, rather than “the first two,” “the fourth,” or “the same construction.” Introduce necessary terms before using them. Match the order and language of gaps, objectives, and chapter roadmap.
- Explain how prior work connects to the Biomimetic Humanoid Project. Avoid unrelated project lists.
- Remove administrative narration and abstractions: saved files, canonical records, “first-class,” “definitions are stressed,” repeated software inventories, and repeated claims of usefulness.
- The latest background notes also identify CPG-specific material inside the computational-neuroscience section. Put that material in the CPG discussion or cut duplication; keep the broader section's scope clear.
- Write direct prose for a mechanical-engineering committee. Use first person for Ben's own work where appropriate, as the advisor requested; this supersedes the older “never I” rule. Leave unrelated paragraphs alone.
- Every paragraph must do a necessary job. Cut the second explanation of the same contribution.

## Torque: reconcile recent work carefully

Re-fetch before this step. Inspect Knee_Torque_Test/ProofFinal/Main.tex and its supporting test/optimization records. Starting points:

- fd1c13f4 (22 October 02:22), Main.tex:540–542: six-point extensor cam-route sweep, with reaction-point geometry and adjoint reduction pending.
- The same manuscript at line 556 still says the full sweep is pending, and line 560 refers only to near-extension points. Resolve those contradictions against test records before adapting anything.
- 7862e4e2 (22 October 04:31): later ExtTest_3_01..03 data, updated Results_table_20mm.xlsx, Knee_Extensor_20mm.m, Ext20mm_compare_20261022.png, BiPam_pulley2_20261022.m, and sweep33/escalation outputs.
- PULLEY_SWEEP_DIRECTIVE_20261021.md covers configuration selection, actuator resizing, moving-via loads, Bowden compliance/packaging, and OpenSim path artifacts. CSV/mat existence does not establish feasibility or dynamic-model integration.

Use the current reaction-point constants and load-cell angle convention. Separate raw measurements, provisional reductions, finalized torque comparisons, and optimizer predictions. Check pulley-surface captions against their generators and documented axes. Adapt only dissertation-relevant changes with exact provenance for every changed number. Do not substitute the manuscript for Chapter 4.

## Deliver small, reviewable patches

Prefer one section or one to three paragraphs per patch. Each batch needs directive/comment IDs, file and section, exact BEFORE and AFTER text, a unified diff against the fresh live baseline, a context anchor and baseline hash, a concise reason, and before/after word counts. Log every removed/moved paragraph, numerical claim, figure block, and label. Shortening is authorized; transparent conservation and relocation accounting are required.

Use honest status labels: proposed locally, applied locally, verified, or awaiting evidence. Do not resolve Overleaf comments or accept/reject tracked edits. If restructuring SNS needs a larger section patch, give a move map and preserve paragraphs that still work. Whole-chapter replacement is a last option: explain why smaller patches fail and leave it as a proposal for Ben.

## Figures and captions

Revised schematics have already been uploaded. Do not redo them. PPTX/EPS sources and alt text remain in chapter components folders. See Dissertation/Notes/figure_polish_20261022/REPORT.md and the later synchronization notes in CHATGPT_REPORT.md. The two combined figures now have shorter captions and 0.94-textwidth panels to retain text above 10 pt while fitting pages.

**Do not change any Chapter 3 isometric-force figures or MATLAB figures.** Avoid broad figure factories. If recommending a removal or move, name the figure and its argumentative cost; preserve its source and repair cross-references. The advisor's graphical-abstract question does not authorize silent deletion of protected Chapter 3 assets.

Captions should usually be 25–60 words: panels, symbols, essential conditions, and one necessary limitation. Use more only when a complex figure requires it. Do not repeat the results or methods paragraph. Use a short optional caption for front-matter lists. Check final print size. No TikZ schematics; preserve native PowerPoint sources.

## Verify and finish

Preserve bibliography bytes (Latin-1 byte-safe editing only), Oxford commas, engineering diameter notation, Modern Robotics conventions, and “model,” not “plant.” No machine/campaign narration in dissertation prose. Do not silently change numbers. Do not kill the user's applications or other sessions to fix a build lock.

For each coherent batch, run the full MiKTeX sequence in a fresh D:\temp build with absolute BIBINPUTS/BSTINPUTS; save the verified PDF in a new build folder. Render changed pages and check captions, figure fit, labels, and references. Inspect warnings in touched sections. Tectonic cannot handle existing Chapter 3 EPS files; do not convert protected assets merely to satisfy the app preview.

Use independent technical-writing review of changed sections and visual review of altered figure pages. Judge actual text and pages, not checklist compliance alone. Read the result as a committee member: is the purpose, contribution, evidence, and limit clear on the first pass? Fix failures without adding unnecessary text.

Finish with targeted patch files, the verified PDF, before/after length measurements, a per-comment/directive compliance table with line evidence, unresolved factual questions, and exact files for Ben to commit. Keep all reports/prompts local. Do not upload proposed prose while the advisor is reviewing without Ben's later explicit authorization.
