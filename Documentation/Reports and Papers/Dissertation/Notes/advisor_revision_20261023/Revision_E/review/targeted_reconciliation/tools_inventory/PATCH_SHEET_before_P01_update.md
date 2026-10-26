# Targeted Overleaf edits — current coordination sheet

**Overleaf is the working document. Do not upload or paste entire local E chapters.** Local E/PDFs are reference drafts that contain proposed prose, newer figures, and supplemental evidence. They do not include the live review history or every author change. Whole-chapter replacement is not the application method.

## Coordination

- Ben owns the live document and author decisions. Preserve his tracked deletions and wording changes. A deleted passage stays deleted unless he explicitly reopens it.
- ChatGPT supplies bounded prose patches, linked to the advisor’s exact comment and existing text, with a reason and evidence limit.
- ZCode supplies source-verified modeling results and implementation details. Its ongoing experiments should arrive as one evidence packet, rather than parallel chapter rewrites. This file is a local handoff; no message was sent to ZCode.
- Apply each approved patch in place with tracking retained. Do not resolve comments or accept/reject existing changes as part of a prose patch. After application, check the associated thread still exists and its new location makes sense. Replacing a commented span can move its anchor; in-place editing is not a guarantee that an anchor remains unchanged.
- Match exact existing text, section label, and comment author/date. If the match is absent or changed, stop that patch and compare the new author text. Line numbers are supplementary navigation only.
- A large section can be a bounded patch when necessary, with explicit start/end text and a list of affected threads. A whole chapter is not the default patch unit.

## Three issues raised now

**Deletion:** the live tracked deletion is the sentence beginning “This suspended stepping result preceded…”. The remainder of the hand-tuning/outcome paragraph is still present. Ben replied “deleted.” The local E draft incorrectly put back a paraphrase, “This suspended result concerns…”. That paraphrase has now been removed locally. P00 records the instruction; nothing needs to be reinserted into Overleaf.

**Connections:** “framework” names the overall toolset. It is not a replacement for a particular signal exchange or mechanical coupling. Use “synaptic connections” for circuit edges, “interfaces” for exchanges between software/model components, and “coupling” for their operation together. In the criticized paragraph, replace the vague two concluding sentences with P01. Preserve the rest of the live paragraph for the next bounded review; this patch does not silently endorse or rewrite its five-system inventory.

**Capability table:** the current Overleaf source does contain a table labeled `tab:sns_capability`, and the Software Tools paragraph references it. The earlier paragraph calls it “the capability table” without a reference. P01 fixes that sentence while retaining the live table. Local E deliberately removed the table. Do not mix the table-free E paragraph with the table-present Overleaf passage. Removing or redesigning the table, if chosen later, must be a separate coordinated patch affecting its mentions and label together.

## Exact patches

### P00 — KEEP DELETED; no new text

**Comment:** ajh26, 23 October 10:51 am: Huh? Just delete all this?; Ben 26 October 5:20 am: deleted

**Find exactly:**

```tex
 This suspended stepping result preceded the later feedback and bilateral-RG extensions; it does not demonstrate walking with those additions.
```

Deletion verified live. Do not restore this sentence or a paraphrase. Retain the rest of the outcome paragraph unless Ben directs otherwise.

### P01 — Replace only these final two sentences in Knowledge Gaps and the New Tools

**Comment:** ajh26, 23 October 11:21 am: need a better word; second marker *

**Find exactly:**

```tex
I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.
```

**Replace with:**

```tex
I developed the model conversions, actuator-force calculations, Simulink blocks, and controller implementations described below. Table~\ref{tab:sns_capability} identifies the interfaces tested and the model integrations that remain incomplete.
```

Exact text read in current Overleaf. Table and label both present. Keep table while applying this patch; do not combine with the table-free local E wording.

### P02 — Replace only item 4 wording; retain Ben’s block-programming description

**Comment:** Ben 26 October 2:21 am tracked change and 5:15 am reply to huh?

**Find exactly:**

```tex
I built a GUI for block-programming neural controllers that were attached mechanical models
```

**Replace with:**

```tex
I built a GUI for block-programming neural controllers coupled to mechanical models
```

Exact current text read live; grammar/coupling wording only.

### P03 — Rename only the subsection title; preserve label

**Comment:** Ben’s explicit heading instruction in chat

**Find exactly:**

```tex
\subsection{Graphical Circuit Construction}\label{sec:sns_new_tools}
```

**Replace with:**

```tex
\subsection{Block Programming Neural Pathways}\label{sec:sns_new_tools}
```

Exact current heading read live. Preserve the existing label and comments.

## Previously requested edits still need bounded application

Figure 6.1’s (A)/(B) labels, Ben’s original shared-RG image, the later per-side RG chronology, and the conversion-repair appendix are prepared locally. They are not automatically present in Overleaf. Their source and evidence are in `review/figure61_correction/`. Apply the figure block, the chronology paragraph, and the appendix insertion as separate patches, preserving the existing document and review threads. The main-text appendix reference and appendix section must be applied together to avoid an unresolved reference. No chapter swap is required.

The larger clarity proposals in `review/clarity_pass/` remain proposals. A row marked “editorial treatment supplied” meant the local draft was edited, not that the live advisor comment was resolved. Those rows must not be used as a certificate of live compliance. This sheet takes precedence for application strategy.

No Overleaf source or review action was changed during this reconciliation. `live_verification.json` records what was actually observed. The replacement `.tex` files here contain fragments only, not complete chapters.

## Bounded fragments for the already requested figure/history/conversion edits

- **P04_figure61_block.tex:** replace only the figure environment containing `\label{fig:animatlab_model_architecture}`, from its own `\begin{figure}` through `\end{figure}`. Upload the supplied `connor_ben_shared_rg_original.png` asset to the existing figure directory before compiling. Preserve the label. This is not a source-file replacement.
- **P05_controller_history.tex:** replace only the paragraph beginning “I subsequently expanded the neural model into left- and right-side circuits” and ending with the Appendix simulation-verification reference. The replacement preserves Ben’s later RG/PF and contact-feedback additions and corrects the chronology. If that paragraph has since changed, compare it before applying.
- **P06_appendix_conversion_section.tex:** insert this section into the supplementary-verification appendix before “Model Construction and Recording,” only if that appendix exists and `sec:app_conversion_repairs` is absent. Do not replace the appendix file or create a second copy of its chapter.
- **P07_conversion_summary.tex:** replace only the paragraph beginning “MyoConverter supplies the initial OpenSim-to-MuJoCo translation” and ending with the 2 ms time-step sentence. This expands that paragraph into the concrete conversion summary. It stops before the existing BPA-interface paragraph. Apply P06 with it, or omit the new appendix-reference sentence until the target section exists.

These fragments are prepared for targeted review/application; none has been applied to Overleaf. P04/P05 follow Ben’s explicit figure and chronology correction. P06/P07 carry the source-supported repair details. Preserve and recheck each affected advisor thread after editing its span. This still does not constitute approval or completion of the remaining larger clarity proposals.

## Verification after reconciliation

The local reference draft was rebuilt after removing the reintroduced qualifier: 269 pages, zero errors, and zero unresolved references/citations; 26 document-wide overflow warnings remain. The compiled file is Revision_E_reference_reconciled.pdf and is a reference only, not a chapter-replacement instruction. P01–P07 remain proposed fragments for the live document. No Overleaf edits were made. P00 is already present as Ben's tracked deletion and needs no new insertion or deletion online.
