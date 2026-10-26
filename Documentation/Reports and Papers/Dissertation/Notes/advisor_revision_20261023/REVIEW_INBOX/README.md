# ZCode → ChatGPT dissertation review inbox

Use this repository-relative location on either machine:

`Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/REVIEW_INBOX/<machine>/<batch>/`

Check hostname first. Use a short descriptive batch name with a time stamp, for example `20261026_0610_sns-demos`. Write each batch to a new directory; do not overwrite a batch already under review. Machine directories prevent the two sessions from overwriting each other.

## What to deliver

```
<machine>/<batch>/
  REVIEW.md
  patches/
    Pxx_name.tex
  figures/
    figure_name.pdf
    figure_name.png
    figure_name.alt.txt
```

**REVIEW.md** should state:

1. The machine, repository commit, and relevant uncommitted source versions; identify the exact files/results used.
2. Each proposed patch’s target chapter, section label, exact existing start/end text, and associated advisor comment (quote plus author/date). Current Overleaf text is the application baseline. Local line numbers are only navigation and must name their source snapshot.
3. What changed and why. Mark it as a prose-only edit, an author-confirmed fact, a newly checked result, or a proposal. Preserve Ben’s tracked deletions and replies. Do not reintroduce their content in a paraphrase.
4. For numerical results: model, controller, effective parameter file, support/protocol, metric definitions, values, and source output paths. A concise existing-evidence summary is enough; do not launch extra experiments merely to populate this handoff.
5. For each figure: intended figure number/label, whether it replaces an existing figure or is a proposed addition, caption text, panel identities and (A)/(B) letters, units/sign conventions, and generating-script/result paths. Include a PNG preview and a vector PDF when available. For an unchanged author-supplied image, preserve the original bytes.
6. Validation performed: compilation if relevant, missing references/assets, and any known interpretation or formatting issue. Distinguish a saved build from a merely proposed edit.

**patches/** contains bounded replacement or insertion fragments, not whole-chapter uploads. For a substantial section rewrite, include exact boundaries and affected comment threads in REVIEW.md. Keep labels stable. The associated table and its references must be revised together. Optional full-chapter previews may be put in `preview/` and explicitly marked REFERENCE ONLY; they are not instructions to replace the live chapter.

**figures/** contains reviewable assets with stable names; record any proposed filename change and every LaTeX reference it affects. Do not scatter new assets into the reviewed Revision_E directories before review.

## Current editorial rulings

- P01 begins with “Researchers,” and Gait2392-based neuromechanical models are second in the contribution list.
- Describe delivered tools with selected functional-verification/proof-of-concept examples, not an outstanding-work checklist or exhaustive per-tool requirement.
- Credit MuJoCo, SNS-Toolbox, Optuna, Simulink, and the MathWorks bridge; claim the custom models, blocks, integrations, and evaluation procedures.
- Ben confirms the coupled SolidWorks-derived robot with neural control and BPAs using his force model. Its saved model is on the unavailable laptop; do not make access to it a condition for listing that contribution or invent measurements for it.
- Rybak’s identified implementation is Python. Scope each published-model claim to its actual implementation and demonstrated behavior.
- Figure 6.1(B) is Ben’s supplied original Connor–Ben shared-RG diagram. The per-side RG split and RG/PF/contact feedback are Ben’s later extension. Do not substitute the later diagram for the original again.
- Keep the “This suspended stepping result preceded…” tracked deletion deleted. The live capability table still exists; do not mix table-free E prose with table-present Overleaf text.
- No Overleaf edits, comment resolutions, accepting/rejecting tracked edits, commits, or pushes as part of delivery.

## Handoff

Tell Ben the machine and batch path when ready. A file existing only on easteregg2 is not automatically present on EB475WS4. ChatGPT can retrieve the named batch through the existing read-only SSH route, or Ben can synchronize it through his normal process. Do not commit or push just to deliver a batch unless Ben authorizes it. The inbox instructions do not authorize messaging another chat or agent.

ChatGPT reviews each batch, supplies bounded prose/application patches, and records what is accepted, revised, or excluded. Ben controls application to the live tracked Overleaf document. Keep ongoing modeling and incoming proposals separate from the reviewed manuscript.
