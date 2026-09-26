# Dissertation placeholder audit - staged 2026-09-25

This audit treats live Overleaf as authoritative. It does not assume that either
`ProofFinal/` or `ZCode_drafts/` is current. No authoritative chapter was edited.

## Submission-blocking placeholders

1. `ZCode_drafts/chapters/02-abstract.tex`: two `FLAG-BEN` values for the human
   knee torque result. Do not paste this paragraph until the measured range-of-motion
   percentage and RMSE are supplied.
2. `ZCode_drafts/chapters/30-results.tex`: the entire Human Knee Torque Validation
   subsection still contains `FLAG-BEN` values and a missing final comparison figure.
   It is not paste-ready.
3. `CPG_spinal_section_draft.tex` contains no active `\fillme{...}` calls. Its
   remaining `\fillme` strings are comments documenting slots that were filled. The
   draft still states that the curriculum did not converge at all three stages, which
   is a scientific limitation to retain or revise deliberately, not a placeholder.

## Placeholders cleared by this staging package

1. `ZCode_drafts/chapters/15-background.tex` contains an empty panel-C box in the
   proposed Steele-knee figure. The staged `steeleknee.pdf` removes that empty panel,
   and `15-background-steele.tex` supplies the matching two-panel caption.
2. The older Xi figures were dropped after review. The staged `xiProjection.pdf` and
   `xiFrameTransform.pdf` implement Ben's replacement specification without restoring
   `xiBalance` or `xiWrapLoss`.

## Stale or misleading future-work framing

1. Research Objective 3 still says "Define a path to neuromechanical control" and
   points only to Chapter `ch:futurework`. Replace it with the staged completed-tools
   objective and the Methods/Results/Discussion cross-references.
2. `ZCode_drafts/chapters/10-introduction.tex` still assigns the working
   neuromechanical pipeline to the Future Work chapter. The organization paragraph
   needs a separate surgical rewrite after the new CPG sections are placed.
3. `ProofFinal/chapters/15-background.tex`, `20-methods.tex`, `30-results.tex`,
   `40-discussion.tex`, `50-futurework.tex`, `60-conclusion.tex`, and
   `93-AppendixB.tex` contain older "planned" or "remaining neural layer" language.
   Do not bulk-replace it. Update each sentence only after confirming the live
   Overleaf wording and the final placement of Sections `sec:spinal-cpg`,
   `sec:results-cpg`, and `sec:disc-cpg`.

## Proposed text that should not be pasted as written

1. `ZCode_drafts/chapters/20-methods.tex` references `xiFrameGeo`, `xiBalance`, and
   `xiWrapLoss`. Those references and floats are superseded by the staged two-figure
   package.
2. The same Methods draft says the extensor's extra correction follows because a
   bracket forms an arch and extends 75 mm in Z. This does not clearly state that the
   BPA wraps and loses usable length; retain the existing equation but rewrite that
   rationale before pasting.
3. The Abstract claims that the complete SNS was verified against published tonic-
   stimulation and deletion phenomena and that lesion/stimulation experiments can be
   run "today." Recheck those claims against the final CPG Results section before use.
4. The Background proposal says the controller "can be verified" before tuning. That
   is procedural/future framing and should be replaced with the tests actually run.
5. The flexor evaluator source is not uniform: `minimizeFlxPin.m` computes the
   two-rotation matrices but its active transform lines use the one-rotation frame,
   whereas the two-rotation evaluator uses $R_zR_y^T$. The staged transform text is
   therefore explicitly scoped to the two-rotation evaluator. Do not describe it as
   the universal implementation without reconciling those active code paths.

## Safe, immediately reviewable items

1. `10-introduction-objectives.tex`: exact replacement for objectives 2 and 3.
2. `15-background-steele.tex` plus staged `steeleknee.pdf`: removes the visible
   placeholder and points the knee contribution to Methods, Results, and Discussion.
3. `20-methods-xi-figures.tex` plus the two staged Xi PDFs: replacement figures and
   insertion text only; no optimizer values or equations are changed.
