# Dissertation staging package - 2026-09-25

Prepared for Ben's review and subsequently merged into the local `ProofFinal/`
working copy. Overleaf and the review-comment state were not changed.

## Staged changes

- `10-introduction-objectives.tex`: replaces only Research Objectives 2 and 3.
- `15-background-steele.tex`: replaces the B9 paragraph and placeholder figure block.
- `20-methods-xi-figures.tex`: inserts the corrected projection and frame-transform
  figures without restoring the discarded Xi figures.
- `xiProjection.pdf` and `xiProjection.svg`: a tall portrait construction with two 3-D bracket frames,
  opposed force vectors, local bracket displacements written with $\eta$, clearly
  offset force-path projections, and the
  tendon elongation at bracket 2. The PDF metadata and the LaTeX block carry matching
  alt text; colors are color-vision-deficiency safe and redundant with labels and line
  weights.
- `xiFrameTransform.pdf`: two stacked geometric panels showing the fixed-space Z
  rotation toward the sagittal projection and the current-body Y rotation toward the
  full attachment point. The transformation matrix is now in the surrounding text,
  rather than embedded in the figure.
- `Figures/Aim2/steeleknee.pdf`: Steele et al.'s original two-panel mechanism and mechanical-stop
  image (their Fig. 6), retaining its native `(a)` and `(b)` designations with no
  added headings or schematic. Attribution is in the staged LaTeX caption under
  CC BY 4.0. The matching citation entry is staged in
  `steele-2023-reference.bib`.
- `PLACEHOLDER_AUDIT.md`: submission blockers and stale proposal text.

PNG files are visual-review copies. Regenerate the review assets with:

```text
C:\Users\Ben Bolen\.conda\envs\myo\python.exe build_staged_figures.py
```

## Merge status

The approved objective wording, Xi methods text, Steele background paragraph,
figure captions, and bibliography entries have been merged into the local
`ProofFinal/` tree. The rejected legacy Xi figures and the unfinished
`FLAG-BEN` human-torque subsection were not merged.
