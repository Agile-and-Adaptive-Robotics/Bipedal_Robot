

## Visual verification (ink-bounds analysis, 110 dpi renders; visual-judge unavailable, fallback per instructions)

Control pages (untouched chapters) define the normal text-block ink envelope:
ink x[110,826] px of a 935 px page, right margin 109 px (approximately 72 pt margins).

| Page (PDF) | Content | Result |
|---|---|---|
| 26-27 | Intro gap + objectives rewrite | Fits envelope |
| 142 | Torque nine-point sweep paragraph | Fits envelope |
| 143, 145, 146 | SADb opener, methods, Explorer figure | Fits envelope |
| 153-154 | tab:sadb_afferents + tab:sadb_pathways | Sideways tables (rotated, centered) — correct in baseline and proposal; NOT defects |
| 156, 158 | SNS opener, gaps section | Fits envelope |
| 164-165 | tab:sns_capability | DEFECT FOUND on first render (ink to page edge, 934/935 px; ~100 cell overfull boxes up to 397 pt) → FIXED: `llX` natural-width columns replaced by fixed ragged-right `p{0.28/0.24}+X` at `\small`. Re-render: ink x[110,826], identical to controls |
| 179 | Verification section | Fits envelope |
| 196 | tab:spike_animatlab | PRE-EXISTING 24 pt overflow (baseline PDF p.203 identical, spans to x=564.8 pt) → FIXED with single `\small` (natural width 493→450 pt; no text or number changes). Now x=542.1 pt, 1.4 pt residual |
| 201 | Controller comparison | Fits envelope |

Overfull-hbox width multisets (final pass): baseline 148, proposal 146
(2 baseline widths removed by rewrites; 0 added).

Rebuild after fixes: 284 pages, 0 errors, 0 undefined references/citations,
0 missing figure files. Proposal build log = `tmp/build_live4.log`; renders
`v3_p164.png` / `v3_p196.png`.
