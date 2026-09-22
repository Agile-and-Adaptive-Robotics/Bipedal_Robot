# CPG dissertation update — copy/paste map for Ben

> **STATE ADDENDUM 2026-09-21 (laptop ZCode, after the EB475WS4 pulls):**
> The draft below was updated on the other machine (commits 4086acb, a83b865,
> a4738b9): `\fillme` count is down 12 → **8**; a new reader-level overview
> figure `fig:cpg-overview` (`circuit_literature.pdf`) now precedes the
> edge-audit `circuit_dengstyle.pdf`; air-stepping numbers were corrected
> (2.108 s period, knee −106.1°, cycle-mean extensor duty 0.691); the
> Curriculum Results subsection carries real stage-1 numbers (137.254 at
> trial 78, 29 cycles, drive 3.154 nA, τh 170 ms) and states honestly that
> the curriculum has not converged at all three stages; new ground-walk
> figures `curr3c_gait_isb.png` / `curr3c_gait_contact.png` are referenced
> and exist. **Every figure the draft references now exists in
> `CPG_airstepping_figs/`** — the old "overlay STALE-ABSENT" item is dead
> (replaced by the curr3c/ISB figures). The PART A–D paste map below is
> unchanged and still valid; the remaining blockers are the 8 fillmes and
> PART D's key-collision check. The draft's paste targets say
> `ProofFinal/`, but ProofFinal is stale (Sep 8) vs Overleaf — paste into
> Overleaf, not ProofFinal.

Everything lives in **`CPG_spinal_section_draft.tex`** (Dissertation
root, next to this file). It is organized as four labeled blocks.
The `\fillme{...}` slots are result numbers I fill when the curriculum
finishes — grep the tex for `fillme` before you paste anything final.

## Where to paste what

| Block in the draft | Paste into (ProofFinal/) | Where exactly |
|---|---|---|
| **PART A — METHODS** (one `\section` + 5 subsections) | `chapters/20-methods.tex` | Append at END of file (after the *Simulink SNS Library and Reflex Topology* subsection) |
| **PART B — RESULTS** (one `\section` + 4 subsections + 2 figures) | `chapters/30-results.tex` | Append at END of file (after *Follow-Up Identification and Route Redesign*) |
| **PART C — DISCUSSION** (one `\section` + 3 subsections) | `chapters/40-discussion.tex` | Append at END of file (after *Synthesis*) |
| **PART D — BIB** (17 entries) | `thesis.bib` | Append at END; **check for key collisions first** — `rybak2006`, `mccrea2008`, `nourse2023`, `derrusoijspeert2023`, `rybak2015` may already exist from earlier chapters. If they do, keep YOUR existing keys and fix my `\citep` calls to match |

## Figures referenced (files already in `CPG_airstepping_figs/`)

| Draft label | File | Status |
|---|---|---|
| `fig:cpg-circuit` | `circuit_dengstyle.pdf` | DONE (NaP architecture, edge-verified) |
| `fig:nap-air` | `hindlimb_style_nap_air.pdf` | DONE (untuned air walk) |
| `fig:cpg-overlay` | `opensim_overlay_gait_cycles.png` | STALE-ABSENT — regenerates from the stage-3 winner run (post-curriculum deliverables step); I will refresh it + the caption numbers |
| (ground walk) | GIF/frame export TBD | I generate after the curriculum; LaTeX can't inline GIFs — I'll export a frames strip or you keep the GIF as supplementary media |

Also available if you want them: `sns_diagram_panels.png` (toolbox
subnetwork panels — nice for an appendix), `air_nap.gif` +
`muscle_force_compare.png`.

## Bibliography cautions

- Entries marked **VERIFY** in their `note` field: I wrote them from
  the PDFs' extracted text + your Zotero keys, but could not confirm
  volume/pages from inside the PDFs — 30 seconds each in Zotero:
  `shinohara2025`, `shevtsova2026`, `perreault2011` (NC7ETD46),
  `jankowska2010` (8RDG5YUZ), `zhang2022` (KIG9DEKJ),
  `dominguez2020` (YANK2JE7), `rybak2024`, `rybak2025`,
  `talpalar2013`.
- The old draft's v-class correspondence table (`tab:vclass`) is NOT in
  the rewrite — the V3/c1 rows changed and the table should be rebuilt
  from `spinal\DESIGN.md` (2026-09-16 section) if you still want it.
  Say the word and I'll regenerate it as a proper LaTeX table.

## What I still owe you (after the curriculum lands)

1. Fill every `\fillme{...}` with the stage-1/2/3 scores + final
   ground-walk numbers (I'll do this in the same file, then tell you
   it's ready to paste).
2. Regenerate `opensim_overlay_gait_cycles.png` + ground-walk figure
   from the stage-3 winner, refresh captions.
3. Regenerate `circuit_dengstyle.pdf` with the winner's gains (the
   figure auto-loads from `curriculum_stage3.json`).
4. Commit the results (second commit) and remind you about pushing.

— ZCode session, 2026-09-16
