# Figures to make — notes for Ben (2026-09-07; status pass 2026-09-23)

**STATUS PASS 2026-09-23 (ZCode, easteregg2):** two items are now DONE by existing
assets (2 and 6A), one partially absorbed (4), one downgraded (5). Still genuinely
missing and high-value: **item 1 (roadmap) and item 3 (pipeline)**. Note the Xi-figure
files below live ONLY on EB475WS4 (untracked) — they must be committed/pushed and
uploaded to Overleaf before the text that references them compiles there.

The monograph currently has figures for Chapters 3–5 (reused from the two papers)
but **none for Chapters 1, 2, 6, or the appendices**. Below are specs for the
figures the text needs or would strongly benefit from, in priority order. Once
you make them, upload to Overleaf (suggest `figs/New/`) and tell me — I'll write
the `\begin{figure}` blocks with captions/labels and insert them at the marked
locations. Specs include suggested format; EPS or PDF preferred for consistency
with existing figures, PNG acceptable.

---

## 1. Dissertation roadmap (Chapter 1) — HIGH priority — **STILL OPEN**

- **Where:** end of §1.5 Dissertation Organization (p. 6–7). Label `fig:roadmap`.
- **Content:** left-to-right block diagram of the program:
  `Human biomechanics benchmark (OpenSim Gait2392)` → `Actuator characterization (Ch. 3–4)` →
  `Joint torque validation (Ch. 4)` → `MuJoCo plant (Ch. 6)` → `SNS controller: 2-layer CPG (Ch. 6)` →
  `Treadmill robot (Ch. 6)`.
  Under the first three blocks: "validated in this dissertation"; under the last
  three: "future work". Optionally a feedback arrow from the robot block back to
  the benchmark labeled "biological insight".
  2026-09-23 note: the "future work" labels predate the completed-work framing rule
  (summary.md ground rule 4) — relabel as what the tools let researchers do now vs
  what they define the path toward, per the I5 organization paragraph.
- **Form:** clean vector boxes, 1 arrow style, no color needed for print (black
  outlines, minimal fills). Fits single column (\textwidth, ~4–5 cm tall).
- **Why:** committees read §1.5 first; this is the one-figure summary of the whole
  dissertation.

## 2. Two-level CPG architecture (Chapter 2 / Chapter 6) — HIGH priority — **DONE for Ch. 6; residual for Ch. 2**

- **Status 2026-09-23:** built and wired into `CPG_spinal_section_draft.tex`:
  `CPG_airstepping_figs/circuit_literature.pdf` (reader-level overview, `fig:cpg-overview`,
  correct shapes per Ben's convention) and `circuit_dengstyle.pdf` (structure-driven,
  every displayed edge asserted against the compiled network). Both exist in the repo.
- **Residual:** Background §2.6 (entry B1 in summary.md) may still want the classic
  "adapted from \citet{rybak_modelling_2006}" borrow for the background chapter —
  cheap, cited, and distinct from our own circuit figures. Ben's call.

## 3. Simulation-to-robot pipeline (Chapter 6) — HIGH priority — **STILL OPEN**

- **Where:** §6.1 A Neuromechanical Simulation Pipeline. Label `fig:pipeline`.
- **Content:** pipeline diagram:
  `Modified OpenSim model (.osim, Robot2392 routing)` → `MyoConverter` → `MuJoCo MJCF
  (muscle actuators w/ activation state)` ⇄ `SNS-Toolbox network (Python)`.
  Inputs/outputs between MuJoCo ⇄ SNS: "Ia, Ib, II, heel/toe contact → SNS" and
  "per-muscle activation → MuJoCo". A second arrow continues: `→ embedded stack:
  Teensy (per-joint loops) + Jetson (SNS) → treadmill robot`.
  Mark the phases: "model conversion", "closed-loop simulation", "hardware deployment".
- **Form:** vector, full \textwidth; matches roadmap figure style from item 1.
- **Why:** §6.1–6.4 is a verbal pipeline today; one figure makes the whole Future
  Work chapter concrete.

## 4. Sensory feedback detail (Chapter 6) — MEDIUM priority — **PARTIALLY ABSORBED**

- **Status 2026-09-23:** the afferent classes and their projection targets are drawn
  in the circuit figures of item 2 (Ia/Ib/II/heel-toe edges in `circuit_literature.pdf`).
  A dedicated single-leg anatomical panel (femur/tibia + spindle/GTO/contact pads) is
  now OPTIONAL — build only if a reviewer wants the anatomy view.

## 5. Optimization-method schematic for the improved torque model (Chapter 3) — LOW priority (downgraded 2026-09-23)

- **Status 2026-09-23:** the force-balance solve and the wrap-loss mechanism now have
  purpose-built figures (`xiBalance.pdf`, `xiWrapLoss.pdf`; see summary.md 2026-09-21
  sections — files on EB475WS4, untracked). The remaining unmet piece is only the
  optimizer-loop flowchart (variables → evaluators → gamultiobj → GoF); build it only
  if the Methods section still reads confusing after M6's restructure lands.

## 6. Appendix A/B illustrations — LOW priority

- **A: DONE 2026-09-21** — `xiFrameGeo.pdf` panel B is the generic bracket-frame
  construction figure this item asked for (plus in-situ frames in panel A). File on
  EB475WS4 (untracked), wired into `ZCode_drafts/chapters/20-methods.tex`.
- **B: OPTIONAL** — a one-panel neuron/synapse schematic matching Table B.1; raw
  material exists (`sns_diagram_panels.png` in `CPG_airstepping_figs/`, SNS_Library
  figures in `Code\Matlab\SNS_Simscape\figures\`) but nothing dissertation-formatted.

---

### Already-have list (updated 2026-09-23 — no action needed)
- Test stand photo/figure: `01_testJigs1.pdf` (Ch. 3)
- Knee test jigs, ICR plot, bracket frames, all torque results figures (Ch. 3–4)
- Model comparison figure (Ch. 5)
- Steele knee paragraph figure: `steeleknee.pdf`/`.png` (Dissertation root, B9; panel C still awaits Ben's test-stand photo)
- Xi methods trio: `xiFrameGeo`/`xiBalance`/`xiWrapLoss` .pdf+.eps (EB475WS4 only, untracked)
- CPG circuit set + air/ground walk + FSA synergy figures: `CPG_airstepping_figs/`
  (`circuit_literature`, `circuit_dengstyle`, `hindlimb_style_nap_air`, `curr3i_*`,
  `fsa_*`)

### Caption/style reminders (PSU)
- Every new figure needs a List-of-Figures-quality caption (number, title, page).
- Schematics adapted from published figures must say "adapted from" + citation in
  the caption.
- Keep 12 pt equivalent text size inside figures where possible.
- Project-wide figure standards (Ben, 2026-09-21): Tol palette from
  `Code\Matlab\Colors.m`, Arial only, 10 pt floor, no italics, CVD-safe
  markers/line styles, alt-text file per figure set — full rules in `AGENTS.md`.
