# Figures to make — notes for Ben (2026-09-07)

The monograph currently has figures for Chapters 3–5 (reused from the two papers)
but **none for Chapters 1, 2, 6, or the appendices**. Below are specs for the
figures the text needs or would strongly benefit from, in priority order. Once
you make them, upload to Overleaf (suggest `figs/New/`) and tell me — I'll write
the `\begin{figure}` blocks with captions/labels and insert them at the marked
locations. Specs include suggested format; EPS or PDF preferred for consistency
with existing figures, PNG acceptable.

---

## 1. Dissertation roadmap (Chapter 1) — HIGH priority

- **Where:** end of §1.5 Dissertation Organization (p. 6–7). Label `fig:roadmap`.
- **Content:** left-to-right block diagram of the program:
  `Human biomechanics benchmark (OpenSim Gait2392)` → `Actuator characterization (Ch. 3–4)` →
  `Joint torque validation (Ch. 4)` → `MuJoCo plant (Ch. 6)` → `SNS controller: 2-layer CPG (Ch. 6)` →
  `Treadmill robot (Ch. 6)`.
  Under the first three blocks: "validated in this dissertation"; under the last
  three: "future work". Optionally a feedback arrow from the robot block back to
  the benchmark labeled "biological insight".
- **Form:** clean vector boxes, 1 arrow style, no color needed for print (black
  outlines, minimal fills). Fits single column (\textwidth, ~4–5 cm tall).
- **Why:** committees read §1.5 first; this is the one-figure summary of the whole
  dissertation.

## 2. Two-level CPG architecture (Chapter 2 / Chapter 6) — HIGH priority

- **Where:** §2.6 Neural Control of Locomotion (referenced conceptually now), reused
  again at §6.2.1 Architecture. Label `fig:cpg` (or `fig:cpg_arch`).
- **Content:** schematic of the Rybak two-layer CPG for ONE leg plus the
  left–right link:
  - Rhythm Generator (RG) half-center: RG-F / RG-E populations with reciprocal inhibition.
  - Pattern Formation (PF) layer: PF-F / PF-E populations with reciprocal inhibition.
  - Motoneuron pools (MN-F, MN-E) → two muscles (flexor/extensor BPA symbols).
  - Commissural pathway to a ghosted mirror RG for the contralateral leg (left–right coordination).
  - Sensory inputs drawn in color or dashed: Ia, Ib, II onto PF/MN, heel/toe contact onto phase switching.
- **Form:** redraw after Rybak et al. 2006 / McCrea & Rybak 2008 (their classic
  figure) but with OUR notation (RG/PF/Σ1–Σ2 names as in those papers), captioned
  "adapted from \citet{rybak_modelling_2006}". Full \textwidth.
- **Why:** the CPG architecture is described twice in words; the committee needs
  the picture. This is the single most important missing figure.
- **Data source:** none — schematic.

## 3. Simulation-to-robot pipeline (Chapter 6) — HIGH priority

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

## 4. Sensory feedback detail (Chapter 6) — MEDIUM priority

- **Where:** §6.2.2 Sensory Feedback Channels. Label `fig:sensory`.
- **Content:** single leg (sagittal view of femur/tibia is fine as simple line
  segments) with the four afferent classes drawn from muscle spindles (Ia, II),
  GTO (Ib), and heel/toe contact pads, with their projection targets in the
  PF/MN circuit of figure 2. Could be combined with item 2 as a two-panel figure
  (A: architecture; B: feedback detail) — your call.
- **Data source:** schematic; muscle geometry can be traced from the test-stand
  CAD screenshots (e.g., `Figures/bktFrame` sources).

## 5. Optimization-method schematic for the improved torque model (Chapter 3) — MEDIUM

- **Where:** §3.7 Improved Torque Model, after eq. (complianceforcebalance). Label `fig:optloop`.
- **Content:** flowchart of one optimizer evaluation: route geometry →
  l_LMT, r_k → force-balance solve (δ) → torque prediction → GoF (RMSE, FVU, max
  residual) → gamultiobj/surrogateopt decision variables (χ0, χ1, χ2, χ3).
- **Why:** the χ terms' meaning is easier to show than tell, and it documents the
  surrogateopt upgrade you're planning (you can add a branch "GA (used here) /
  surrogateopt (in progress)").
- **Data source:** `minimizeFlxPin.m` / `buildKneeFlexorContext20mm.m` structure
  (I extracted the full data flow on 2026-09-07; ask if you want the call graph).

## 6. Appendix A/B illustrations — LOW priority (nice-to-have)

- **A:** small diagram of the bracket frame (x̂_br, ŷ_br, ẑ_br) and force
  direction û with the projection onto the line of action — supports the
  compliance appendix. Could reuse/redraw `bktFrame.pdf` annotation.
- **B:** one-panel neuron/synapse schematic (leaky integrator with Ia/Ib/II and
  graded synapse symbols, matching Table B.1 values).

---

### Already-have list (no action needed)
- Test stand photo/figure: `01_testJigs1.pdf` (Ch. 3)
- Knee test jigs, ICR plot, bracket frames, all torque results figures (Ch. 3–4)
- Model comparison figure (Ch. 5)

### Caption/style reminders (PSU)
- Every new figure needs a List-of-Figures-quality caption (number, title, page).
- Schematics adapted from published figures must say "adapted from" + citation in
  the caption.
- Keep 12 pt equivalent text size inside figures where possible.
