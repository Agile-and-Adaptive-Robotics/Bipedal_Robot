# U14 + U21 — ZCode execution report (2026-10-23)

Working copy: fresh Overleaf ZIP (Bolen_Dissertation.zip, downloaded 2026-10-23 08:25)
extracted to `D:\Github\Bipedal_Robot\tmp\overleaf_live_20261023\`. All 12 differing
chapters are content-identical to Revision_C (line-endings only), so the live Overleaf
source is the batch baseline. Local verification build:
`D:\temp\diss_u14u21_20261023\main.pdf` (256 pp; pdflatex + bibtex from the output dir
+ pdflatex ×2 — 0 errors, 0 undefined references/citations, 0 missing files, 0
float-too-large warnings). Overleaf-live page count differs by engine (established
2–3 pp drift), so 256 is not a shortening result.

## U14 (verbatim)

> Fig. 4.11 panel C showing torque magnitude across RoM: ZCode-owned; do not duplicate
> concurrent plot work.

Target = `fig:FlxGrp` (Fig. 4.11), `chapters/30-torque.tex` (11th figure environment of
ch:torque; cross-checked against U09=bktFrame 4.4 and U15=cam photo 4.12).

Edits (final state):

1. `30-torque.tex:475-477` — appended to the panel-B paragraph (wording per Ben's
   2026-10-23 correction: "data points", not "holds"):
   "Panel \textbf{C} shows the outcome for the redesigned two-BPA flexor route of
   Section~\ref{sec:redesign_physical}: across its fifteen data points from \ang{-119}
   to \ang{5}, the measured torque magnitude met or exceeded the human
   biceps-femoris short-head maximum at every point, with margins from
   \qty{0.42}{\newton\meter} at \ang{-105} to \qty{17.0}{\newton\meter} at \ang{-14}."
2. `30-torque.tex:479-489` — figure block: FlxGrp.eps include moved to
   `width=0.72\textwidth`, then `\par\smallskip`, then a centered 0.30\textwidth
   minipage carrying `\textbf{(C)}` above `Flx20mm_meetsHuman.pdf` (copied
   byte-identically from `Documentation\Reports and Papers\Knee_Torque_Test\Figures\`
   into `figs/30-torque/`). Panels A/B caption text preserved byte-for-byte. Caption
   (C) sentence, finalized per Ben's disposition (he applied the stub
   "\textbf{(C)} The redesigned two-BPA flexor route" on Overleaf; this is the short
   completion): "\textbf{(C)} The redesigned two-BPA flexor route
   ($\diameter$\qty{20}{\mm}, $l_{0}=\qty{56.2}{\cm}$;
   Section~\ref{sec:redesign_physical}) with the corrected \qty{620}{\kPa} prediction
   (magenta dashed) and the maximum human isometric torque (black dotted); indigo
   circles are measured torque, sized by gauge pressure." Label `fig:FlxGrp` survives
   exactly once.
3. Layout history (evidence of iteration): 1.0\textwidth stack → "Float too large for
   page by 195.02pt"; 0.46/short-caption → 135.54pt; 0.72/0.30 → 0.396pt; 0.72/0.30 +
   smallskip → fits, warning gone. Class metrics probed: \textheight 548.5pt,
   \textwidth 469.75pt. Rendered page (printed p.101 = physical 120) visually checked:
   A/B intact and legible, (C) placed and labeled, full caption on the page.

Number sources of record: generator `tmp\paper_bridge\figscripts\
gen_Flx20mm_3panel_20261023.m` re-run headless this session (MATLAB R2025a; stdout:
15/15 meet-or-exceed; smallest margin +0.42 N·m at −105.0°; largest +17.00 at −14.0°;
rest = 0.5615 m → 56.2 cm; corrected [min max] = [−91.91, −6.64]; rigid peak −190.34)
and the torque-paper caption of record `Knee_Torque_Test\ProofFinal\Main.tex:544,561`
(fig:Flx20mm). Asset regenerated today by that run (mtime 08:56; earlier 05:58 copies
superseded in place — same script, same data). χ values in the run match Table 4.4
(+3.93737 mm, 3.998e4, 1.473e4, 0.158289).

Layout caveat for Ben: the 548.5pt text height cannot hold A/B at 1.0\textwidth plus a
floor-legal (C) panel plus the preserved caption — A/B print ≈28% smaller than the
live figure. If you want A/B full-size instead, the alternatives are (a) regenerating
FlxGrp as a true three-panel asset from `Figure components\FlxGrp\FlxGrp.fig`
(plot work is ZCode-owned per U14), or (b) moving (C) to a companion figure
(renumbers 4.12+). Current state is the minimal-churn solution; say the word and I
will re-cut it.

## U21 (verbatim, abridged)

> ZCode verification: Table 6.1 split verification/limits; … external MuJoCo–Simulink
> bridge attribution and direction … resolved welded joint, actual muscle integration
> status. Do not guess answers. … Do NOT revive welded-joint claim.

Verification results (all source-traced, none guessed):

| Question | Finding | Evidence |
|---|---|---|
| Bridge attribution | Official MathWorks `mujoco-simulink-blockset` (mathworks-robotics), compiled against MuJoCo 2.3.7, with an AARL local patch (keyframe-0 start) + sensor-patched model copy + a/b/c test suite | `Code\Matlab\SNS_Simscape\mujoco_bridge\BRIDGE_REPORT.md:3,109-152` |
| Bridge direction | MuJoCo model embeds INSIDE Simulink: S-function steps `mj_step`, ctrl in (width = n actuators), sensor bus out | BRIDGE_REPORT.md:32-68 |
| "92 actuators, 5 s closed-loop" claim | WRONG as written — that run is E2 `e2_cpg_open_loop.m` (all 92, 5 s, open-loop, no afferent return); the closed-loop experiment is E1 (one motor pool, sensors→afferents→MN) | BRIDGE_REPORT.md:219-253 |
| Table 6.1 row 4 (mechanical imports) | Import achieved is true; the honest limit is the humanoid import's ankle DOF (rigidly connected follower ports; 59 world-tie cuts needed); knee DOF present in the pivoted BA_003 leg-rig import (3 revolute) | `Documentation\Project_Digests\digest_2026-09-30.md:127-141`; `dev\` logs add_ankles11/rig_rebuild5/surgery3 |
| W2L rhythm-core numbers (§6.4.2) | Verified exact: 1.000 s periods, 10/10 bursts, r ≥ 0.99996 at −2 ms, latch without heel contacts | `campaigns\20260930\CROSS_PLATFORM_REPORT.md` §1.3 + `sns_w2l.log` quote |
| Appendix D claims | Deng port 1.938 s vs 1.94 NumPy, r −0.86 ✓; full-net verify 4.2e-06 mV ✓; 1,446 blocks/149 joints ✓; AnimatLab 11/11 ✓; BPA march 70% harness/880/1372 N/372 kPa ✓; wean −161.57→−145.07 ✓ (matches anchor −161.56754…) | Appendix D vs CROSS_PLATFORM_REPORT §1, README_SNS_Simscape, AGENTS.md records |
| Controller ports / demos / curriculum | Existing text accurate, no edit needed: generated network = 376 external-input ports + 92-wide output (export_network_json.py / sns_build_from_json.m); Appendix D `sec:sns_ports` + `sec:sns_tuning_methods` cover port checks, TPE staging, 1% basin checks; spiking twin comparisons back row 2's "controller comparisons" | `95-AppendixD.tex`; README_SNS_Simscape |
| Muscle integration | Existing text accurate: 92 muscle-tendon actuators in converted MuJoCo; BPA quasi-static interface (bpa_muscle.py Festo port); Simscape BPA blocks remain provisional (BPACPGLegDemo 0 sustained cycles at stock params — honestly stated in campaign records, not overstated in ch. 6) | §6.4.1; CROSS_PLATFORM_REPORT §1.4 |
| Authored vs existing tools | Table row names already distinguish them; §6.3-6.4 attributions verified accurate (MyoConverter initial translation + Ben's repair/audit; library retains SNS-Toolbox formulation; bridge = MathWorks blockset) | `60-sns.tex:29,53,59` |
| Welded joint | Stays deleted — zero "weld" hits in any live chapter | findstr over `chapters\` |

Edits made (final state):

1. `60-sns.tex:46` — Table 6.1 row 4 Limit cell:
   BEFORE "Standalone reflex demonstrations do not establish control of imported models"
   AFTER "The imported lower-body mechanism's ankle joints remain rigid, and standalone
   reflex demonstrations do not establish control of imported models"
2. `60-sns.tex:109` — §6.4.3 paragraph:
   - attribution/direction clause added: "The MuJoCo--Simulink bridge, built on the
     official MathWorks blockset that embeds the MuJoCo model inside Simulink as the
     mechanical model, was tested against…"
   - factual fix: "…drove all 92 actuators during a \qty{5}{\second} open-loop run,
     without afferent return; the closed-loop reflex experiment exercised commands and
     sensory measurements through the same interface."
   - ("official", not "open-source": the blockset license carries a field-of-use
     clause — supervisor Gate 3 escalation, adopted.)

## Gates

Four-gate supervisor run: **GO** (no gate failed) on the 09:10 state; dispositions:
directed own-manuscript number reuse (panel C data is the manuscript's own, publication
note intact); zero external verbatim runs; slop sweep clean; scope = exactly the four
targeted hunks, no drive-by rewrites. After that audit, three further mechanical
changes landed: figure include widths + skip (layout only, no text), and the
"open-source"→"official" wording the supervisor itself recommended. Gates 1/2/4 are
text-invariant under those changes; Gate 3's numbers unchanged. Final-state build and
renders (pages 101 and 129 printed) pass. No upload has been made — gates re-run
before any upload per the supervisor's condition.

## Files for Ben

Changed (durable copies beside this report; live edits in
`tmp\overleaf_live_20261023\`):
- `chapters/30-torque.tex` (U14)
- `chapters/60-sns.tex` (U21)
- New asset for Overleaf `figs/30-torque/`: `Flx20mm_meetsHuman.pdf`
  (source of record: `Documentation\Reports and Papers\Knee_Torque_Test\Figures\`)

Verification PDF: `D:\temp\diss_u14u21_20261023\main.pdf` (256 pp).
Page renders: `D:\temp\diss_u14u21_20261023\page120_final.png`,
`page141.png` (Table 6.1, from the pre-wording pass — cell text identical),
`page148_final.png`.

## Open / for Ben

1. Layout choice for Fig 4.11 (see caveat above) — current 0.72/0.30 is minimal-churn;
   a three-panel regeneration from FlxGrp.fig is available on your word.
2. Upload needs your explicit go (working copy is the fresh zip; nothing has been
   uploaded; no commits/pushes made).
3. `fig:sns_curriculum` caption says "five-stage"; the curriculum code family documents
   stages 1–3 (air-deaff/air-aff/ground). The figure predates this pass and was
   conserved by earlier sweeps — flagged, not changed.
