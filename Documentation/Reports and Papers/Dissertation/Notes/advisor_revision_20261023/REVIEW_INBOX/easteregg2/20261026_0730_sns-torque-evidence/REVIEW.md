# Review batch

- Machine: easteregg2
- Repository commit / uncommitted source version: d19d4ae2 (working tree
  clean). Prose source snapshot: **Revision_D working copy of 2026-10-23**
  (`D:\GitHub\_zcode_backups\diss_context_20261024\Revision_D\`, from the
  EB475WS4 zip), which contains this session's edits; the Overleaf
  baseline has since moved to Revision_E, so every boundary below must be
  re-checked against current Overleaf text before application. All
  bounded texts are quoted from the Revision_D snapshot.
- Purpose: deliver this session's SNS-chapter evidence additions and
  torque-chapter corrections as bounded proposals: the perturbation/
  pathway-deletion subsection + figure, the capability-table rewrite per
  Ben's rubric, the conversion-checks paragraph answering the advisor's
  2026-10-23 question, the pulley-section correction, and the torque
  figure set (4.9/4.10 merge, 4.11C shared axes, 4.13/4.14 journal
  regeneration without captions).

## Patch manifest

| Fragment | Target chapter / section label | Exact existing start/end text | Advisor comment author/date/quote | Change and reason |
|---|---|---|---|---|
| P01_sns_conversion_paragraph.tex | 60-sns, `\subsection{Model Conversion and Actuator Coupling}\label{sec:sns_existing_tools}` ¶1 | STARTS "MyoConverter supplies the initial OpenSim-to-MuJoCo translation. I added procedures…" ENDS "…exchange signals at a \qty{2}{\milli\second} step." (¶2 "For the BPA-routed robot model…" untouched) | Advisor 2026-10-23 11:31 (asked what problems the conversion checks actually found and corrected, rather than implying defects by calling the procedure an audit); Ben Overleaf replies 14:03/14:05/14:07 (pruning, ISB frame, patella, robotbody-first import) | Newly checked result: names the five concrete problems found/fixed (robotbody infeasible-parameter failure and human-first sequence; 3-per-side rotators rendered inert; axis convention preserved via body-frame remap; patella follower limits + rect_fem reroute; contact-pair artifact). Evidence dossier: `_zcode_backups\diss_campaign_20261024\B1_conversion_results.md` |
| P02_sns_capability_table.tex | 60-sns, Table `tab:sns_capability` (four-row "Model/interface" table) | STARTS `\begin{table}[htbp]` containing caption "Implemented tools, representative checks, and limits." ENDS its `\end{table}`; the two framing sentences at "Table~\ref{tab:sns_capability} distinguishes…" (§ sns_gaps) and "…separates the implemented interfaces…" (§ preliminary_sim_methods) are part of the patch | Ben 2026-10-26 (mobile): "Table 6.1 sucks. We don't need the last column that makes it look like i have incomplete work. I spelled out what was missing and it should be in the chatgpt report. You're trying to make the old table work instead of rewriting it." | Table rewrite per the rubric recorded in CHATGPT_REPORT.md (gap / why it mattered / tool developed / demonstrated outcome — no deficiency column). Row 5 follows the current editorial ruling: the coupled SolidWorks-derived robot with neural control and BPAs is author-confirmed (saved model on the unavailable laptop; no invented measurements). |
| P03_sns_perturbation_subsection.tex | 60-sns, new subsection inserted before `\section{Capabilities and Remaining Limits}\label{sec:sns_discussion}` | INSERTION POINT: immediately after the `\subsection{Closed-Loop Model Coupling}` paragraph ending "…remains a limit on transferring tuned controller parameters between implementations." | Ben 2026-10-25/26 directives: perturbation/deletion evidence on a *neuromechanical* model, Deng-poster figure vocabulary (stimulus to RG/PF, external force, pathway deletions) | Newly checked result: perturbation-robustness + pathway-deletion subsection with the corrected (fully-wired) numbers; see Evidence below for the two-audit correction trail. Prose-only claims scoped: supported progression, no unsupported walking. |
| P04_torque_pulley_section.tex | 30-torque, `\subsection{Mono- and Bi-Pulley Tendon Routing}\label{sec:disc_pulley}` ¶1 and the "subsequent transmission study" sentence | STARTS "At unity transmission, the knee-flexor pulley route increased…" ENDS "…Physical verification of the pulley configurations remains future work." (whole subsection body; the following `\subsection{Physical Tests…}` untouched) | Ben 2026-10-26 (mobile): "Torque section mono & bi-pulley talks about unity transmission which is nonsense. The point of the reverse pulley is to trade force for range of motion. The mono and bi just signify if there is one or two dominant DoF that are being calculated." | Author-confirmed fact (terminology/semantics): "unity transmission" removed (only the statistical FVU-below-unity usages elsewhere remain); reverse pulley stated as force-for-range-of-motion trade; mono/bi defined as one vs two dominant DoF; "preserves net actuator force" corrected. |
| P05_torque_fig4_9_10_merge.tex | 30-torque, figures `fig:Ext10mm_52cm` and `fig:Ext_passthru` | STARTS `\begin{figure}[htbp]` containing `\includegraphics[width=1\textwidth]{Ext10mm_52cm.eps}` ENDS the `\end{figure}` of the `KneeExt_10mm_passthru.png` figure; the in-text sentence "Fig.\ref{fig:Ext_passthru} shows the difficulty…" is part of the patch | Ben 2026-10-26 (mobile): "Figure 4.9 should merge with 4.10." | Merge into one (A)/(B) figure: (A) torque plot, (B) routing photo; `fig:Ext_passthru` label retired, its in-text reference repointed to panel B of the merged figure. |
| P06_torque_figs_4_13_14.tex | 30-torque, figures `fig:extensor_physical_final` and `fig:flexor_physical_final` | STARTS `\begin{figure}[p]` containing `\includegraphics[width=0.97\textwidth]{extensor_20mm_physical.png}` ENDS the `\end{figure}` of the flexor block | Ben 2026-10-26 (mobile): "Figure 4.13 and 4.14 don't need captions. Make sure they are journal level quality." | Captions removed (labels retained for cross-references); both figures regenerated as vector PDFs in a matched journal style from Ben's unmodified analysis scripts (see figures/). LaTeX switches from .png to .pdf assets. |

## Figure manifest

| Asset | Intended number / label | Replace or add | Caption and panel identities | Generator and data |
|---|---|---|---|---|
| w2l_perturbation_panels.pdf/.png/.alt.txt | new figure, `fig:w2l_perturbation` (SNS ch.; number per final placement) | Add | (A) control knee flexion L/R + heel-contact shading; (B,C) 1-s inhibitory tonic stimulus to the rhythm layer, mild/strong; (D,E) 40 N fore/aft pelvis pushes at t=8 s (COM x + tilt); (F) strides-per-leg table for six pathway deletions. Knee flexion positive; knee angle MuJoCo sign; forces N. | `spinal\reports_20261025_persyn\_poster_panels.py` (committed, corrected) + forwarding-fixed `test_w2l_ground.py`; raw logs `rawlog_*.npz` in the same folder; deterministic across two fresh processes |
| extensor_20mm_physical.pdf/.png | ch.4 `fig:extensor_physical_final` (was 4.13) | Replace asset | Caption removed per directive; content: measured (indigo) vs back-calculated (gold) vs optimized 620 kPa (magenta solid) vs human target (dotted); marker area ∝ pressure | `Mesh_Optimization\Results\Knee_Extensor_20mm.m` (unmodified) + journal wrapper; 600 dpi PNG + vector PDF |
| flexor_20mm_physical.pdf/.png | ch.4 `fig:flexor_physical_final` (was 4.14) | Replace asset | Same conventions as extensor twin; identical axis window (−125…10 deg), matched tick styling | `Knee_Flexor_Data_20mm.m` (unmodified) + wrapper |
| FlxGrp_20261023.eps | ch.4 `fig:FlxGrp` panel C (was 4.11C) | Replace asset | Panel C now shares panel B's y-range and x-range (both 20 mm flexor vs the same human target) and carries the matching legend; panel B legend white-backed at northeast | `tmp\paper_bridge\figscripts\gen_FlxGrp_3panel_20261023.m` (corrected, committed); QC png included as preview |

## Evidence

- Model / controller / effective parameters: corrected W2L biped
  (bilateralrg2, 96/204) in MuJoCo 2.3.7 at 2 ms; drive pinned te=3/tf=4;
  round_w25.json parameter set with all per-target/per-synapse keys
  actually applied (post forwarding fix — see correction trail below);
  harness wsup 0.25 + sagittal constraint; stimulus 1.5/3.0 nA 1 s on the
  RG drive ports; push 40 N, 0.15 s, pelvis.
- Protocol / support / stimuli / perturbations: as in P03; deletions =
  zeroed pathway classes (crossed heel, Ia chain, Renshaw, heel cascade,
  commissural, Ib afferent).
- Results and definitions: primary heel-strike strides (>20 N rising
  edge, 0.45 s merge) per leg over a 16 s window. Control 17L/16R; stim
  mild 18/17, strong 19/17; push fwd 18/18, back 18/16; deletions
  16/16, 13/13, 16/14, 17/15, 13/12, Ib 0/0. Mass-change study
  (33–38 strides/leg, pelvis ≥0.90 m in all 10 conditions; 9/10 exceed
  the 45 deg tilt gate) reported with the posture caveat. Wean ladder
  70→10 % with fresh verification per rung (6-family-tuned
  configuration; see correction trail).
- Saved output paths: `spinal\reports_20261025_persyn\`
  (deletion_results.json, changed_body_test_w25.json [6-family config],
  changed_body_TRUEWINNER_w25.json, rawlog_*.npz, round_w*.json);
  `spinal\reports_20261024_condsolve\` (class-level study + verdict);
  supervisor logs in `_zcode_backups\diss_campaign_20261024\`.
- Correction trail (material to interpreting the numbers): two audits
  corrected this evidence set. (1) The panel script originally duplicated
  left-side stride counts into the right column and mislabeled a stored
  knee range (external audit 2026-10-26; fixed, rerun, raw logs saved).
  (2) The evaluation path originally discarded 55 of 61 conductance keys
  and the runner did not forward per-target gains, so the first
  per-synapse search was effectively 6-dimensional and the first
  Ia/Renshaw/heel "deletions" were no-ops (external audit 2026-10-26;
  runner fixed; committed panel script rerun; deterministic). A full
  per-synapse search redo with live forwarding is RUNNING at delivery
  time and is not part of this batch's claims.
- Scope / author-confirmed facts: pulley semantics (P04) are Ben's
  rulings; the coupled SolidWorks robot row (P02) is author-confirmed per
  the inbox's current editorial rulings; Rybak claims elsewhere remain
  scoped to the Python implementation.

## Validation

- Build or checks performed: the fully patched Revision_D working copy
  compiles clean (fresh `build_polish_20261026`, MiKTeX full path,
  absolute BIBINPUTS/BSTINPUTS): 197 pages, 0 figure errors, 25 overfull
  boxes (unchanged from baseline). Figures passed the visual-judge gates
  after repair rounds (panel-F table clipping; minor-tick mismatch;
  legend occlusion — all recorded in
  `_zcode_backups\diss_campaign_20261024\SUPERVISOR_polish_20261026.md`).
  This is a saved build of the patched working copy, NOT of Overleaf.
- Known issues: boundaries must be re-verified against Revision_E/
  Overleaf (the snapshot predates it); the editorial ruling on the
  "suspended stepping result" tracked deletion and the live capability
  table is respected (P02 replaces the table's content, it does not
  remove the table); intended figure numbers shift with final placement.

## Preserved author decisions

- Existing tracked deletions, replies, wording, and figure identities
  retained: no tracked deletion is reintroduced or paraphrased; the
  Connor–Ben shared-RG diagram identity is untouched; Ben's pulley
  terminology rules are applied verbatim in P04; figure 6.1(B) original
  is not substituted; no Overleaf edits, no comment resolutions, no
  commits/pushes performed for this delivery.
