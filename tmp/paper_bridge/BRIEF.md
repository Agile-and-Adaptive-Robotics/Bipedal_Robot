# BRIEF — Isometric Torque Paper bridge session (read this FIRST)

Working context for every agent in this run. Paths are workspace-relative
(workspace root = `D:\GitHub\Bipedal_Robot`). Machine = easteregg2.

## Deliverable

Journal manuscript for **Frontiers in AI and Robotics**, working title
**"Design Factors for Isometric Robotic Knee Torque Using Braided Pneumatic
Actuators"** (short/running title: "Design Factors for Isometric Robotic Knee
Torque"). The local project mirror is
`Documentation/Reports and Papers/Knee_Torque_Test/ProofFinal/`
(monolithic `Main.tex`, `FrontiersinHarvard` class, natbib `\citet/\citep`,
`BolenFrontiers22.bib`, siunitx loaded, figures in `ProofFinal/Figures/`).
The current draft carries an older title and an early-register abstract/intro;
the current PDF/zip in the parent folder show the state Ben exported from
Overleaf on 2026-10-17.

## Strict folder rules (non-negotiable)

- `.tex` edits go ONLY into `Knee_Torque_Test/ProofFinal/`.
- New `.bib` entries go ONLY into `Knee_Torque_Test/Bibliography/BolenFrontiers22.bib`
  (and the mirrored copy at `ProofFinal/BolenFrontiers22.bib` stays in sync).
- Tables are created in the `.tex` AND exported as `.xlsx` into `Knee_Torque_Test/Tables/`.
- Do NOT create, move, rename, or restructure anything else inside
  `Knee_Torque_Test/` (leave the PDFs, zip, jpg, Gantt xlsx, and the local
  `Figures/`, `Figure components/` untouched — that folder is Ben's).
- READ-ONLY everywhere else: ALL dissertation sources (`Documentation/.../Dissertation/`),
  all code, all data, all recorded results. No git commits (Ben commits via
  GitHub Desktop). Scratch work goes ONLY in `tmp/paper_bridge/`.

## Editorial contract (journal register)

- Structure like the Actuators force paper (bolen_2026, *Actuators* 15(5):230;
  local source folder `Documentation/Reports and Papers/Isometric_Force_Test/`):
  Introduction, Methods and Materials, Results, a LONG Discussion, Conclusion.
- The LONG Discussion carries the follow-up system-identification results
  (Xi0–Xi3 correction factors, two-bracket flexor/extensor identification) and
  the route-redesign results (mesh-optimization attachment redesign meeting
  human isometric torque targets) as ADDITIONAL RESULTS.
- Fictive-strain framing: correction terms are expressed through the unloaded
  620 kPa contraction epsilon_620 (the "fictive strain" of the free actuator),
  as the dissertation chapter does.
- NEVER mention campaign dates, file names, pick/trial numbers, or tooling
  session artifacts. The paper reads as if it always knew.
- Oxford comma, always. American spelling. Past tense for procedures; present
  tense for model behavior and interpretation. "We/our" for this work.
- Units through siunitx (`\qty`, `\qtylist`, `\SIrange`). No bare "10mm".
- Diameter: engineering symbol `$\diameter$` (add `\usepackage{mathabx}` to the
  preamble if the current preamble lacks it; flag this in your notes).
- Banned AI-slop: "delve", "it is important to note", "plays a crucial role",
  metaphorical "landscape", moreover/furthermore/additionally chains, em-dash
  padding (max ~1 per 500 words, never clause-splicing), uniform paragraph
  rhythm. Calibrate voice against
  `Documentation/.../Dissertation/Notes/revision_campaign_20261015/style-contract.md`.
- Every number needs a named source (results ledger row, dissertation chapter
  location, or cited reference). No source, no number. Never harmonize
  disagreeing numbers silently.
- Call the physical leg/test articles by what they are (robot leg, test stand,
  robot model). The word "plant" is banned for the robot or its dynamics.
- Do NOT change the author list, affiliations, or correspondence block.
  Flag (do not fix): `Main.tex` currently lists 4 authors while the
  dissertation publication note (30-torque.tex:5) names 7 contributors for
  bolen_isometric_2025 — Ben's call.
- Citations: only keys already present in `ProofFinal/BolenFrontiers22.bib`.
  A missing needed citation gets REPORTED (claim + suggested source), never
  invented, never silently dropped.

## Source material

- Primary text source (READ-ONLY, committed state = 0d23da03, 2026-10-17;
  chapter still under revision until 2026-10-21 — do not reconcile against
  uncommitted edits, there are none right now):
  `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/30-torque.tex`
  with figures in sibling `Figures/30-torque/`. Its publication note names the
  manuscript-overlap sections (bolen_isometric_2025).
- Archived test results (READ-ONLY): `Testing_Data/2022_02_Festo/`
  (minimizer family, `Dig_out/` result mats + logs, plot scripts, `.fig`
  exemplars), `Testing_Data/2026_06_Festo/` (Results_table xlsx tabs — read
  the tables, not the long raw files), `Code/Matlab/Mesh_Optimization/Results/`
  (route-redesign), `Code/Matlab/Robot_Data/` + `Code/Matlab/Functions/`
  (MonoPam classes/helpers).
- Values of record prefer the dissertation chapter and its appendices
  (committed state) over any individual run artifact; disagreements get
  flagged, not averaged.

## Figure specification (for regenerated result figures)

- Master format exemplar: `Testing_Data/2022_02_Festo/FlxPin_group.fig`.
- Script conventions: `Plot_KneeFlx_all.m` (Arial, FontSize 12 bold axes,
  LineWidth 2, minor ticks both axes, TickLength [0.025 0.05], legend 8 pt
  SouthWest, XLim [-120 15], XTick [-120 -90 -60 -30 0]) and
  `Plot_KneeFlxPin10mm_42cm.m` (hex palette list c{1..8}, gold #FFD700 first).
- PlotPrompt (decoded copies: `tmp/paper_bridge/PlotPrompt_rtf.txt` and
  `PlotPrompt_utf8.txt`): measured = indigo filled scatter (MarkerFaceAlpha
  0.75), hybrid/corrected model = gold c{1}, original uncorrected model =
  gray dashed [0.4 0.4 0.4] LineWidth 2, corrected prediction = c{5}
  LineWidth 2.5; tiled panels labeled (A)–(D) at 12 pt bold Arial with
  EdgeColor none; shared axes labels; GridLineStyle none.
- QC rules from `revision_campaign_20261015/figure_qc.md`: every caption
  states "negative values denote knee flexion"; fonts must stay legible at
  print scale (~6.2 in width); no hue-only legend separation (vary dash or
  marker too); labeled, unit-bearing axes.
- New/refreshed figures land in `ProofFinal/Figures/` (same filename when
  replacing a stale figure so `\includegraphics` stays stable; EPS like the
  existing set). Plot scripts live in `tmp/paper_bridge/figscripts/` and are
  run headless with `matlab -batch` (MATLAB R2025a is on PATH; addpath
  `Code/Matlab/Functions` and `Code/Matlab/Robot_Data` as the existing plot
  scripts require; load data from `Testing_Data/2022_02_Festo`).

## Toolchain facts (easteregg2)

- Python: `D:/Anaconda/python.exe` (openpyxl 3.0.10 available) — use for
  `.xlsx` table exports and any scripted checks. Call by full path.
- MATLAB: `matlab -batch "..."` (R2025a).
- TeX (installing): `D:/MiKTeX/miktex/bin/x64/pdflatex.exe` and `bibtex.exe`.
  May not exist yet — treat missing binaries as "compile pending", never as a
  manuscript failure. The project also carries `frontiers_SupplementaryMaterial.tex`
  (`\externaldocument` via xr) — compile it once before Main if xr complains
  about a missing aux.
- Git: available, but NO commits, NO branches, NO pushes.
