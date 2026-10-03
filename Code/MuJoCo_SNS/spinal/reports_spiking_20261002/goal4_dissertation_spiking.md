# Goal 4 — dissertation update: spiking-neuron versions of the models (2026-10-02)

Host easteregg2, repo `D:\GitHub\Bipedal_Robot`, branch `KneeTestSetup_BenBo_stw`.
Working copy of record: `Documentation\Reports and Papers\Dissertation\ProofFinal\`
(+ the Dissertation `Figures\` tree). No git state was changed (status/diff only);
no Overleaf interaction (owned by the other chat per the ask).

## 1. What was edited (complete list, additive only — `git diff` = 283 insertions, 0 deletions)

| File (ProofFinal\chapters\) | Edit |
|---|---|
| `20-methods.tex` | NEW subsection `\subsection{Spiking-Neuron Mirrors of the Simulation Models}\label{sec:spiking_methods}` appended at end of file (after the `simulink_reflex_topology` sidewaysfigure, inside section `sec:preliminary_sim_methods`): hybrid doctrine, the MuJoCo mirror build, the script-calibration approach, the three Simulink spiking blocks + twin method, the AnimatLab conversion procedure, and the gate set. 6 fenced blocks marked `% [ZCODE 2026-10-02 Spiking-mirror campaign BEGIN ... END]`. |
| `30-results.tex` | NEW section `\section{Spiking-Neuron Mirror Results}\label{sec:spiking_results}` appended at end of file with three subsections (MuJoCo / Simulink / AnimatLab), 3 tables (`tab:spike_ground`, `tab:spike_knee`, `tab:spike_beer`, `tab:spike_animatlab`) and 3 figures (`fig:spike_rhythm`, `fig:spike_ground`, `fig:spike_knee`). |
| `40-discussion.tex` | NEW section `\section{What Spiking Changes and Does Not}\label{sec:disc_spiking}` appended after `sec:synthesis`: plant-side insulation under the hybrid doctrine, period/cadence deltas with identified knobs, basin-robustness gain vs the analog build's near-bifurcation marginality, the rate-code re-tune requirement (AnimatLab negative + ablation), and the bilateral-gait lead flagged as one run with no statistics. |
| `93-AppendixB.tex` | NEW section `\section{Spiking-Mirror Implementation, Calibration, and Verification}\label{app:spiking_details}` appended at end of file: adapting-LIF equation `(eq:app_lif)` + per-spike conductance equation `(eq:app_spk_syn)`, cell-parameter table `tab:app_spk_cells`, calibration-of-record table `tab:app_spk_cal`, LANDSCAPE full gate table `tab:app_spk_gates` (sidewaystable), LANDSCAPE AnimatLab variant table `tab:app_spk_animatlab_full` (sidewaystable), conversion rules + the per-connexion `<Type>` trap, extra figure `fig:spike_animatlab`, and the honest scope notes (no GUI verification of the .asim copies; the two deliberately not-built conversions). |

NOT touched: `02-abstract.tex`, `10-introduction.tex`, `60-conclusion.tex` (left out on
purpose — the spiking results are methods/results/discussion material; nothing in them
rises to a claim the abstract or conclusion currently makes about the non-spiking work,
and the ask says when in doubt, leave out), `15-background.tex`, `50-futurework.tex`,
`92-AppendixA.tex`, `94-AppendixC.tex`, `main.tex`, `thesis.bib` (only existing citation
keys were used: `szczecinski_perspective_2023`, `nourse_sns_2023`).

## 2. Figures added (all in Dissertation figure dirs; pdf + png + `_alt.txt` beside each)

| File (Dissertation\Figures\) | Content | Data source |
|---|---|---|
| `30-results/spiking_mirror_rhythm.{pdf,png}` + `_alt.txt` | (A) non-spiking RG E−F readout (constant drive, 20 s), (B) spiking mirror raw RG_E/RG_F cumulative spike counts | `reports_spiking_20261002/goal4_gate2_traces.npz`, dumped this session by `tools/goal4_dump_gate2_traces.py` (verbatim GATE-2 recipe; printed periods 0.399/0.621 s match `logs/gate2_rhythm2.log` exactly — verified this session) |
| `30-results/spiking_mirror_ground_gait.{pdf,png}` + `_alt.txt` | knee + hip r/l, non-spiking s3k winner vs untuned spiking mirror (16 s ground runs) | `spinal/spinal_run_spkbase_s3k.npz` + `spinal/spinal_run_spkbase_spiking.npz` (goal-1 artifacts, no re-run) |
| `30-results/simulink_knee_spiking_twin.{pdf,png}` + `_alt.txt` | KneeReflexDemo vs KneeReflexDemo_Spiking: full 5 s + settle window | `goal4_knee_traces.mat`, exported this session via MATLAB R2025a (`logs/goal4_knee_export.log`; printed rise/settle values match `goal2_knee_spiking.mat` exactly) |
| `93-AppendixB/animatlab_spiking_conversion.{pdf,png}` + `_alt.txt` | AnimatLab W2L hip+knee: graded baseline / spiking conversion / thresholds-only ablation (new figure dir created) | chart `.txt` byproducts `models/W2L_modern/`, `models/W2L_modern_spiking/`, `models/bisect/` (ranges verified against `metrics_bisect_W2L_neuronsonly.json` etc.) |

Figure standards met and how it was checked:
- **7.5×10 in usable area**: all four figures are 7.5×6.0 / 7.5×5.2 / 7.5×5.2 / 7.5×7.0 in
  (measured from the 300-dpi PNGs: 2250×1800, 2250×1560, 2250×1560, 2250×2100 px).
- **Arial ≥10 pt, no italics**: PDF font inventory is exactly `ArialMT` with
  `/ItalicAngle 0` for every figure (regex over the embedded font descriptors);
  no mathtext is used in any label.
- **Colors.m Paul Tol palette** (`Code/Matlab/Colors.m`): indigo `#0000FF` (left /
  baseline / non-spiking) and orange `#FFB14E` (right / twin / spiking), pink `#EA5F94`
  for the third variant — solid vs dashed as line-style redundancy, matching the
  existing `animatlab_phase1_joint_angles` convention.
- **Alt-text file beside every figure**: `_alt.txt` written next to each pdf/png AND
  `% ALT TEXT:` comment pointers left at each `\includegraphics` (the xi-figure
  precedent); every number in the alt texts was verified against the plotted arrays
  (two alt texts were corrected during this pass when the data showed my first
  wording was wrong — rhythm panel A swing is 2.3–7.1 mV sustained, not 0–8.3).
- **Visual pass**: run personally on all four figures via the remote vision backend
  (per-figure inspection: legibility, clipping, trace distinguishability, content
  match). The primary local vision plugin returned HTTP 401 (down), so the pass used
  the alternate backend + programmatic checks (ink-margin bounds ≥41 px on all four
  sides of every figure, font inventory, italic scan). Two real findings were fixed:
  the rhythm caption/alt-text misdescribed the sustained swing, and the AnimatLab
  caption now states the 0.5–5 s chart window (the charts genuinely start at 0.5 s).
- **Shape-coded synapses**: not applicable — these are trace/data figures; no circuit
  diagram was drawn.

## 3. Number provenance (every number in the new tex traces to a gated artifact)

- MuJoCo family: `goal1_walkers_spiking.md` §3/§5/§7/§9 + `logs/baseline_s3k{,_postedit}.log`,
  `logs/gate2_rhythm2.log`, `basin_gate_spiking.json`, `logs/gate3_basin.log`,
  `air_smoke_{ns,spiking}.json`, `logs/ground_eval_spiking.log`, `spiking_calibration.json`.
- Simulink family: `goal2_sns_simscape_spiking.md` §1–§4 + `goal2_baselines.mat`,
  `goal2_knee_spiking.mat`, `goal2_beer_spiking.mat`, `goal2_rgmn_spiking.mat`
  (all re-read this session; knee rise/settle re-verified by the fresh MATLAB export).
- AnimatLab family: `goal3_animatlab_spiking.md` §3–§5 + `metrics_baseline_*.json`,
  `metrics_*_spiking*.json`, `metrics_bisect_W2L_neuronsonly.json`
  (all re-summarized this session from the jsons directly).
- All headline results stated in the tex are the SHIPPED (supervisor-gated) claims
  verbatim; the two gate FAILs (period +55.6 %, air cadence) and the one PARTIAL
  (RG→MN alternation) are stated as failures/partial in the tex, matching the reports.
- Two numbers in my tex come from the goal-2 supervisor's audit rather than the goal-2
  report body, both corrections in my favor of accuracy: the library is **15** blocks
  after the additive change (report says 14, audit corrected), and the AnimatLab
  baseline metric files number **5** (report §3 says 7, §8 says ×5).

## 4. Checks run this session (exact commands)

1. MATLAB trace export (goal-4 figure data):
   `"D:\Program Files\MATLAB\R2025a\bin\matlab.exe" -batch "run('D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal4_export_knee_traces.m')"`
   → `logs/goal4_knee_export.log`: base rise 15.79, settle 43.5 (41.6–44.9); spk rise
   15.79, settle 44.4 (41.9–46.7) — identical to `goal2_knee_spiking.mat`. PASS.
2. Gate-2 trace dump (figure data + gate re-verification):
   `cd Code\MuJoCo_SNS\spinal && CONDA_PREFIX=D:\Anaconda\envs\myo D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\goal4_dump_gate2_traces.py`
   → `logs/goal4_gate2_traces.log`: non-spiking 0.399 ± 0.164 s, spiking 0.621 ± 0.002 s
   (E 4.9 / F 3.2 Hz), +55.6 % — every line matches `logs/gate2_rhythm2.log`. PASS
   (the trace data behind `spiking_mirror_rhythm` is gate-verified).
3. Figure QA (programmatic): inline PIL/numpy/re script over the four PNGs+PDFs —
   sizes, ink margins, font inventory `ArialMT`, `/ItalicAngle 0`. PASS.
4. Visual pass: alternate vision backend, one call per figure (see §2). PASS after
   the two caption/alt-text corrections.
5. Static LaTeX checks: `D:\Anaconda\envs\myo\python.exe reports_spiking_20261002\tools\goal4_tex_check.py`
   → no unresolved `\ref`s, no citation keys outside `thesis.bib`, all referenced
   figure files exist, environments + braces balanced in all four edited files,
   no (R)/(TM)/lowercase matlab in the new blocks. PASS.
   (NOT RUN: an actual LaTeX compile — easteregg2 has no local TeX
   (`latexmk`/`pdflatex` absent; `D:\MiKTeX` is EB475WS4-only) and Overleaf is out of
   scope for this goal. The Overleaf chat must compile once after merging.)
6. `git status --porcelain` re-run at the end: my modifications are exactly the four
   chapter `.tex` files; new untracked files are the figure sets, this folder's
   `goal4_*` files, and the sibling goals' pre-existing changes (none touched).
   One stray file from a failed first shell invocation
   (`spinal/reports_spiking_20261002logsgoal4_gate2_traces.log`, mangled redirect
   path) was removed by me this session.

## 5. Honest scope notes / things left deliberately alone

- **Pre-existing defect found, NOT fixed (out of additive scope)**:
  `30-results.tex` contains `\section{Follow-Up Identification and Route
  Redesign}\label{sec:ongoing}` **twice** (lines ~331 and ~358) — the first carries
  superseded Xi values (+8.9 mm...), the second the ZCODE-2026-09-26 audited values
  (+3.94 mm...). LaTeX will warn about the duplicate label and typeset both sections.
  Ben/the Overleaf chat should delete the first block (the audited one is the
  superset), but deleting existing content was outside this goal's additive mandate.
- **Pre-existing staleness noted**: `20-methods.tex` still says "seven-block Simulink
  SNS library" and `93-AppendixB.tex` says the network is "406-neuron" — both predate
  the 2026-09-22 library redesign (12 blocks) and later network additions (410
  populations). My new text states the current counts (15 blocks; (410, 376, 1186))
  without rewriting the older sentences.
- **ProofFinal vs Overleaf divergence**: the CPG spinal-network sections
  (CPG_DISSERTATION_UPDATE_NOTES.md PART A/B/C) were pasted to Overleaf only, so
  ProofFinal has no tuning-curriculum section. My additions are written to stand
  alone in ProofFinal, referencing only `sec:preliminary_sim_methods`,
  `sec:preliminary_sim_results`, `sec:prior_walker`, `app:neuron_equations`, and
  `ch:futurework` — all of which exist in ProofFinal. When the Overleaf chat merges,
  the two additions are complementary, not conflicting (mine never references the
  curriculum sections).
- **No new tuning, no new studies** — per the ask; every "future work" flag in the
  tex (fresh curriculum, w2lvar/syn6 mirrors, full Simulink conversion, GUI check of
  the .asim copies) mirrors the gated reports' own scoping.
- The untuned bilateral ground gait is stated in the tex exactly as gated: "a single
  run with no statistics, reported as a lead... not as an established gait result."

## 6. Files created by this goal (all under `reports_spiking_20261002\` except the figures)

`goal4_export_knee_traces.m`, `goal4_knee_traces.mat`, `goal4_gate2_traces.npz`,
`tools/goal4_fig_style.py`, `tools/goal4_fig_ground_gait.py`,
`tools/goal4_fig_rhythm.py`, `tools/goal4_fig_knee_twin.py`,
`tools/goal4_fig_animatlab.py`, `tools/goal4_dump_gate2_traces.py`,
`tools/goal4_tex_check.py`, `logs/goal4_knee_export.log`,
`logs/goal4_gate2_traces.log`, this report — plus the four figure sets under
`Documentation\Reports and Papers\Dissertation\Figures\{30-results,93-AppendixB}\`.
