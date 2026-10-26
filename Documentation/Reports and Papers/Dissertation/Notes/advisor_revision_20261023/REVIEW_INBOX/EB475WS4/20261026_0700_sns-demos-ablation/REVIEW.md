# Review batch

- Machine: EB475WS4
- Repository commit / uncommitted source version: HEAD 025e3897 (Ben's
  2026-10-26 merge/pull, branch KneeTestSetup_BenBo_stw). Baseline
  chapter snapshot for anchors: `Revision_D\chapters\60-sns.tex` AS
  COMMITTED (includes the easteregg2 polish pass: knowledge-gap table +
  sec:w2l_perturbation). Uncommitted sources used: figure assets and
  results under `Code\Matlab\SNS_Simscape\results\` (generation below).
- Purpose: deliver (1) replacement demo figures for the SNS demos
  subsection — Ben's 2026-10-26 rejections addressed: leg present,
  cyclic response, cross-software comparison, two-state beer mechanism
  with visible muscles, postsynaptic markers, elevation sign convention;
  (2) the new feedback-ablation pair result (a newly checked result, not
  author-confirmed) as figure + evidence.

## Patch manifest

| Fragment | Target chapter / section label | Exact existing start/end text | Advisor comment author/date/quote | Change and reason |
|---|---|---|---|---|
| P01_sns_demo_figures.tex | ch. 60-sns, insertion after the figure env ending `\label{fig:simulink_reflex_topology}` (sec:sns_simulink_demos) and after `\label{fig:w2l_perturbation}`'s env (sec:w2l_perturbation) | No text replaced — pure insertions; anchor lines are `\end{figure}` closing those two floats in the named snapshot | Ben (session directives 2026-10-26): knee reflex "no leg ... no comparison to the demo response vs literature response (or other software's response)"; CPG "no leg model doing any movement"; beer "muscles ... empty cup and a full beer cup"; markers "at the end of the synapse"; sign: figure "appear[s] ... the glass is raised higher" | PROPOSAL: three demo figure floats + one ablation float with captions; replaces no existing figure (existing `fig:simulink_reflex_topology` retained) |

## Figure manifest

| Asset | Intended number / label | Replace or add | Caption and panel identities | Generator and data |
|---|---|---|---|---|
| fig_kneereflex_v4.pdf/.png | fig:sns_demo_kneereflex | ADD | A kinematics + numpy twin; B leg mechanism + BPA routing; C reflex pathway (markers at endings, afferents dashed) | `demos\sns_demo_figures_v4.m`; data `results\sns_demo_results.mat` + `results\knee_twin_ref.mat` (twin: `dev\knee_reflex_numpy_twin_20261026.py`, exact constants from the builders) |
| fig_kneecpg_v4.pdf/.png | fig:sns_demo_kneecpg | ADD | A cyclic window; B one-cycle leg strip; C architecture + four-corner core | same generator; `results\sns_cpg_results.mat` |
| fig_beer_v4.pdf/.png | fig:sns_demo_beer | ADD | B1 empty-cup start, B2 full-cup end, A elevation (0 = level, negative below; ON −2.00° above OFF −7.85°), C pathway | same generator; `results\sns_beer_results.mat` (runs 1=ON, 2=OFF) |
| ablation_pair_metrics_20261026.pdf/.png | fig:sns_ablation_pair | ADD | bar panels: cycles L (33/1), cycles R (28/0), forward (1.16/0.09 m); intact vs ablated | `dev\plot_ablation_metrics_20261026.py` from `results\ablation_{intact,ablated}_20261026.json` |

## Evidence

- Model / controller / effective parameters: W2L corrected connectome
  (`bilateralrg2`, census 96/204) in MuJoCo 2.3.7 via the wean-campaign
  evaluator; frozen knobs both conditions: te 3.0, tf 4.0, comm 1.0,
  ctn_amp 4.0, cref 50.0, ibnA 1.0, ch_hc 0.08, ch_ank 1.0,
  ch_contra 0.08; protocol legcoll 0, sink 0.06, rig 0.25,
  walker_yield 1.5. Ablated condition zeroes ch_hc, ch_ank, ch_contra,
  ctn_amp, ibnA ONLY.
- Protocol / support / stimuli / perturbations: 30 s fresh-process evals
  per condition; harness/rig protocol identical; no retuning in either
  condition. Predictions pre-registered in
  `SNS_Simscape\reports_20261024\PREREG_REPRO_ABLATION_20261026.md`.
- Results and definitions: survive=true/false and objective from the
  campaign metric (verbatim `_w2l_wean_eval`); cycles = primary heel
  strikes; forward = pelvis travel. INTACT: obj 90.0, 33/28 cycles,
  phase 0.492 (n=30), +1.163 m, pelvis min 0.905 m. ABLATED: obj 5.8,
  1/0 cycles, phase undefined, +0.091 m. Prediction P5 (phase degrades
  first, cycling continues) was WRONG — cycling itself collapses.
- Saved output paths: `SNS_Simscape\results\ablation_{intact,ablated}
  _20261026.json`; CSV rows tagged `ablation_*_20261026` in
  `spinal\reports_20261016\w2l_wean\wean_trials.csv`; full report
  `SNS_Simscape\results\ABLATION_PAIR_REPORT_20261026.md`. Traces npz
  did NOT persist (harness redirect) — metrics only; rerun available.
- Scope / author-confirmed facts: Ia/II chains are not harness-knobbed
  and remain at template defaults in BOTH conditions — this is a
  contact+Ib ablation, not total deafferentation (stated in the report).
  Beer/knee/CPG demos use placeholder rig parameters (documented); the
  numpy twin is a cross-software check, not a literature comparison.

## Validation

- Build or checks performed: figures generated headless R2025a today;
  no LaTeX build run (insertions not applied; floats reference assets
  by relative name assuming figures/ co-location per project
  convention). Ablation runs are the saved builds of record (JSONs +
  CSV rows above).
- Known reference/asset/layout/interpretation issues: anchors taken
  from the committed Revision_D snapshot — RE-VERIFY against live
  Overleaf (application baseline) before applying; existing
  `fig:simulink_reflex_topology` kept, so the subsection gains floats —
  layout/page budget to be checked at application; ablation figure is
  metrics bars (traces not persisted); knee demo comparison is vs our
  numpy twin (other software), explicitly NOT a literature-response
  comparison.

## Preserved author decisions

- Existing tracked deletions, replies, wording, and figure identities
  retained: `fig:simulink_reflex_topology` and `KneeReflex_circuit.pdf`
  untouched; the easteregg2 polish-pass prose (knowledge-gap table,
  sec:w2l_perturbation) untouched; no tracked content reintroduced;
  reviewed chapters not modified (patches are proposals only).
