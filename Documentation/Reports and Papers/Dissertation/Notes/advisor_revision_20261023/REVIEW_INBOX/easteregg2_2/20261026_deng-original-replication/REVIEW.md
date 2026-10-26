# Review batch

- Machine: easteregg2 (this batch folder: `easteregg2_2`)
- Repository commit / uncommitted source version: HEAD `d19d4ae2` (merge of
  KneeTestSetup_BenBo_stw), working tree CLEAN for all of
  `Documentation/Reports and Papers/Dissertation/` (verified with
  `git status` before batch assembly). Quoted chapter text is from
  `ProofFinal/chapters/60-sns.tex` as it exists at this commit (file mtime
  2026-10-23 04:59, last chapter-touching commit `39e3f707`). Egg's copy may
  be older than the live Overleaf/Revision_E text: **the Overleaf text is the
  application baseline**; if the anchor sentences below have moved or changed
  there, re-anchor before applying and treat this batch's quotes as the
  easteregg2 snapshot only.
- Purpose: deliver a newly checked result from 2026-10-26 (Order B1): the
  ORIGINAL archived Deng 2019 AnimatLab project (ModelDB 264844) run headless
  in its original engine with the published stimulation protocol of Deng et
  al. 2023, extending the chapter's original-engine replication series. One
  prose insertion, one bibliography insertion, one proposed figure.

## Patch manifest

| Fragment | Target chapter / section label | Exact existing start/end text | Advisor comment author/date/quote | Change and reason |
|---|---|---|---|---|
| `patches/P01_deng_original_replication.tex` | Ch. 7 (`60-sns.tex`), `sec:preliminary_sim_results` | INSERT between the paragraph ending `"...reproducing the published behavior of the original model in its original engine."` and the paragraph beginning `"The replication series has a defined extension, and its current boundary is stated plainly."` (char offset 50672 in the snapshot file; the two paragraphs are currently separated by one blank line) | No open advisor thread names this passage. Governing ruling from the review round (REVIEW_INBOX README, "Current editorial rulings"): "Rybak's identified implementation is Python. **Scope each published-model claim to its actual implementation and demonstrated behavior.**" | NEWLY CHECKED RESULT, proposed addition. One paragraph reporting the original-project rerun and its response to the published 2 nA stimulation protocol, scoped: original engine (not the program's ports), antagonist-only archived project (not the 2023 biarticular model), period-level recovery numbers only. Reason: the chapter's replication series currently stops at the lineage's own archived models; the published original has now been run, and the response family matches what the 2023 paper itself reports for this predecessor. |
| `patches/P02_bib_deng_biarticular_2023.tex` | `Bibliography/thesis.bib` | INSERT a new entry near `deng_neuromechanical_2019` (after `deng_biomechanical_2022`) | none (supports P01's citation) | PROPOSED bib entry for the 2023 Living Machines paper. Venue verified by web search 2026-10-26 (Biomimetic and Biohybrid Systems / Living Machines 2023, Springer LNCS); **volume, pages, and DOI intentionally left blank for Ben's Zotero record** — do not invent them. |
| `patches/P03_figure_env_deng_original_perturbation.tex` (OPTIONAL) | Ch. 7, `sec:preliminary_sim_results`, after the P01 paragraph | NEW float after the P01 insertion point | none (supports the proposed figure) | PROPOSED figure environment (`fig:deng_original_perturbation`) with the full caption, siunitx units throughout. Apply only if the figure is accepted. |

No existing sentence is deleted or rewritten by this batch; both fragments are
pure insertions, and the referenced tracked deletions and replies in the
surrounding sections are untouched.

## Figure manifest

| Asset | Intended number / label | Replace or add | Caption and panel identities | Generator and data |
|---|---|---|---|---|
| `figures/deng_original_perturbation_20261026.pdf` + `.png` + `.caption.txt` + `.alt.txt` | proposed NEW figure, suggested label `fig:deng_original_perturbation`, Section `sec:preliminary_sim_results` next to P01 | ADD (no existing figure replaced) | Full caption in `.caption.txt`. Panels: (A) left RG flexor, stimulus at stimulus onset ("rapid short step"); (B) left RG extensor, same run, rhythm continues; (C) left RG extensor under the knee-PF stimulus, indistinguishable from baseline; (D) left knee extensor motoneuron under the knee-PF stimulus held 2.0–2.5 s, local suppression. Units: mV vs s; voltages clipped at 0 mV by the integrate-and-fire element; +current = depolarizing. Panel letters are margin letters left of each tile per house convention; Arial; Tol-bright blue/red for baseline/stimulated. | `Code/MuJoCo_SNS/spinal/reports_20261026_dengoriginal/make_figure_review_20261026.py` over the saved chart files `deng_rg_mn.txt` (runs `run_B0_baseline`, `run_P1b_rgflexor_2p1`, `run_P2b_kneepfflexor_2p1`, `run_P2_kneepfflexor`) |
| `figures/REFERENCE_ONLY_twin_suppression_eb475ws4.png` | none | REFERENCE ONLY — not proposed for the manuscript | The Simulink twin's banked perturbation result (eb475ws4, 2026-10-26), included because the B1 comparison suggests its "extensor suppression" reading is plausibly the same phase-delay family seen through a 0.5 s window on a 1.94 s cycle. No manuscript claim is made from it here. | eb475ws4 `Code/Matlab/SNS_Simscape/results/pictures/deng_perturbation_twin_20261026.png` |

### Revision history of this batch

- REV 1 (2026-10-26, first pass): as described above.
- REV 2 (2026-10-26, after supervisor round 1 — overall BLOCK on two
  undispositioned claims): **F1 fixed** — the rapid-short-step/rhythm-
  continuation attribution now names the paper's biarticular model (its Fig. 4
  subject), and the no-effect knee-PF response is attributed to the paper's
  antagonist-only predecessor, the model class the archived project belongs
  to; the magnitudes disclaimer was kept. **F2 fixed by tightened wording** —
  the ports clause now reads "has received this current only at an extended
  stimulus window, not at the paper-exact \SIrange{2}{2.1}{\second} timing",
  which is exactly true of the Simulink twin run
  (`SNS_Simscape/dev/run_deng_perturbation_20261026.m`, 2.0--2.5 s window);
  no characterization of the twin's response is made in the patch. M1--M5
  fixed (caption siunitx + window wording; P02 anchor-by-neighbor + sentence
  case; P03 figure-environment fragment added). Open for Ben regardless:
  (a) the twin's "extensor suppression" vs phase-delay reinterpretation
  (REFERENCE_ONLY figure) — a twin re-run at the paper-exact window would
  settle it; (b) P02 volume/pages/DOI from Zotero before application.
  Supervisor round 2 verdict: SUPERVISOR_verdict.md (final GO/BLOCK there).

## Evidence

- Model / controller / effective parameters: the ORIGINAL archived project
  `K_A wirh RG march 21\K_A synergy.aproj` (ModelDB 264844; 90 non-spiking
  neurons, 12 muscles, 196 synapses), UNMODIFIED except: standalone export via
  AnimatLab's own serializer (`DengExport.exe` in-process route; the export is
  `deng2019_original_export.asim`, neuron/connexion/adapter counts and GUIDs
  verified against `deng2019_mine_20261025.json`), extended simulation window
  to 10.1 s, and added/extended voltage charts. All conditions are logged
  edits of that export (`run_*/edits.json`). Physics engine: Bullet, the
  project's own `AnimatModule`. No parameter of the network or body was
  changed in any run.
- Protocol / support / stimuli / perturbations: the published protocol of
  Deng et al. 2023 — external \qty{2}{\nano\ampere} tonic current on the
  rhythm-generator flexor population, and on the knee pattern-formation
  flexor population, each from \SIrange{2}{2.1}{\second} (paper-exact runs
  P1b/P2b); longer 2.0–2.5 s variants (P1/P2) and the project's own disabled
  −10 nA extensor-side block enabled as delivered (P3) also run for
  comparison. The project's as-delivered \qty{0.05}{\newton\meter} torque
  stimulus stays enabled in every condition, including the baseline.
- Results and definitions (burst = upward crossing of −60 mV; period from
  burst onsets; phase offset vs nearest baseline onset wrapped to ±half
  period): baseline free-runs at \qty{0.2316 \pm 0.0509}{\second} (left RG
  flexor, all-run) with legs in phase (L-vs-R flexor lag 0.007 s) and E/F
  alternating within each leg (lag 0.088–0.094 s); P1b extra flexor burst at
  t = 2.000 s (baseline next onset 2.433 s), extensor advanced to 2.124 s
  (baseline 2.589 s), post-window period 0.2238–0.2242 s vs baseline
  0.2212–0.2214 s; P2b post-window RG phase $-2.3 \pm 2.4$ percent and period
  0.2214 s vs 0.2214 s; P2 (0.5 s window) in-window left knee extensor
  motoneuron mean −86.7 mV vs baseline −64.5 mV. Numbers quoted in P01 are
  exactly these.
- Saved output paths (all in `Code/MuJoCo_SNS/spinal/reports_20261026_dengoriginal/`):
  `DENG_ORIGINAL_20261026_REPORT.md` (campaign report + verdict),
  `run_*/deng_rg_mn.txt` + `run_*/analysis.json` (per-condition stats),
  `phase_table_20261026.txt`, `summary_stats_20261026.json`,
  `compare_phase_20261026.py` (phase-offset definitions),
  `deng2019_original_export.asim`. Prereg:
  `D:\GitHub\_zcode_backups\diss_campaign_20261024\PREREG_B1_deng_original_20261026.md`
  (predictions registered before any run; two predictions falsified and
  reported there, including the original's 0.23 s period versus the Simulink
  rebuild's published-table 1.938 s period at `60-sns.tex:181` — the archived
  file and the published parameter tables are different parameterizations).
- Scope / author-confirmed facts: the ModelDB↔paper mapping (264844 = the
  Deng 2019 two-layer-CPG rat-hindlimb model's AnimatLab project) follows
  Ben's 2026-10-26 order and the project's internal structure (knee–ankle
  pattern-formation layers, 90 neurons); the 2023 paper's predecessor-vs-
  biarticular distinction is from its full text (quarantined until after the
  preregistered predictions were written). The paper's own 15.44%/15% and 5%
  phase numbers are NOT claimed for the archived project and are not quoted
  in P01.

## Validation

- Build or checks performed: NO compilation — this batch is proposals only;
  the reviewed chapters were not touched (per Ben's instruction), so no build
  of record exists for it. The figure generator ran clean over the saved
  chart files; the patch fragment is plain LaTeX (siunitx `\qty`,
  `\SIrange`, `\qty{2.3 \pm 2.4}{\percent}`) consistent with the chapter's
  existing macros; `\SIrange` needs no new package. diss-supervisor four-gate
  check over the fragments: see SUPERVISOR_verdict.md in this batch.
- Known reference/asset/layout/interpretation issues: (1) P01's anchors are
  from the easteregg2 snapshot; re-anchor against Overleaf if the paragraphs
  moved. (2) The figure is proposed, not placed; numbering/label are
  suggestions. (3) P02's bib entry is incomplete by design (volume/pages/DOI
  pending Ben's Zotero). (4) The twin reinterpretation is discussed as a
  question for ChatGPT/Ben (see REFERENCE_ONLY figure), not asserted in any
  patch. (5) The 0.23 s period is the archived project's own; it differs from
  the published-table rebuild's 1.94 s and P01 states the run identity
  explicitly so the two cannot be conflated.

## Preserved author decisions

- No tracked deletion, reply, or existing sentence is altered: both patches
  are insertions between existing paragraphs / bib entries.
- Figure 6.1(B) and all existing figure identities untouched.
- The "This suspended stepping result preceded…" deletion and all other
  README-listed rulings respected; nothing from a tracked deletion is
  reintroduced.
- The limits paragraph at `60-sns.tex:439` ("the systematic stimulation
  series … are not yet reproduced") was deliberately left untouched: it
  concerns the program's own networks, and P01's closing sentence carries the
  new boundary instead.
