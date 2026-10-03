# Goal 3 — AnimatLab headless baselines + spiking-neuron copies (2026-10-02)

Machine: easteregg2. AnimatLab 2 at `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab`.
All runs headless via `AnimatSimulator.exe`. Every model file Ben owns is UNTOUCHED — all work
happened on copies inside this report directory (`git status` shows only untracked additions in
`Code/MuJoCo_SNS/spinal/reports_spiking_20261002/`; the tracked-file modifications visible in
`git status` belong to the sibling SNS-Simscape/MuJoCo spiking campaign sharing this directory,
not to this AnimatLab task).

Runner used for every simulation (quote the path — the model folders contain spaces):

```
"D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe" "<asim>"
```

packaged by `tools\run_all.cmd` which writes `logs\<tag>_<name>.log` per run.

## 1. What "spiking" means in these files (verified, not assumed)

- The AnimatLab public source (`IntegrateFireSim/Neuron.cpp`, fetched 2026-10-02 from
  animatlab.com/SDK) confirms the `.asim` has ONE serialized `<Neuron>` class and **no spiking
  flag**: `m_bSpike=(m_dMemPot>m_dThresh)`. Spiking vs non-spiking is a parameter regime —
  the non-spiking family ships `InitialThresh` at +50/+200 mV (unreachable); the Li/2023
  spiking family ships −55 mV (IN/PF/RG/contact), −59 mV (MN), and keeps 12 graded "MV" output
  cells at rest −100 mV.
- Synapse kinds live in THREE places that must agree: the `<SynapseType><Type>` string, the
  section the type sits in (`<SpikingSynapses>` / `<NonSpikingSynapses>` / `<ElectricalSynapses>`),
  and **each `<Connexion>`'s own `<Type>` enum (0 = SpikingChemical, 1 = NonSpikingChemical,
  2 = Electrical)** — established empirically here from the Li/phase1 connexions (all 0) vs the
  W2L-family connexions (all 1). Missing the third place kills the run (§6, trap T1).

## 2. Model inventory (from the .asim/.aproj XML, read-only)

| model | neural style | neurons | synapse types (S/NS/E) | connexions | baseline asim run |
|---|---|---|---|---|---|
| `Walker_2_Layer_CPG_Standalone_modern.asim` | all non-spiking (thr 50/200) | 80 | 8/17/2 | 150 | yes, 5.1 s |
| `Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim` | all non-spiking | 96 | 8/20/2 | 208 | yes, 5.1 s |
| `Biped_2xCPG_wSubs_Standalone.asim` (STALE export) | all non-spiking | 122 | 8/50/2 | 243 | yes, 5.1 s |
| `walk new new tester added 2 axis_Standalone.asim` ("Li 2023") | spiking (32×−55, 12×−59, 12 MV) | 56 | 8/3/2 | 84 | yes, 10.02 s |
| `walk new new tester added 2 axis_phase1.asim` (instrumented) | spiking (same 56-neuron set as Li 2023 — ID sets identical, verified) | 56 | 8/3/2 | 84 | yes, 10.02 s |

Findings on lineage (checked by md5 + XML):
- The `Li Model` folder's two asims are **byte-identical** to the `Walker_2_Layer_CPG` folder's
  copies (md5 `ce756e…`, `5e167f…`); the folder's unique artifact is `walk tester rearranged.aproj`
  (44 Spiking + 12 NonSpiking neurons, its `walk tester rearranged.asim` does not exist in the
  repo). The "Li model asim" baseline is therefore the shared 2023 spiking walker, run from
  `models\Li_2023\Li_2023.asim`.
- `Biped_2xCPG_wSubs.aproj` today holds 122 NonSpiking neurons, 0 Spiking (GUI classes), and its
  `Biped_2xCPG_wSubs_Standalone.asim` predates the 9/26 synapse restores (the ask's STALE trap).
  No fresh asim export exists (GUI export only), so the Biped family is bracketed by BOTH the
  stale non-spiking export AND the instrumented phase-1 spiking asim — both were run as baselines.
- The phase-1 asim is the 2023 spiking-lineage instrumented export (DataTool_8 charts both legs'
  joint angles in degrees plus RG/PF/MN voltages), NOT an export of the current aproj.

## 3. Baselines (all exit 0, "Simulation stopped" at the configured end)

| metric | W2L modern (graded) | BilateralRG (graded) | Biped stale export (graded) | Li 2023 (spiking) | Biped phase1 (spiking) |
|---|---|---|---|---|---|
| RG burst period | 0.443 s (10 bursts/5 s) | 0.462 s (10–11/5 s) | none — latched (0–1 onsets) | n/a (not charted) | 1.305 s (8 bursts/10 s) |
| RG burst freq | 2.26 Hz | 2.16 Hz | dead | n/a | 0.77 Hz |
| L/R antiphase | n/a (single L RG charted) | folded lag −0.42 cycle, r=0.869 | n/a | n/a | folded lag +0.494 cycle, r=0.703 (RG_L vs RG_R stance) |
| hip range | −15.0..23.2° | −6.9..23.3° | −15.1..18.4° | n/a | L −15.1..23.1°, R −15.1..23.2° |
| knee range | −1.0..60.3° | −0.3..37.2° | n/a (not charted) | n/a | L −0.6..55.5°, R −0.4..55.9° |
| ankle range | −20.4..−4.0° | −20.3..−4.1° | n/a | n/a | −20.9..0° both |
| contacts / travel | n/a | n/a | n/a | L foot active 55.5%, 8 episodes; R foot 44.3%, 7; height 0.952..1.020 m; distance −3.46 → +2.91 m (walks) | same chart set |

Cross-checks vs record: W2L 2.26 Hz and BipRG 0.462 s / near-half-cycle lag reproduce the
2026-09-24 session numbers (2.22 Hz joints, 0.464 s, +0.474 cycle); phase1's 1.305 s reproduces
the 2026-09-25 `Li Model\DataTool_7.txt` period of record. Evidence files:
`metrics_baseline_*.json` (7 files) + the chart byproducts under `models\<model>\`.

## 4. Spiking copies — what was built

Converter `tools\make_spiking.py` (pure string surgery on a copy; originals untouched):

- every neuron with `InitialThresh ≥ 0`: threshold → spiking regime, `RelativeAccom` → 0
  (Li-native). v1: −55 mV (−59 for MN-named cells). v2 calibrated: per-class thresholds from the
  measured baseline graded peaks (`baseline_graded_ranges.json`): RG/PF/IN −53.5 (peaks ≈ −52),
  MN −44 (peaks ≈ −42.4), RE −44.5 (peaks ≈ −43).
- every `NonSpikingChemical` type: re-emitted as `SpikingChemical` (same ID/Name/Equil/SynAmp;
  SynAmp now = per-spike conductance increment), moved into `<SpikingSynapses>`;
  Decay 10 ms (v3 probe: 30 ms), MaxRelCond = max(5, 2·SynAmp), Hebbian/facilitation off.
- every connexion referencing a converted type: `<Type>1</Type>` → `<Type>0</Type>`.
- IDs unchanged everywhere, so connexion wiring, stimuli, and chart-column GUIDs stay valid;
  no cloning, so no fresh GUIDs were needed.

| copy | neurons converted (MN) | types | connexions flipped | validate | run |
|---|---|---|---|---|---|
| `models\W2L_modern_spiking\W2L_modern_spiking.asim` (v1) | 80 (24) | 17 | 150 | PASS | 5.1 s clean |
| `models\BipRG_modern_spiking\BipRG_modern_spiking.asim` (v1) | 96 (24) | 20 | 208 | PASS | 5.1 s clean |
| `models\Biped_standalone_stale_spiking\Biped_standalone_stale_spiking.asim` (v1) | 122 (22) | 50 | 243 | PASS | 5.1 s clean |
| `…_spiking_v2\…` (calibrated thresholds, ×3 models) | same | same | same | PASS | 5.1 s clean each |
| `models\W2L_spiking_v3_decay30\…` (Decay 30 ms probe) | 80 | 17 | 150 | PASS | 5.1 s clean |
| `models\bisect\W2L_neuronsonly.asim` (ablation: thresholds only, graded synapses kept) | 80 | 0 | 0 | PASS | 5.1 s clean |

`tools\validate.py` per copy: whole-file XML parse, per-page CDATA parse (vacuous — the .asim
carries no `DiagramXml`; drawings are .aproj-only), neuron-ID / connexion-tuple / synapse-type-ID /
chart-column / stimulus sets identical to the original, no thresholds ≥ 0 remain, no
NonSpikingChemical types remain, every connexion `<Type>` enum matches its synapse class, every
SynapseTypeID resolves, SimEndTime (5.1) exceeds every chart EndTime (5). The one WARN
(stretch-receptor chart targets "unknown") is a checker limitation — those are `<RigidBody>`
elements, identical in the originals, and their charts write data.

## 5. Spiking-copy results vs baselines (honest outcome: the rhythm does not survive)

| variant | RG spiking | rhythm verdict | joint ranges (L) |
|---|---|---|---|
| W2L v1 (thr −55/−59, decay 10) | 2 / 0 spikes in 5 s | dead — network settles to fixed point | hip 6.1..6.9°, knee 10.9..12.2° (baseline: −15..23°, −1..60°) |
| W2L v2 (calibrated thr) | 2 / 0 | dead | hip 0.5..0.8° |
| W2L v3 (decay 30) | 2 / 0 | dead | hip −3.8..−3.8° |
| BipRG v1 | L 6–9 Hz; R 442–448 Hz | asymmetric latch — R half-centers saturate at max rate, L nearly silent | hip pinned 23.0°, knee ~0°, ankle −20° |
| BipRG v2 | flx cells 450 Hz, ext silent | latch | hip 23.0°, knee 32.1..32.6° |
| Biped v1 | 58 bursts on BOTH HCs @ 0.085 s (11.8 Hz) | **synchronous co-bursting, no alternation** — L RG ext-vs-flx spike-count envelopes are IDENTICAL (58 spikes each; xcorr r = 1.000 at lag 0.00 s; convention and number persisted in `metrics_Biped_standalone_stale_spiking.json` → `same_leg_xcorr`, tool `tools\sameleg_xcorr.py`) — a locked 11.8 Hz tremor replacing the latched baseline | hip −15.3..18.5° |
| Biped v2 | 0 spikes | dead | hip −15.1..18.6° (passive settle) |
| ablation (thresholds only, graded synapses) | 1 / 0 spikes | **rhythm preserved**: hip −15.1..23.2°, knee −0.8..47.9°, ankle −20.2..−4.5° vs baseline −15..23°, −1..60°, −20..−4° | — |

Mechanistic read (from the measured traces): the graded half-centers sit at −61..−52 mV with
continuous mutual inhibition (RG-Inhibit SynAmp 2.749 µS always on during a plateau). Converting
to per-spike conductances at the same amplitude either starves the loop (W2L/Biped-v2: a plateau
crossing the threshold fires only 0–2 spikes, far too little inhibitory charge to switch the
other half-center) or saturates it (BipRG: threshold 3–7 mV below the graded peak → 450 Hz max-rate
firing). The ablation isolates the cause: with spiking-regime thresholds but graded synapses the
network still oscillates through subthreshold voltages (graded synapses read Vm continuously), so
**the synapse conversion alone breaks the rhythm**; a spiking version needs its synaptic
amplitudes, decays, and AHP re-tuned as a rate code (Li's model was tuned that way from birth —
its RG/PF cells use GMaxCa 5–6 with tonic 5–6 nA and stock IPSP/EPSP types). This is a finding,
not a hidden failure; the pipeline, validators, and copies are in place for that tuning campaign.

## 6. Traps hit (new ones for the record)

- **T1 (NEW): `<Connexion><Type>` is a second class tag.** Converting only the SynapseType makes
  AnimatSimulator die at load with only `Critical error occurred: A critical simulation error has
  occurred.` (no detail; `Simulation stopped. Time: 0`). The enum must be flipped to 0 for
  SpikingChemical. Evidence: `logs\spiking_*.log` (crashed, pre-fix) vs `logs\spiking2_*.log`
  (clean, post-fix). `tools\validate.py` now checks enum/class consistency.
- **T2: shared chart filenames collide.** Chart `OutputFilename`s are relative and identical
  across the W2L-family models ("Rhythm Generator.txt", …); running several models from one
  directory overwrites each other's byproducts. Fix: one directory per model under `models\`.
- **T3: trailing zero-fill rows** in chart .txt (the skill's documented trap) re-confirmed —
  9 all-zero rows at t ≥ 5.0002 inflate maxima to exactly 0.0 mV and add a fake crossing;
  `tools\analyze.py` strips them.
- **T4: a run that crashes at load TRUNCATES the chart .txt files** it opened — the crashed
  synapses-only bisect zeroed the ablation's byproducts; the ablation had to be re-run
  (`logs\bisect2_*.log`).
- **T5: bash `\$m` inside Windows-path strings** collapses to literal `$m` — use forward-slash
  paths when scripting the SNS python from git-bash.
- Honored from the standing list: no original model touched (copies only), chart GUIDs left
  intact by keeping all IDs, SimEndTime > chart EndTime preserved, no GUI ever opened so nothing
  was saved from a GUI, and no .aproj was handed to Ben (none created — see §7).

## 7. GUI-verification status + scope notes

- All artifacts created here are `.asim` files (runnable headless — that is how they were
  verified) plus scripts, logs, and metrics. **No `.aproj` spiking mirrors were built**, so no
  `.aproj` needs the GUI-open rule; if Ben wants GUI-editable spiking versions, the same
  conversion must be re-implemented against the .aproj schema (Value/Scale/Actual triplets,
  `SynapticTypeID`, page CDATA drawings) and then GUI-opened per the standing rule.
- **PENDING BEN GUI CHECK** — none of the spiking-copy .asim artifacts below has been opened in
  AnimatLab2.exe (no GUI in this headless session); they are verified by headless running only:
  `models\W2L_modern_spiking\W2L_modern_spiking.asim`,
  `models\BipRG_modern_spiking\BipRG_modern_spiking.asim`,
  `models\Biped_standalone_stale_spiking\Biped_standalone_stale_spiking.asim`, the three
  `…_spiking_v2\…` copies, `models\W2L_spiking_v3_decay30\…`, and the two
  `models\bisect\W2L_*.asim` files. If Ben accepts headless-run-verified .asims without the
  label, record that ruling here and drop the marker.
- The Biped family's spiking copy is built on the STALE `Biped_2xCPG_wSubs_Standalone.asim`
  export by necessity (no fresh export exists; headless sessions cannot export). Re-exporting the
  current aproj and re-running `make_spiking.py` is a 2-minute refresh once Ben exports.
- The sibling campaign sharing this directory owns `goal2_*.m`, `goal2_baselines.mat`,
  `logs\baseline_{s3k,syn6,w2lvar}.log`, AND the following `tools\` scripts (their headers
  reference SPIKING_MIRROR_PLAN.md / build_network_spiking.py gates):
  `baseline_eval.py`, `basin_gate_spiking.py`, `check_selfsustain_spiking.py`,
  `probe_build_spiking.py`, `probe_hybrid.py`, `probe_lif.py`, `smoke_air_spiking.py`,
  `topology_mirror_check.py`, `tune_rg_pair.py`, plus `tools\__pycache__`. Every other file
  here is this task's.

## 8. Files created by this task (all under `Code/MuJoCo_SNS/spinal/reports_spiking_20261002/`)

- `tools\`: `inventory.py`, `analyze.py`, `make_spiking.py`, `validate.py`, `bisect_variant.py`,
  `sameleg_xcorr.py`, `run_all.cmd`
- `models\<model>\` ×5 baselines (asim copy + chart .txt byproducts): `W2L_modern`,
  `BipRG_modern`, `Biped_standalone_stale`, `Li_2023`, `Biped_phase1`
- `models\` spiking copies (asim + byproducts): `W2L_modern_spiking`, `BipRG_modern_spiking`,
  `Biped_standalone_stale_spiking`, `…_v2` ×3, `W2L_spiking_v3_decay30`, `bisect\`
  (`W2L_neuronsonly.asim` + byproducts, `W2L_synonly.asim` kept as the T1 crash evidence)
- metrics: `metrics_baseline_*.json` ×5, `metrics_*spiking*.json` ×7,
  `metrics_bisect_W2L_neuronsonly.json`, `baseline_graded_ranges.json`
- logs: `baseline_*.log` ×5, `spiking_*.log` ×3 (crashed, pre-T1-fix), `spiking2_*.log` ×3,
  `spikingv2_*.log` ×3, `v3_*.log`, `v3b_*.log`, `bisect_*.log` ×2
- this report. Total ≈229 MB (mostly chart byproducts; the keepers are the asims, JSONs, logs,
  and this report if Ben wants the directory slim before committing).

Pre-existing files modified: **none** (verified `git status --porcelain`; only untracked
additions under this directory from this task).

## 9. Reproduce

```
:: baselines
"D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatSimulator.exe" ^
  "D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002\models\W2L_modern\W2L_modern.asim"
:: (same for BipRG_modern, Biped_standalone_stale, Li_2023, Biped_phase1 in their folders)

:: build + validate + run a spiking copy
D:\Anaconda\envs\myo\python.exe tools\make_spiking.py models\W2L_modern\W2L_modern.asim models\W2L_modern_spiking\W2L_modern_spiking.asim
D:\Anaconda\envs\myo\python.exe tools\validate.py  models\W2L_modern\W2L_modern.asim models\W2L_modern_spiking\W2L_modern_spiking.asim
D:\Anaconda\envs\myo\python.exe tools\analyze.py   models\W2L_modern_spiking\W2L_modern_spiking.asim --mode spiking --json metrics_W2L_modern_spiking.json
```
(run from `Code\MuJoCo_SNS\spinal\reports_spiking_20261002`; v2 adds `--calibrated`, v3 adds
`--decay 30`; baselines analyze with `--mode graded`, Li/phase1 with `--mode spiking`.)
