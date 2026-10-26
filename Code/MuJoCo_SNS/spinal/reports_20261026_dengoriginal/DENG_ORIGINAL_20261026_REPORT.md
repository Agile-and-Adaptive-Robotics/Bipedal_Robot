# Order B1 — ORIGINAL Deng 2019 project: headless run + 2 nA perturbation (2026-10-26, egg)

**Machine:** easteregg2 (hostname checked). AnimatSimulator 2.x runs LOCALLY
here (the prompt assumed eb475ws4 — easteregg2 has its own install; the SSH to
eb475ws4 was used only to scp the original files over, read-only). No
simulation ran on eb475ws4.

## VERDICT (three-way)

**The original reproduces the paper's PERTURBATION FAMILY in the original
project file, and the paper's own account of this 2019-era model is confirmed
— but the effect is PHASE RESET, not burst suppression, which reinterprets our
Simulink twin's banked "suppression" result.**

1. **Original vs paper (Deng/Hunt/Chiel/Quinn 2023, "Biarticular muscles
   improve the stability of a neuromechanical model of the rat hindlimb"):**
   with the paper-exact protocol (2 nA, 2.0–2.1 s, RG flexor) the original
   shows exactly the paper's Figure-4 phenomenology — a **rapid short step at
   stimulus onset** (an extra flexor burst at t = 2.000 s vs baseline 2.433 s,
   followed by an early extensor burst at 2.124 s) and **rhythm continuation at
   the baseline period** (post-window period 0.224 s vs baseline 0.221 s, +1%).
   The paper's *quantitative* phase numbers (15.44% next-cycle delay, then 15%)
   belong to the 2023 BIARTICULAR model, which is NOT this file: in this 2019
   ModelDB model the ipsilateral side re-phases (offset scatters; slow 1%
   period drift accumulates) while the CONTRALATERAL side settles at a stable
   **+23.5% phase delay**. Same family (persistent phase shift with period
   recovery), not the same number — expected, since the paper says its model
   differs from this predecessor.
2. **Original vs the paper's account of the 2019 model:** confirmed. The paper
   states that in the earlier antagonist-only model "perturbing the pattern
   formation network of a single joint did not affect the neuron activations
   of the rhythm generator" — the paper-exact knee-PF run (P2b) shows RG phase
   −2.3 ± 2.4% (no effect). The 2023 model's 5% phase ADVANCE from a knee-PF
   kick is NOT present in this 2019 file, exactly as the paper's
   predecessor-vs-current contrast implies.
3. **Original vs our Simulink twin (banked):** the twin reported "extensor
   bursts fully suppressed during the 2.0–2.5 s window (0 in-window), then
   period 1.9377 s vs 1.938 reference". The original does NOT suppress
   extensor bursts — it phase-resets. **Leading reinterpretation (flagged for
   Ben, verifiable on eb475ws4):** the twin's 0 in-window extensor bursts is
   consistent with a 15–25% phase DELAY of its ~1.94 s cycle (0.29–0.49 s
   shift pushing the expected in-window burst past 2.5 s) — i.e., the twin and
   the original may agree in family (phase delay + recovery), and "suppression"
   may be a window-sampling artifact of the twin's slow cycle. The twin's
   period recovery (0.02%) and the original's (+1%) are the same phenomenon.

**Prereg scorecard (PREREG_B1_deng_original_20261026.md, written before any
run):** B0 "self-sustained no-drive rhythm" HELD; "period order 1–3 s"
FALSIFIED (actual 0.23 s — the original runs ~8x faster than the twin's
parameterization); "antiphase L/R" FALSIFIED (left and right fire in-phase,
lag 0.007 s — a synchronous two-hindlimb pattern; E/F ALTERNATION within each
leg held, lag 0.088–0.094 s ≈ 40% of cycle). P1 "extensor bursts interrupted
or suppressed in-window" FALSIFIED (they continue); "flexor phase perturbed /
phase-shifted with recovery" HELD; "recovery within 10% of baseline period"
HELD (+1%). P2 "strong knee-MN effect, small-or-zero RG shift, recovery" HELD
(paper-exact window: RG −2.3%, knee MN local transient; long window: knee
extensor MN mean −22.3 mV suppressed in-window, RG delayed +28%). P3 "brief
extensor-side event with recovery" HELD (phase jump in-window, post −2%).

## The pipeline (first-of-its-kind for this program)

The ModelDB 264844 tree contains NO .asim anywhere (project + Back-up + K&A
synergy 1 variants + the pristine zip all checked) — only GUI-dialect
.aproj. AnimatSimulator rejects aproj ("Element: \Simulation" missing). The
runnable standalone simulation was produced by the AnimatLab GUI's OWN
serializer, driven headless: DengExport.exe (compiled into the AnimatLab bin)
loads AnimatGUI.dll in-process, replicates AnimatLab2.exe's startup
(Server.Initialize + StartApplication on a message-pumped STA thread with a
built-in error-dialog dismiss-er), then LoadProject + ExportStandAloneSim.
Gotchas that cost the iterations (all documented in gui_export/): PowerShell
hosting fails (relative module loads resolve against the entry exe), x86
compile required (native 32-bit sim DLLs), an exe.config with the
DllExclusionList section is REQUIRED (CatalogPluginModules NREs without it),
work must be queued into StartApplication's modal loop after handle creation
+ settle, and ExportStandAloneSim writes the sim XML to the passed path
filename. Export fidelity verified: 90 Neurons / 196 Connexions / 42 Adapters
/ 56 SynAmps, stimulus GUIDs, and all key neuron GUIDs byte-checked against
deng2019_mine_20261025.json.

Conditions (build_conditions_20261026.py, every edit logged in each
run_*/edits.json): B0 baseline (stimuli exactly as delivered — the original's
own `Stimulus_1` torque nudge 0.05 N·m stays ON in every condition; the two
authors' TonicCurrent blocks stay disabled); P1 RG-flexor +2 nA 2.0–2.5 s
(the order's window); P2 knee-PF-flexor (LH_K&A PF F) +2 nA 2.0–2.5 s; P3 the
authors' own RG_Stin block enabled as delivered (−10 nA on LH_RG **E**,
2.0–2.1 s — note: the authors' apparatus targets the EXTENSOR side, not the
2023 paper's flexor); P1b/P2b paper-exact 2.0–2.1 s windows. All runs 10.1 s,
1 ms chart interval, Bullet engine (the project's own AnimatModule), each in
its own folder with console.log + chart files + analysis.json.

## Numbers of record

Baseline B0 (10.1 s): RG E and F each ~42–43 bursts; period LH 0.232±0.051 /
RH 0.226±0.011 s; L-vs-R flexor lag 0.007 s (in-phase); within-leg E-vs-F lag
0.088 (LH) / 0.094 (RH) s; knee extensor MNs burst ~44 cycles, knee flexor
MNs essentially silent (0–1 crossings); voltages saturate at 0 mV (IntegrateFire
clipping).

| run | stimulus | in-window effect | post-window period (LH_RG_E) | post phase vs B0 |
|---|---|---|---|---|
| B0 | none | — | 0.2214 s | — |
| P1b (paper-exact) | +2 nA LH_RG F, 2.0–2.1 | extra flexor burst at 2.000 (B0: 2.433); extensor advanced to 2.124 (B0: 2.589) | 0.2238 s (+1.1%) | ipsi re-phased; contra +23.5±10% delay |
| P1 (order's window) | +2 nA LH_RG F, 2.0–2.5 | flexor depolarized +2.5 mV, extra burst at 2.000; rhythm continues (no suppression) | 0.2224 s (+0.5%) | ipsi +25.3±4.2% delay; contra +13.6% |
| P2b (paper-exact) | +2 nA LH_K&A PF F, 2.0–2.1 | knee-PF local only (KAPF_F −6% phase); RG untouched | 0.2214 s (+0.0%) | −2.3±2.4% (no effect) |
| P2 (order's window) | +2 nA LH_K&A PF F, 2.0–2.5 | knee extensor MN suppressed −22.3 mV in-window; KAPF_F +2.2 mV | 0.2236 s (+1.0%) | ipsi +28%, contra −13% |
| P3 (authors' block) | −10 nA LH_RG E, 2.0–2.1 | phase jump (E dips to −68.6 mV), flexor freed at 2.002 | 0.2213 s (−0.05%) | −2.0±7% (recovered) |

Sign convention confirmed empirically: +CurrentOn depolarizes (P1 in-window
LH_RG_F mean −58.26 mV vs baseline −60.79 mV).

## Honest gaps / caveats

- The prompt's protocol window "2.0 to ~2.5 s" is wrong vs the paper: the
  paper says 2.0–2.1 s (both conditions). Both windows were run; the
  paper-exact runs (P1b/P2b) are the comparison of record.
- The 2023 paper's model is the BIARTICULAR extension; the ModelDB 264844
  file is the monoarticular predecessor. The paper's 15.44%/15% and 5%
  numbers are NOT expected to match this file — the paper's own text says the
  predecessor behaves differently (and that predecessor behavior is what we
  reproduce). Any claim of "reproducing Deng 2023's Figure 4/5" must be
  scoped to: original-project response matches the paper's description of its
  own perturbation protocol family.
- The poster (poster_p1.png) was not visually parsed in this session (image
  preview unavailable here); the comparison uses the full text only.
- The twin comparison is against Ben's banked numbers
  (results\deng_perturbation_twin_20261026.mat/.png on eb475ws4); the twin was
  NOT re-run here. The reinterpretation of its "suppression" as phase delay is
  flagged, not proven.
- First 1.5 s of every run are transient-free (charts start at 0; the rhythm
  is running from t≈0 — no kickoff needed, as prereg'd).
- Files DengExport.exe / DengExport.exe.config / DengExport.exe.lic now sit in
  the AnimatLab bin dir (outside the repo) — left in place for reuse; delete
  freely, they are regenerable from gui_export/DengExport.cs.

## Files

- `deng2019_original_export.asim` — the authors'-serializer export of record
- `build_conditions_20261026.py`, `analyze_conditions_20261026.py`,
  `compare_phase_20261026.py`, `make_figure_20261026.py`
- `run_B0_baseline/`, `run_P1_rgflexor/`, `run_P2_kneepfflexor/`,
  `run_P3_authorsblock/`, `run_P1b_rgflexor_2p1/`, `run_P2b_kneepfflexor_2p1/`
  — each: condition asim + edits.json + console.log + chart .txt
  (deng_rg_mn.txt = 11 voltage columns @ 1 ms) + analysis.json
- `deng_original_perturbation_20261026.png` — 4-panel evidence figure
- `originals/` — the scp'd pristine assets (aproj variants + mine json + dev
  scripts), untouched
- `gui_export/` — the in-process export toolchain (DengExport.cs + IL-dump
  scripts + the exported source copy); `probe_aproj/` — the failed
  direct-aproj probe log
- Prereg: `D:\GitHub\_zcode_backups\diss_campaign_20261024\PREREG_B1_deng_original_20261026.md`

## Reproduce (egg)

```
D:\Anaconda\envs\myo\python.exe build_conditions_20261026.py
cd run_B0_baseline && "/d/Program Files (x86)/NeuroRobotic Technologies/AnimatLab/bin/AnimatSimulator.exe" <abs path>\B0_baseline.asim > console.log 2>&1
D:\Anaconda\envs\myo\python.exe analyze_conditions_20261026.py run_B0_baseline --json run_B0_baseline\analysis.json
D:\Anaconda\envs\myo\python.exe compare_phase_20261026.py
```
(one folder at a time for the simulator; re-export from the original aproj:
compile gui_export\DengExport.cs (x86) into the AnimatLab bin, copy
AnimatLab2.exe.config → DengExport.exe.config, run `DengExport <aproj>`)

## Not run / not claimed

- No gait/forward-progression claims (the model air-walks in place in its own
  setup; this order concerns circuit perturbation response).
- No GUI open-check of the condition asims in AnimatLab2 (headless pipeline
  only); the exported asim loads and runs in AnimatSimulator cleanly, which
  is the acceptance that matters here.
- Follow-on (Li contact-feedback recipe) NOT started — budget went to the
  export pipeline; it is queued, not silently dropped.

— B1 executor, 2026-10-26
