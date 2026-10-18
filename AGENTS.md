# AGENTS.md — Standing context for AI coding sessions in Bipedal_Robot

Loaded automatically at session start. Keep it current; keep it lean.
**Cleaned 2026-10-18 at Ben's mandate** (he ruled the stale results, old mat-file
narratives, incorrect Xi values, HX711 detail, and old neuromechanics campaign logs
were confusing sessions). Everything cut is preserved verbatim in this file's git
history (diff against any pre-2026-10-18 revision) and in the referenced reports:
`Code\MuJoCo_SNS\spinal\DESIGN.md` (spinal stack history), `spinal\reports_*\`
(per-campaign reports), `Dig_out\` logs.

## Machines

**Machine-ID guard: check `hostname` FIRST and read only that machine's bullet.**
Paths in the wrong machine's bullet do not exist elsewhere — that is a different
machine, NOT staleness.

- **EB475WS4 (desktop, this repo at `D:\Github\Bipedal_Robot`)** — 12 logical cores,
  32 GB RAM (not a heavy-runs box; cap `parpool` ~6–8). Anaconda base
  `C:\ProgramData\anaconda3` (NOT writable). Spinal SNS env = **`myo`** at
  `C:\Users\Ben Bolen\.conda\envs\myo\python.exe` (py3.10.21, numpy 1.22.4 — scipy
  1.9.3 needs the 1.22 ABI, mujoco 2.3.7, sns-toolbox 1.5.2, optuna 5.0.0, torch cpu,
  matplotlib 3.7.5). Second env `MyoSuite` (py3.9, mujoco 3.1.2) at
  `D:\Users\Ben Bolen\.anaconda3\envs\MyoSuite`. **No git CLI on PATH** — use GitHub
  Desktop's bundled git: `C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-*\
  resources\app\git\cmd\git.exe` (glob the app-* version; works for log/diff/show).
  **Inline `python -c "multi-line"` gets EATEN by this shell — write script files.**
  MATLAB R2025a U1 at `D:\Program Files\MATLAB\R2025a` (on PATH; Optimization,
  Simulink, and gamultiobj/surrogateopt/patternsearch verified RUNNING headless;
  `license('test')` is unreliable here — trust a runtime test). **Simscape Multibody
  `smimport` is BLOCKED** (license server reserves all SimMechanics seats) — the
  URDF→smimport route runs on the laptop, not here. SOLIDWORKS 2025 SP3 at
  `D:\Program Files\SOLIDWORKS Corp` (SW2URDF add-in registered; SW COM fixed
  per-user after the C:-drive rebuild — real CLSID path in HKCU, pywin32 pip-installed
  in `myo`). MiKTeX user-scope at `D:\MiKTeX` (pdflatex
  `D:\MiKTeX\miktex\bin\x64`, on user PATH; all-users install staged at
  `D:\temp\miktex_setup\elevated_install2.ps1` if ever wanted). AnimatLab headless
  runs here: `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\
  AnimatSimulator.exe`. SCONE 2.4.4 at `C:\Program Files\SCONE` (Hyfydy trial
  licensed through 2026-10-24).
- **DESKTOP-5Q16KE9 (laptop, main workstation)** — repos under
  `C:\Users\Ben\Documents\GitHub\` (no D: drive). MATLAB R2025b, SOLIDWORKS 2025
  SP4.1. Conda `C:\Users\Ben\.anaconda3`; SNS env `myoconv`
  (`C:\Users\Ben\.anaconda3\envs\myoconv\python.exe`, py3.10.21, mujoco 2.3.7,
  sns-toolbox 1.5.2, myoconverter editable). **MuJoCo↔Simulink bridge BUILT + PROVEN
  here** (R2025b; rebuild driver `mujoco_bridge\matlab\laptop_build_bridge.m`).
  **Simscape Multibody license WORKS here for smimport AND hand-built models** —
  Ben's cylinder-elbow and CAD-import work belongs on this machine. Modest hardware
  (6 cores, 16 GB) — cap `parpool` at 6. Not SSH-set-up (as of 2026-09-28).
- **easteregg2 (heavy box, 10 cores, 128 GB RAM)** — repos under `D:\GitHub\`.
  **Heavy parallel optimization runs belong here.** Python = Anaconda at `D:\Anaconda`,
  NOT on PATH (bare `python` = Microsoft-Store stub): call
  `D:\Anaconda\envs\myo\python.exe` (py3.10.21, mujoco 2.3.7, sns-toolbox 1.5.2,
  opensim 4.4.1, optuna 5.0.0 — shared `optuna_walk.db` schema-matched with EB475WS4).
  **git CLI 2.54 on PATH** (unlike EB475WS4). SSH-reachable from the other lab
  machines (aliases `easteregg2` / `eb475ws4`; setup details machine-local only —
  this repo is PUBLIC). MATLAB R2025a. Ghostscript at `D:\Anaconda\envs\gs\Library\bin`.
- Custom skills are version-controlled in the `ZCode_Skills` repo (EB475WS4:
  `C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills`; easteregg2: same path — moved
  there 2026-10-17 from `D:\GitHub\ZCode_Skills` at Ben's request, junctions
  repointed, git history intact at 62b9ef7; laptop: `Documents\GitHub\ZCode_Skills`). On EB475WS4
  `C:\Users\Ben Bolen\.zcode\skills\` holds directory junctions into that repo
  (animatlab, latex-overleaf, matlab, mujoco, neuroscience, opensim, python,
  sns-toolbox, solidworks, scone, dissertation). easteregg2 uses the same junction
  pattern at `C:\Users\Ben Bolen\.agents\skills\`. Edit the repo copy; Ben commits.

## Project purpose (priority order)

**Ben's prose rules:** Oxford comma, always. Call the simulated bodies **models**
(human model, W2L robot model) — never "plants". Minimal em-dashes.

1. **Dissertation first.** Deadline: **submitted 2026-10-21; defended 2026-11-06**.
   LaTeX source in `Documentation\Reports and Papers\Dissertation\ProofFinal\`
   (+ sibling `Figures\<chapter>\` trees). Topic-chapter structure per the advisor's
   C12 directive: 10-introduction, 15-background, 20-muscles (Ben's Actuators 15(5):230
   paper), 30-torque (WIP manuscript bolen_isometric_2025), 50-database, 60-sns
   (tools-framed), 70-overall ("Discussion"), 80-conclusion, appendices 92/93/94; no
   ankle chapter. Chapters 3+ carry Methods/Results/Discussion sections.
   **In any dissertation session, FIRST read** the `dissertation` skill, then
   `Notes\revision_campaign_20261015\BEN_TODO.md`, `Notes\DISS_STATUS.md`, the newest
   `Documentation\Project_Digests\digest_*.md`, and the campaign folder
   (`style-contract.md`, `RULINGS_20261016.md`, `SNS_DISCUSSION_REDESIGN_20261017.md`,
   `FLOW_REVIEW_20261017.md`, `CAMPAIGN_REPORT.md`). Notation: Cxx/Txx/Rxx = advisor
   comment threads / tracked changes / resolved, all in `CHATGPT_REPORT.md` — gloss
   them for Ben. Standing rules: execute Ben's directives VERBATIM; no silent
   deletions or number changes (flag everything); every revision pass ends with a
   clean build (MiKTeX pdflatex full path, fresh `build_*` dir, bibtex with absolute
   BIBINPUTS/BSTINPUTS) + a per-directive compliance check with line evidence;
   NEVER round-trip thesis.bib through default-encoding PowerShell (Latin-1 byte-safe
   only — legacy Zotero mojibake); no machine-voiced comments in .tex; no dates or
   campaign provenance in prose; Modern-Robotics notation; engineering diameter
   symbol (`mathabx \diameter`), never `\phi`; "faster than real time" is an approved,
   literature-supported claim ONLY in its scoped form (what actually runs fast);
   minimal AI self-narration in repo docs. Campaign artifacts in
   `Notes\revision_campaign_20261015\`; defense-deck pipeline in
   `Dissertation\Defense\deck_builder\`.
   **2026-10-17/18 pass (this campaign):** journal-register writing pass over all
   chapters + FLOW_REVIEW TOP-10 fixes landed; open rulings applied (χ4 named in the
   Discussion, SADb naming aligned, "pool"→"fit set" in App. C + Ch. 4, stale
   commented numbers removed, easteregg2 wean results folded into chs. 6–7 — see
   the spinal section below for those numbers); failed torque figures regenerated
   (`*_r2.png`) + OpenSim-model isometric render staged in `Figures\30-torque\`
   (Ben dispositions wiring); Overleaf re-push staged. Verification build of record:
   ProofFinal `build_r7` (or the latest build_r* the campaign reports).
2. Design and control of bipedal humanoid robot legs with artificial muscles
   (PAMs/BPAs) controlled by a synthetic nervous system. Lab: AARL, PSU.
3. **Xi-correction-factor program** — the design factors χ0–χ3 (candidate χ4) that
   remove the rigid-body idealization from the BPA joint-torque model.
   **Values of record = the dissertation's Appendix C** (Tables `tab:app_xi_values`,
   `tab:app_opt_bounds`) and Chapter 4's adopted set — read those, NOT older numbers
   from git history or retired mats. Identification: `gamultiobj` Pareto fronts +
   holdout validation (the surrogate-vs-gamultiobj question is settled: gamultiobj for
   identification, surrogateopt+patternsearch for route redesign). One (χ1,χ2) pair
   shared across flexor and extensor configurations (advisor requirement, enforced).
   Open: a hip-frame second bracket on the flexor route (unit direction + stiffness
   projected on the force path), final 20 mm torque sweeps (Ben is machining/assembling
   the final parts), Routing27 v2 route-builder crash (12/27 actuators — lumbar +
   distal foot; post-defense).
4. **Mesh_Optimization rework** (`Opt_run.m` flexor / `Opt_run_Ext.m` extensor, with
   Morrow): redesign mechanics so each configuration (2× 20 mm BPAs per route) meets
   or exceeds the magnitude of select human monoarticular muscles over the knee RoM.
   Human torque/moment-arm targets from OpenSim (Gait2392; `OpenSim_Vasti_Results.txt`,
   `OpenSim_Bifem_Results.txt`). Results in `Mesh_Optimization\Results\`; the cam
   extensor design of record = `Vas_Pam_20mm_Result_cam_20261014_0416.mat` (min margin
   +30.83% vs required envelope, +39.77% vs plain target; pressure 269–507 kPa; bend
   radius 96.9 mm; clearance 43.1 mm; moment arm 59.4–102.8 mm). Plots must be
   journal-ready (Tol palette from `Code\Matlab\Colors.m`, Arial, 10 pt floor,
   CVD-safe markers AND line styles — never hue/dash alone, alt-text file per figure).
5. Run tests on the human leg; show predictions from 4 meet or exceed human torque
   values. First extensor 20 mm points: ExtTest_1_01 (1.163 / 7.849 / 23.437 N·m at
   10°/11°/10°, hung loads 3.7/25.2/76 N, 388.5/388/592 kPa); full sweeps pending
   hardware.

## Directory map

- `Code\Matlab\Mesh_Optimization\` — active optimization workspace (Routing27 v1/v2,
  ExtCam solver set: `Opt_run_ExtCam.m`, `objective_/nonlcon_/camToExtDesign_/
  buildDistalRingLocation_ *20mm.m`, `save_cam_result_20261014.m`).
- `Code\Matlab\Functions\`, `Code\Matlab\Robot_Data\` — shared helpers; MonoPam
  classes/data. `Code\Matlab\Human_Data\`, `Code\Matlab\Bone_Mesh_Plots\` — human/bone
  data and plotting. `Code\Matlab\minimizers\`, `Code\Matlab\Previous Optimization
  Code\` — legacy; don't build on.
- `Code\Matlab\HX711 v3.0\HX711 v3.0\` — the HX711/BPA load-cell app, and the ONLY
  home of it. **Live app = `HX711_BPA.m`**, a programmatic uifigure class; launch via
  `Start_HX711_BPA` from anywhere (the launcher addpaths its own folder). Calibration:
  original Tare→Scale-Factor flow or known-factor entry; a TARE UPDATES ONLY THE ZERO
  OFFSET, never the scale factor. PID pressure controller + step-test cal on the
  Pressure Ctrl tab. Factors persist in `hx711_bpa_last_cal.mat` next to the app
  (gitignored). Offline self-test: `test_HX711_BPA_offline` (hardware-free, run after
  pulling). Full detail, workflows, and the firmware notes live in
  `README_HX711_BPA.md` — read that, not this file. Companion DAQ firmware:
  `Testing_Data\2026_06_Festo\ValveDataAcquisition\ValveDataAcquisition.ino`
  (5-column output; sketch must match the app version).
- `Code\MuJoCo_SNS\` — MuJoCo + SNS-Toolbox pipeline (see its README.md): custom BPA
  muscle (`bpa_muscle.py`, port of the Festo force model), MuJoCo glue,
  `add_bpa_to_mjcf.py`. Converted models in
  `Solid_Models\OpenSim\Gait2392_Robotbody\mjc\`; `gait2392_robot.osim` there is
  GENERATED by `repair_robotbody_muscles.py` — don't hand-edit it.
  - **`spinal\` — the spinal-cord network stack. Full history + gotcha archive:
    `spinal\DESIGN.md` (read before touching runner.py / build_network.py /
    params.py). Per-campaign reports: `spinal\reports_2026{0923,0924,0925,1016}\`.**
    The stack: two-level RG+PF SNS networks (410-neuron production build: 92 MN
    pools + Ia/II/Ib afferents, per-leg RG half-centers, phase-shifted PF groups,
    stance-gated Ib load sharing, static-opt standing posture bias) coupled to the
    converted MuJoCo model at a **2 ms neural step = physics step**, forces from the
    BPA model inside the loop. Entry points: `runner.py` (stand→walk→stand; flags
    `--eval --fitted --bestN --view --scope --vest ...`; `HowToRunCode.md` for
    Spyder), `check_rhythm.py`, `audit_signs.py` (run DYNAMICALLY), `kine_ref.py`
    (acceptance metric), `optuna_walk.py` + `_curriculum*.py` family (dbs
    `optuna_walk.db`, `optuna_w2lvar.db`, `optuna_syn6.db`).
    **Contracts that MUST hold:** default build census (410, 376, 1186) + the
    bit-exact anchor eval **−161.56754173676563** (flavor: s3k stage-4 full_rules=0
    under current physics; the −160.23/−148.69 older baselines do NOT reproduce —
    stage renumber + left-ref fix moved them; do not conflate). Conditional topology:
    new neurons are built only when their gain > 0. JSON RULE: every params.G knob a
    study uses is saved in the winner json AND has an `if key in best` loader branch.
    Basin gate on finalists (±1% on every scalar; closed-loop form for tuned
    walkers). Anti-lottery: per-term diagnostics per trial + gain-sensitivity probe
    + wean ladder as the structural test. Frozen/NaN runs never count as walking.
    **BLAS/summation-order nondeterminism is real** and tuned nets sit near
    bifurcations — verify wean/rung winners in a FRESH process before believing them;
    in-process vs CLI contexts can land on different attractors with identical config.
    Production configs (censuses): `s3k` lineage (410/376/1186), `w2lvar`
    (888/382/7866), `syn6` (794/382/3033), plus spiking mirrors
    (`build_network_spiking.py`, env-selector, MN layer stays graded) and the
    AARL_NET selector for variants. **v10-era `--fitted --bestN` winners no longer
    walk under current physics — only s3*-lineage winners are current.**
    **Results of record (2026-10-16/17 wean campaign, reports_20261016\):** HUMAN
    model: feedback re-tuned at drive scale 0.381 → kine **−145.07** (10+10 cycles,
    duty 0.71/0.65, basin 8/8; `wean_human\round8_best3.json`); wean ladder holds
    bilateral stepping at EVERY drive scale down to **0.0** (8+11 cycles, −157.07 —
    feedback-sustained, un-retuned config dies below 0.85); generalizes across Ong
    references (−159.66 @ 0.75 m/s …); still in-place under the support rig, no
    forward progression; two lottery weights reported (toe_rge, phase_reset_f —
    removal does not degrade). ROBOT (W2L): L/R split drive + tuned feedback →
    antiphase progression 0.0662 m/s, +1.718 m/30 s, 25/23 cycles, drive weaned to
    0.60, basin 12/12 (`w2l_wean\w2l_wean_best_params.json`); below 0.60 dies;
    harness NOT weaned; duty below the strict walking gate. **No harness-free walking
    anywhere in the program's records — state it exactly that way.**
    `spinal\gait_refs\*.npz` — the 43-cycle locomotion reference library (Arnold
    running 8 subj × 4 speeds; Ong walking 0.5–2.0 m/s + self-selected, drop
    `ong_selfsel_Init200`'s degenerate left cycle if sampled; Falisse predicted
    walking) — this is "the many-subject IK references" (staging history:
    `D:\temp\gait_lib_staging`, source thumb drive `F:\Biomechanics data`).
    Knee convention: **flexion NEGATIVE** (converter preserved OpenSim's), joint
    ships `limited="false"` pre-patch. `spinal_run.npz`/`spinal_run.png` get
    overwritten by evals — redirect via AARL_NPZ or git-restore.
    Tools around it: `connectome_block_editor.html` (v3.3, self-contained —
    templates EMBEDDED; Ben's exports land in
    `Neuromechanical_Models\Mujoco_SNS_models\` — look there first; localStorage is
    per-origin: clean TEST_* artifacts after agent testing), `neuro_scope.py`,
    `draw_circuit.py`, `basin_gate.py`/`basin_closed_*.py`, `synergy_model.py` +
    `synergy_basis.npz`. `spinal\replication\` — Shevtsova/Shinohara/Rybak
    circuit-rule translations (per-edge citations; NOT built/verified as running
    models — the dissertation says so honestly). `spinal\w2l_cpg\` (Ben's
    contact-driven bilateral CPG in SNS, 95n/208syn — 91 is stale; W2L_SMOKE gate)
    and `spinal\w2l_mujoco\` (W2L .aproj→MJCF, `w2l_mjcf_fixed.xml`; canonical CLI
    knobs amp=1.5 pulse=0.5 ibnA=1.0; the AnimatLab originals run under a VIRTUAL
    WALKER (WalkingPath body) — harness march is the walker analog, free standing is
    structurally impossible for that geometry). L/R drive knobs (2026-10-16):
    alpha_l/gamma_l/alpha_r/gamma_r, default-0, default-off = bit-identical
    (`drive_defaults_gate.py`). MuJoCo fidelity ruling: keep 2.3.7/implicitfast/2 ms
    through the dissertation; Hyfydy (m2012fast elastic tendon, Hunt-Crossley
    contact, PSEM-TW integrator, ~95× RT) is the post-deadline engine option.
- `Code\Arduino\`, `Code\Festo\` — embedded/valve hardware code.
- `Code\Hardware\walker_io\` — Python walker IO stack (Orin daemon, lab console,
  valves, sensors, afferents; README inside). Belongs with the hardware code after
  the reorg below.
- `tmp\` — session scratch at repo root; disposable by convention.
- **POST-DEFENSE REORG PLAN (Ben, 2026-09-30 — execute only after submission; run
  the repo-hygiene refs-impact protocol first):** consolidate SCONE, AnimatLab,
  SNS_Simscape, and SNS-toolbox model content under `Neuromechanical_Models\`;
  organize `spinal\`'s loose files internally; move `walker_io` next to the other
  hardware code.
- **Xi1/Xi2 semantics (Ben, 2026-09-07):** *effective system-stiffness parameters*,
  not literal bracket beam stiffness — the fitted compliance lumps bracket,
  fixtures, and the test-stand cable winch. Don't over-interpret structurally.
  **CV-mat fold bookkeeping (Ben, 2026-10-18 — read before quoting any
  "validation" test):** in the pinned-knee front mats' `filtered_results`, the
  leading columns are the fold's held-out (validation) test indices.
  Flexor `minimizeFlxPin10_results.mat`: columns 2–3 = the two validation tests;
  "noT3" means the 47 cm test was NOT evaluated in those runs at all — the folds
  live on the other four tests (2 held, 2 trained). Adopted flexor pick row 77:
  validation {46 cm, 40 cm-tendon}, training {48 cm, 41 cm}. Extensor
  `minimizeExt10mmX3_results.mat` (`holdCols=[2 3 4]`): columns 2–4 = the three
  held-out tests of that member's fold over the six used tests; the 42 cm and
  48 cm tests were never used (bracket shift on those early tests). Adopted
  extensor pick row 32:
  validation {40, 40-tendon, 47 cm}, training {43, 43-tendon, 46 cm} — 3
  training / 3 validation / 3 not-used (42, 42-tendon, 48 cm) over the nine
  (Ben, 2026-10-18). Flexor: 2 training / 2 validation / 1 not-used (47 cm).
  The mats' own
  `PICK` fields are legacy — the picks of record are what
  `buildKneeFlexorContext20mm.m` (77) and `buildKneeExtContext20mm.m` (32) load.
  Bracket reference points are Ben's CAD interpretation choices (current: flexor
  insertion Pbri, flexor origin Pbr2 (pinned-flexor only), extensor Pbr — exact
  values in Appendix C of the dissertation and the evaluator headers). Frame
  convention (1- vs 2-rotation) is provably irrelevant for y-symmetric K=[a,b,a]
  arrays; the per-convention K-orderings are Ben's structural preference, carried
  not derived. The pinned-flexor data cannot uniquely determine the stiffness pair
  (flat likelihood valley) — the adopted pair was selected by the external
  requirement that both redesign configurations meet the 5% overshoot target.
- `Testing_Data\` — **important.** `2022_02_Festo\` = the Xi-minimizer family.
  **Evaluator of record: `minimizeFlxPin.m`** (full two-bracket mechanics; MIXED
  convention since 2026-09-13: insertion/tibia bracket one-rotation pitch-only frame
  K=[X1,X2,X1], origin/hip bracket two-rotation K2=[X1,X1,X2]; +5.3° encoder offset
  on kf(3); env `FLX_T5_YMM` for the test-5 insertion-bracket y-shift).
  **Pinned-knee records of record (Ben, 2026-10-18): `minimizeFlxPin10_results.mat`
  (flexor front, row 77 pick) and `minimizeExt10mmX3_results.mat` (extensor front,
  row 32 pick, pair-locked; verified 2026-10-18 to carry the same fronts, picks,
  LOCKSRC, and baselineScores as the dated duplicates). Ben will DELETE the other
  pinned-knee result mats to save space — never build on or resurrect them.**
  `minimizeFlxPin2brk.m` and the other legacy evaluators are STALE — do not re-run
  as-is. Harnesses follow Ben's Collect/Dig scheme (`Dig_*` analyze existing mats,
  `Collect_*` run sweeps; logs/outputs in `2022_02_Festo\Dig_out\`). Extensor
  identification used tests {1,2,5,6,7,8} of nine: adopted fold fitting {5,6,7},
  validation {1,2,8}; tests 3/4/9 (42, 42-tendon, 48 cm) not used — the attachment
  bracket shifted during those early tests. Subfolders Flx/Ext × sizes are
  needed only by certain plotters — add to path transiently and REMOVE afterward.
  `2026_06_Festo\` — long raw files summarized in
  `Results_table_10mm_pinned.xlsx` and `Results_table_20mm.xlsx` (read the tabs, not
  the raw files). 20 mm data-entry tabs `ExtTest20mm_1` + `FlxTest20mm_51cm`;
  measured-torque and pressure back-calculation sections live in
  `Knee_Flexor_Data_20mm.m` / `Knee_Extensor_20mm.m` — UPDATE the d/ang
  reaction-point constants after rig measurement before trusting measured torque.
  **LC-angle convention (Ben, 2026-10-04, authoritative for ALL new tests): the
  load-cell angle is measured FROM THE TIBIA AXIS — torque about t1 =
  F·sin(LC+0.83°)·d; all new consumers convert with 90−LC before the adjoint
  transform. The 2022-era arrays used the OLD convention — do not mix.** Also holds
  the BPA force/pressure characterization structures (some from Muscle_Sensory).
- `Solid_Models\Biomimetics_2022-Knee_Test\` — knee test setup CAD/STLs/point
  clouds. `Solid_Models\OpenSim\Gait2392_Robotbody\` — human reference models.
  Provenance: `ConnorBipedal.osim` (+_Adjusted variants) = Ben's untouched
  master's-era robot (Connor Morrow, 2020); `gait2392_robotbody.osim` = the
  2026-evolved copy (LEFT-leg reference model; right leg carries the custom paths).
  `ConnorBipedal_symmetric.osim` and `gait2327.osim` are GENERATED (mirror /
  XML-surgery scripts) — regenerate, don't hand-edit. MyoConverter outputs are
  regenerable (~116 MB, gitignored) — rebuild with MyoConverter, don't restore from
  git. Bone geometry for packaging: `Code\Matlab\Bone_Mesh_Plots\Open_Sim_Bone_
  Geometry\` point clouds + MyoConverter STLs. Master's data thumb drive =
  `D:\Bipedal humanoid` (on machines that have it).
- `Documentation\Reports and Papers\` — papers; `Dissertation\` = priority 1.
  `Documentation\Research Notes\` — lit reviews with V/L/U provenance tags.
- `Neuromechanical_Models\` — AnimatLab CPG walker models + Ben's connectome
  exports (his downloads land in `Mujoco_SNS_models\`). The instrumented phase-1
  file is `Biped_2xCPG_wSubs\walk new new tester added 2 axis_phase1.asim`; a fresh
  standalone export exists (2026-09-09). The .aproj is mid-surgery (bilateral RH
  wiring, biarticulars added; reorganized model latches under tonic drive —
  conductance tuning is the active step; no locomotor claim beyond verified
  assembly). Canonical connectome: `ben_connectome_canonical.json` (110n/216e).
  Wiring rulings of record (incl. the CORRECTED Deng audit: afferents DO reach
  the RG/PF half-centers in the Deng original — retraction of the earlier "zero
  afferents" claim, Ben-verified 2026-10-03): `spinal\WIRING_RULINGS_20260903.md`.
  Latch lesson: Deng-style persistent-Na half-centers need tau_h FIXED (voltage-
  dependent tau_h quenches burst termination; `spinal\deng_cpg_ode.py`,
  `_tau_h_check.py`). AnimatLab-vs-SNS_Library synapse saturation conventions are
  OPPOSITE — keep straight when porting.
  `SNS_Simscape\` — the Simulink SNS block library (12+3 blocks; README inside;
  MuJoCo bridge proven on the laptop).

## Standing conventions

- Execute Ben's directives verbatim; ambiguities get asked, not inferred. Every
  revision pass ends with verification + a per-directive compliance check.
- No silent deletions or number changes anywhere in the dissertation; every number
  in prose traceable to a Testing_Data file, a generating script/mat, or a citation.
- Ben commits via GitHub Desktop; agents leave the working tree clean for him and
  say what to commit. NO commits/pushes by agents unless Ben says so.
- This repo is PUBLIC — never put credentials, license details, or SSH setup here.
