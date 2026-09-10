# AGENTS.md — Standing context for AI coding sessions in Bipedal_Robot

Loaded automatically at session start. Keep it current; keep it lean.

## Machines

- **DESKTOP-5Q16KE9 (laptop), the main workstation** — newest/premium
  MATLAB + SolidWorks, but modest hardware (6 cores, 16 GB RAM). MATLAB **R2025b** at
  `C:\Program Files\MATLAB\R2025b`; SOLIDWORKS **2025 SP4.1** (33.4.1) at
  `C:\Program Files\SOLIDWORKS Corp`. Its GitHub repos live under
  `C:\Users\Ben\Documents\GitHub\` — **no D: drive on this machine**.
- **easteregg2 (desktop; sessions there read this file too)** — 10 cores, 128 GB RAM; older
  installs (MATLAB R2025a, SOLIDWORKS 2025 SP03). Repos live under `D:\GitHub\`.
  **Heavy parallel optimization runs belong on easteregg2**;
  `parpool(10)` only makes sense there — on the laptop cap the pool at 6.
  **Python = Anaconda at `D:\Anaconda`, NOT on PATH** (why PyCharm says "no python"):
  call env pythons by full path — MuJoCo/SNS env = `D:\Anaconda\envs\myo\python.exe`
  (py3.10, mujoco 2.3.7, sns-toolbox 1.5.2, opensim 4.4.1; this machine's name is `myo` —
  `myoconv` is the laptop's copy of the same toolchain). Other envs: `opensim`
  (py3.11 + opensim 4.6), `d2l`. PyCharm: Add Interpreter → Conda → executable
  `D:\Anaconda\condabin\conda.bat` → existing env `myo` (never let PyCharm create a venv).
- Custom skills (`matlab`, `solidworks`, `latex-overleaf`, `myoconverter`) are version-controlled in
  `Documents\GitHub\ZCode_Skills`. On this machine `C:\Users\Ben\.zcode\skills\` holds
  **directory junctions** into that repo — edit the repo copy, then Ben commits via GitHub
  Desktop. On easteregg2 skills are copied folders under
  `C:\Users\Ben Bolen\.agents\skills\`.

## Project purpose (priority order)

1. **Dissertation first** (deadline: this week, Sept 2026). LaTeX source lives in
   `Documentation\Reports and Papers\Dissertation\` (repo-tracked as of 2026-09-08; the
   `ProofFinal\` folder is the working copy and `upload\` mirror is the canonical Overleaf
   copy). Use the `latex-overleaf` skill for Overleaf work (installed on this machine via the
   ZCode_Skills repo junction).
2. Design and control of bipedal humanoid robot legs with artificial muscles (PAMs/BPAs)
   controlled by a synthetic nervous system. Lab: AARL (Agile and Adaptive Robotics Lab), PSU.
3. **Xi-correction-factor program** — run minimizers against pinned-knee test data in
   `Testing_Data\` to correct the rigid-body + PAM/BPA model:
   - **A. Flexor, pinned knee:** measure torque, compare to prediction (scripts in
     `Mesh_Optimization\Results\` + `Robot_Data\`), solve for Xi0–Xi2. Train on some tests,
     validate on others, discard the rest; keep the "best" result. *Open question: whether to
     revisit — different result pick, different transformation matrix, brackets on both sides
     of the knee, or a new Xi3-style term for stiffness/compliance along the force vector,
     independent of force orientation.*
   - **B. Extensor, pinned knee (muscle wraps over the joint):** solve for Xi0 and Xi3.
     Xi3 = loss of usable length as a function of arc length and the square of the additive
     complement to strain. Two physical brackets exist but the model approximates them as one;
     the bracket not used in calculation is mechanically constrained to limited −x deflection.
     Extensor and flexor Xi0 share meaning but differ in value.
   - **C. Biomimetic humanoid knee, flexor:** validation using Xi0–Xi2 from A. Original
     prediction showed sufficient torque; corrected prediction shows insufficient. One bracket.
   - **D. Biomimetic humanoid knee, extensor:** validation with Xi0–Xi3. One bracket.
   - *Open question: switch the iii solver from gamultiobjective to surrogateopt?*
4. **Mesh_Optimization rework** (`Opt_run.m` flexor / `Opt_run_Ext.m` extensor, with Morrow):
   use the Xi values to redesign the mechanics so each configuration (2× 20 mm BPAs each)
   meets or exceeds the magnitude of select human monoarticular muscles. Human torque and
   moment-arm data come from OpenSim (Gait2392; see directory map). Results land in
   `Mesh_Optimization\Results\` and feed scripts there used for hand-tuning or W.A.G.-ing a
   BPA prediction. Plots must be journal-publication ready. Legacy human calculations kept
   for now. *Open: need higher Xi1 — options are single-bracket frame orientation change, or
   a second bracket on the hip frame (unit direction + stiffness array order specified by Ben,
   still Xi1/Xi2, with the 2nd bracket's compliance projected on the force path u).*
5. Run tests on the human leg; show predictions from 4 meet or exceed human torque values.

## Directory map

- `Code\Matlab\Mesh_Optimization\` — **active workspace** (optimizers, route builders, results).
- `Code\Matlab\Functions\`, `Code\Matlab\Robot_Data\` — shared helpers; MonoPam classes/data.
- `Code\Matlab\Human_Data\`, `Code\Matlab\Bone_Mesh_Plots\` — human/bone data and plotting.
- `Code\Matlab\minimizers\`, `Code\Matlab\Previous Optimization Code\` — **legacy**; don't build on.
- `Code\Matlab\HX711-LoadCell\` — load-cell apps used with test data.
- `Code\MuJoCo_SNS\` — MuJoCo + SNS-Toolbox pipeline (see its README.md):
  custom BPA muscle (`bpa_muscle.py`, exact port of festo4/maxBPAforce/
  balanceX3), MuJoCo glue, SNS demo, `add_bpa_to_mjcf.py`. Env: conda
  `myoconv` (py3.10, opensim 4.5.2, mujoco 2.3.7, myoconverter editable from
  `Documents\GitHub\myoconverter` — carries local `[patched]` robustness
  fixes). Converted models live in
  `Solid_Models\OpenSim\Gait2392_Robotbody\mjc\`; `gait2392_robot.osim` there
  is GENERATED by `repair_robotbody_muscles.py` (mirrors Ben's edited
  right-side muscle paths onto the left, re-enables all 92 muscles) — don't
  hand-edit it, edit the repair script or robotbody and regenerate.
  - **`spinal/` — gait2392 spinal cord network (built 2026-09-09; see
    `spinal\DESIGN.md`)**: two-level RG+PF SNS network (406 neurons) for the
    converted simbody model — 92 MN pools + Ia/II/Ib afferents, per-leg
    half-center RG, 4 phase-shifted PF groups/leg, stance-gated Ib load
    sharing, descending DRIVE + balance inputs, solved standing posture
    (static-opt) injected as per-MN bias. Run with the `myo` env from
    cwd `Code\MuJoCo_SNS\spinal`: `check_rhythm.py` (rhythm alone — VERIFIED:
    0.885 s period, antiphase −0.88), `runner.py` (stand→walk→stand),
    `audit_signs.py` (anatomical sign audit, run DYNAMICALLY — static
    `actuator_moment` misses the pathpoint equality couplings),
    `diag_nan.py`, `fit_synapses.py` (IK/SO activations → NMF synergies +
    per-phase NNLS back-solve scaffold). Model repairs live in
    `apply_harness()` in runner.py, all audit-verified: hip_flexion/
    hip_adduction hinge axes WERE flipped vs OpenSim (negated); rect_fem
    lacked a patella wrap and pulled the knee into flexion (rerouted over
    the vastii's vas_med-P4 patella-tracking via point); quad_fem/gem/peri
    pruned (`PRUNE_MUSCLES`, Ben's list). DO NOT weld the conditional
    pathpoints (freezes moment arms — the earlier attempt did); they stay
    equality-coupled with `boundmass`/`boundinertia`. Gotchas: MuJoCo
    muscle Fmax = `actuator_gainprm[:,2]`; Ia/Ib afferents must be
    pure-signal (resting tone drove constant reciprocal inhibition);
    speed = DRIVE + presynaptic reflex-gain modulation (`params.MOD`).
    **Open blocker:** leg-DoF NaNs during sim (dt=2 ms, rigid pelvis rig) —
    dynamics not kinematics; suspects + order in DESIGN.md (mesh foot
    contacts first). Lit grounding: Rybak/McCrea RG+PF, Bunz 2026 (reflex
    speed control), Ben's Zotero "Sensory Afferent Database" collection;
    read Di Russo/Ijspeert/Bouri 2023 JNE before any novelty claims.
- `Code\Arduino\`, `Code\Festo\` — embedded/valve hardware code.
- **Xi1/Xi2 semantics (Ben, 2026-09-07):** they are *effective system-stiffness parameters*, not
  literal bracket beam stiffness — the fitted compliance lumps in the bracket, fixtures, and the
  cable winch on the test mechanism. "Bending" (Xi2) is modeled as a simple Hooke-law spring in
  N/m with no length/EI dependence. Ben's shorthand "bracket axial/bending stiffness" is
  convenience, not definition; don't over-interpret them structurally.
- Bracket reference points (Pbr, Pbri, Pbr2) are interpretation choices Ben sets from the CAD
  (Onyx FDM parts, not simple beams). Moving a point changes identified stiffness dramatically
  (Pbr2 move swung flexor Xi1 by 16x). Current active points: flexor insertion Pbri
  [-27.5,-107.81,-0.54]mm, flexor origin Pbr2 [-52.61,0,75.06]mm (pinned-flexor ONLY), extensor
  Pbr [-3.84,-46.44,62.5]mm (rib midpoint, medial side). Both 2brk brackets use two-rotation
  (Z then Y) frames — **proven 2026-09-08: frame convention is irrelevant** for y-symmetric
  stiffness arrays K=[a,b,a] (compliance is invariant under the y-axis second rotation;
  1trans/2trans agree to ~1e-17; transMode flag in minimizeFlxPin2brk demonstrates this).
  The 39x Xi1 swings between flexor CV runs (5.4e4 vs 2.1e6, with Xi0 +7.5mm vs ~0) are
  **gamultiobj stochasticity on a flat likelihood valley** — the pinned-flexor data cannot
  distinguish them (both fit mean RMSE ~1.6). Break ties by cross-configuration consistency
  or Lm_p/Lm_h checks, not reruns; consider rng seeding for reproducibility.
  Advisor requirement: one (Xi1,Xi2) consistent across pinned-flexor,
  pinned-extensor, and both biomimetic configurations.
- `Testing_Data\` — **important.** Immediate subfolders `2022_02_Festo\` and `2026_06_Festo\` matter.
  - `2022_02_Festo\` holds the Xi-minimizer family: outer CV drivers (`minimizeFlxPin10mm.m`,
    `minimizeFlxPin10mmX3.m`, `minimizeExt10mmX3.m`) call inner evaluators (`minimizeFlxPin`,
    `minimizeFlxPinX3`, `minimizeExtX3`, `minimizeExt`, `minimizeFlx`) that each carry their own
    `computeForceVector`/`Lok`/`fortz`. New (Sept 2026): two-bracket flexor method —
    `minimizeFlxPin2brk.m` (evaluator; the 6th `transMode` arg / the driver's
    `FLX2BRK_TRANS` env selects the convention — the laptop session's separate-file 1trans
    evaluator was DELETED 2026-09-08 at Ben's ruling: transMode is the vehicle of record.
    Insertion bracket: two-rotation (Z then Y) frame, K = [X1,X2,X1] (Sept-8-morning arm-1
    runs used [X1,X2,X2]). Origin bracket at `Pbr2` (current: [-52.61, 0, 75.06]/1000,
    Ben-set): 1trans pitch-only frame → K2 = [X2,X1,X2]; 2trans two-rotation frame →
    K2 = [X1,X1,X2] (Ben, late 2026-09-08). `USE_BRACKET2` flag).
    **Pbr2 applies ONLY to the pinned-knee flexor configuration** — the extensor evaluators
    (minimizeExtX3, minimizeExt) have their own independent bracket offsets; do not port Pbr2.
    + driver `minimizeFlxPin10mm_2brk.m` (env `FLX2BRK_MODE`
    = smoke|full, `FLX2BRK_SOLVER` = gamultiobj|surrogateopt, `FLX2BRK_TRANS` = 2trans|1trans),
    **compile-loop bug fixed 2026-09-09** (`for i = 1:numBPA` → `1:numel(results_cv)`): with ALLBPA
    overrides the old loop silently dropped fully-computed folds from the pooled front (noT3 runs
    kept 4 of 6 folds, full-5-test runs 5 of 10; harmless in the legacy 2-test era where
    numBPA == #folds). Existing mats can be recompiled post-hoc from stored `results_cv` — no rerun
    needed. Fold-structure findings: `Dig_foldLeverage.m` + log/mat in `Dig_out\` (2026-09-09):
    all folds land on the same score plateau (fold choice doesn't change conclusions); pooled
    fronts are ~80% cross-fold duplicates; predictability 46cm≈41cm < 40cm-tendon < 48cm≈47cm
    (47cm now normal after the encoder fix); Xi2 binds at ub=2e4 in recentered-bounds mats;
    1trans-vs-2trans and old-vs-new-bounds pooled fronts barely overlap in x-space (~0-3%).
    harnesses renamed to Ben's Collect/Dig scheme (2026-09-08): `Dig_crossPredict.m`
    (flexor→biomimetic/extensor cross-prediction), `Collect_ExtPinX3_sweep.m` (high-Xi1/Xi3
    hunt, pool {1,2,5,6,7,8}; tests 3/4/9 EXCLUDED per Ben),
    `Collect_batch.m` (overnight sequence; solver A/B settled 2026-09-07 — gamultiobj kept,
    the env flag is still honored); `Dig_*` = analyze existing mats (Dig_CVpatterns,
    Dig_ExtPinX3_CV, Dig_ExtPinX3_Xi3map, Dig_FlxPin_2brkt); their logs/one-off outputs live
    in `2022_02_Festo\Dig_out\`. `Robot_Data` must be on the MATLAB path for the biomimetic
    evaluators.
    2026-09-08 late rerun (commit 19fea63): driver bounds re-centered on the biomimetic
    hand-tune winner (Xi0 +12 mm, Xi1 5e5, Xi2 1e4) → Xi1 ∈ [3e4,1e6], Xi2 ∈ [5e3,2e4],
    initial population 0.5–1.5 cm / 5e4–5e5 / 7e3–1.5e4. Current CV mats:
    `minimizeFlxPin10_results_20260908_2brkt_{1trans,2trans}_{noT3,noT3noT5}.mat`
    (BPA #3 = 47 cm EXCLUDED from training; #5 = 41 cm also dropped in noT3noT5; flexor
    labels: 1–5 = 48cm, 46cm, 47cm, 40cm-tendon, 41cm). All T3-including fronts (20260907
    + full/noT5 variants) archived in `Dig_out\old_T3_results\`; log
    `Dig_out\encoder_campaign_noT3newXi_20260908.log`. Biomimetic-flexor chain
    `Dig_FlxBio_dubfilt` → `_handtune` → `_refine` (mats/logs in Dig_out; hardcoded D:/
    paths, written on easteregg2). **Scope of the 1trans≡2trans proof:** y-symmetric arrays
    only (K_x=K_z); the per-convention K2 orderings are Ben's buckling argument, NOT covered
    by the proof — that choice is carried, not derived.
  - **Laptop-session mining handoff (2026-09-08, merged and adjudicated):** findings in
    `2022_02_Festo\HANDOFF_laptop_20260908.md` — no flexor test is garbage; Xi2 is
    consistent across configs (~1.1e4) while Xi1 is the flat, undetermined one; Xi0 pins
    at lb=0 on the front (Ben may want lb<0). Kept under the Collect/Dig scheme:
    `Dig_FlxPin_2brkt_picksScan.m` (whole-front pick scoring + shortlist cross-prediction)
    and `Dig_allbpaNumHoldScan.m` (per-test held-out table + E-fold numHoldout spectrum);
    both still UNTESTED end-to-end. Deleted as superseded: the `_1trans` evaluator/driver,
    `nightBatch_20260908.m`, `mine_smoke_20260908.m`. The "_smoke mat" question is CLOSED:
    that file is `minimizeFlxPin10_results_20260907_2brkt_2trans.mat`, now archived in
    `Dig_out\old_T3_results\` (canonical restored-2trans full CV; no TRANSMODE field is
    expected for its vintage) — never resurrect a smoke-named copy.
  - `2022_02_Festo\` subfolders (Flx/Ext × 10mm/10mm_pinned/20mm/40mm): needed only by certain
    plotting functions — add to path transiently and **remove afterward**; they shadow other
    plotting functions if left on path.
  - `2026_06_Festo\` subfolders (Ext_10mm_pinned, Flx_10mm_pinned, Torque_Pressure_Tester V1–V3,
    ValveDataAcquisition): long raw data files, summarized in tabs of
    `2026_06_Festo\Results_table_10mm_pinned.xlsx` — read that instead of the raw files.
  - Also holds large BPA max-force / force-pressure-contraction characterization structures;
    some processed here, some from the Muscle_Sensory repository. Used with `Code\Matlab\Functions`
    and the HX711 apps.
- `Solid_Models\Biomimetics_2022-Knee_Test\` — the knee test setup (CAD, STLs, point clouds).
- `Solid_Models\OpenSim\Gait2392_Robotbody\` — OpenSim human reference models/data (incl.
  Bifemsh/Vastus-adjusted variants used for muscle torque targets). `gait2327.osim` there is
  **GENERATED** (2026-09-08, `D:\GitHub\myoconverter\build_gait2327.py`): stock gait2392_simbody
  with muscle GeometryPaths (points+wraps) ported from robotbody via XML surgery — stock Thelen
  params, all 92 muscles enabled; verified paths == robotbody 92/92. Key finding behind it:
  **robotbody is a LEFT-LEG reference model** — left-leg paths are byte-identical to stock; the
  right leg has 25 muscles with moved paths (15 with different point counts, e.g. ercspn_r 2→9)
  onto shared placeholder points, and 19 right-side muscles `appliesForce=false`. Stock params +
  those right-leg paths = path lengths off 7–480% → OpenSim equilibrium fails (tfl_r) and force
  maps would be garbage. MyoConverter outputs (`Gait2392_Robotbody\myosuite_gait2392_robotbody\`,
  `Solid_Models\OpenSim\myosuite_gait2392_simbody\`) are **regenerable** (~116 MB, Geometry =
  stock OpenSim STL copies) — gitignored 2026-09-08; rebuild with MyoConverter instead of
  restoring from git.
- `Documentation\Reports and Papers\` — papers; `Dissertation\` = the dissertation (see priority 1).
  `Documentation\Research Notes\` — agent-written lit reviews with V/L/U provenance tags
  (board-sports/gymnastics motor control; kickflip neuromuscular hypothesis).
- `Neuromechanical_Models\` — Animatlab CPG walker models. The `_Standalone.asim` exports
  are STALE (predate the subsystem reorg); the instrumented phase-1 file is
  `Biped_2xCPG_wSubs\walk new new tester added 2 axis_phase1.asim`; Ben exported a fresh
  `Biped_2xCPG_wSubs_Standalone.asim` 2026-09-09 (runs clean headless; has RG charts).
  DataTool_*.txt are run byproducts. **2026-09-09 .aproj surgery (GUI-verified clean open;
  ACTIVE WORK — scheduled 8am PT session continues):** RH-side wiring completed (RH_RG
  half-center, RG→PF drive, hip Ia/II + ankle Ib/II coordination — mirrored from LH),
  LH↔RH RG commissural inhibition (4 OffPages on the top page + 2 nA tonic on L RG ext),
  4 biarticular muscles per leg added (Gas/BFlh/Semimem/RF: Gait2392 Fmax 2241/896/1288/
  1169 N applied; full Deng-style chains from the existing PF layers), RH knee/ankle
  muscle-drive adapters fixed (were Renshaw-driven, now MN-driven), grids off, all page
  drawings rebuilt. **Ben added connections manually and reports some are STILL missing —
  audit drawn arrows vs the standard section FIRST.** Also open: RG does not oscillate
  (latches ≤10 nA; sweep 20–40 nA / pulse / Deng Table A2 conductances), 16 placeholder
  attachments (then RestingLength = TSL+OFL), no Renshaw/II on new muscles, TFL needs an
  abduction DOF, LH_HipZ→LH_Hip rename pending. **Toolchain + format traps:**
  `Biped_2xCPG_wSubs\tools\` (build_A2→build_B2→repair_dup_ids2→fix_gas_ia3→
  rebuild_pages_v3→move_links_by_page→set_fmax2; pipeline from pristine git 6437441
  backup — re-run the WHOLE chain, never patch a patched file). Traps: page CDATAs are
  NOT covered by whole-file XML validation (validate each page separately); tag
  depth-counters must skip CDATA and match `<Node attr>` opens; page rebuilds must emit
  `<Version>` + closing tags; cloned blocks need child-object GUIDs re-rolled; **NEVER
  hand Ben an .aproj without opening it in AnimatLab2.exe first** (zero Error dialogs =
  pass; each broken page throws one dialog) and never save from the GUI (taskkill /F).
  Full brief: `CHATGPT_HANDOFF.md` AnimatLab section.
- `Code\Matlab\SNS_Simscape\` — SNS neuron block library (`SNS_Library.slx`: non-spiking RC
  neurons, E/I synapses with AUTO E/I icons, Ia/Ib afferents) + `KneeReflexDemo.slx`
  (antagonist BPA knee reflex demo; runs in plain Simulink). **Diagram conventions
  (Ben, 2026-09-09):** open circle = neuron, white triangle = excitatory, solid black
  circle = inhibitory, ellipse = muscle; tints = Okabe-Ito CVD-safe; markers auto-draw
  from the sign of Esyn. Journal figures regenerate via `sns_draw_circuit.m` (circuit
  redraw) and `sns_function_subnetworks.m` (Szczecinski 2017 arithmetic primitives:
  add/sub/div/mul/diff/integrate) into `figures\`. **CAD route is now URDF** (Ben's
  choice): SolidWorks → sw2urdf add-in → .urdf → `smimport`. The old "Simscape
  Multibody not licensed" blocker was FALSE — the license carries it as legacy feature
  `SimMechanics` (=1) and `sns_urdf_smoke.m` proved smimport works (2026-09-09).
  sw2urdf v1.6.1 NOT yet installed (installer in Downloads; official build targets
  SW2021, so on SW2025 watch for the vanishing-dialog issue, issue #147); fallback:
  Simscape Multibody Link IS installed+registered (disabled — enable in SW Tools >
  Add-Ins). `import_simscape_when_ready.m` takes URDF or XML. Tendon parts
  (`Tendon_Extensor/Flexor.SLDPRT`) are built but NOT yet inserted into `09_BA_003`;
  GUI steps in the README. Gotchas: mask icon drawing commands take numbers only
  (no LineSpec/name-value; `color('black')` not RGB), and keep masked blocks square
  so circle icons stay round.
- Repo root: `CHATGPT_HANDOFF.md` (brief for other AI assistants when ZCode is unavailable)
  and `CHATGPT_REPORT.md` (their report back; 2026-09-08 edition covers Overleaf
  manuscript-status edits) — keep both current when work is handed off.

## OpenSim / MyoConverter / SNS-Toolbox on easteregg2 (Sept 2026)

Distinct from the laptop's `myoconv` env (above): easteregg2 uses conda env **`myo`**
(py3.10). Working pin set: conda `opensim=4.4.1=py310np121` (opensim-org channel) +
pip `mujoco==2.3.7` + pip `vtk==9.2.6` + pip `trimesh==3.23.5` + pip
`numpy==1.21.6` (**must be the pip/OpenBLAS build** — conda's MKL numpy crash-planes
with 0xc06d007e in np.dot once opensim+mujoco are loaded) + `sns-toolbox` (--no-deps,
plus CPU torch --no-deps + filelock/typing-extensions/sympy/networkx/jinja2/fsspec +
graphviz). Run python with `CONDA_PREFIX=D:\Anaconda\envs\myo` (mujoco reads it).
Second env `opensim` = py3.11 + pip opensim 4.6 wheel (PyPI wheels are cp311+ only on
Windows) for model editing / opensim-cmd follow-ups. IK output:
`Documents\OpenSim\4.6\Models\Gait2392_Simbody\subject01_walk1_ik.mot` (workflow:
Scale first — the bundled IK setup consumes its `subject01_simbody.osim` output).

- MyoConverter clone: `D:\GitHub\myoconverter` (API: `O2MPipeline(osim, geometry,
  out, **kwargs)`). Drivers: `convert_gait2392.py` (stock→MuJoCo, ~90 min, DONE),
  `build_gait2327.py`, `sns_gait2392_cosim.py` (SNS↔MuJoCo co-sim, proven),
  `compare_pathpoints.py`/`verify_gait2327.py`. Reference CPG: 
  `D:\GitHub\Two_layer_CPG_SNS_Toolbox` (Nourse 2023; `build_net(dt)` runs verbatim).
- **Windows pitfalls (each cost real time — do not relearn):** MyoConverter step 3
  uses multiprocessing spawn → every driver script needs `if __name__ == "__main__":`;
  guard-less children re-run the whole pipeline and hold the output log open
  (PermissionError WinError 32) — kill zombie python children before rerunning.
  Classic conda solver hangs 20+ min on the upstream yml — use libmamba, and prefer
  pip for mujoco/vtk/numpy so conda only carries opensim.
- **SNS-Toolbox 1.5.2 wheel ships two parallel class trees**: use `Network` from
  `sns_toolbox.networks` (has `.compile(backend='numpy')`) with neuron/synapse classes
  from legacy top-level `sns_toolbox.neurons` / `sns_toolbox.connections`
  (`resting_potential`, `e_lo`/`e_hi` kwargs). `design.neurons` instances FAIL the
  isinstance check inside design/networks.py.
- Co-sim pattern (proven, ~4x real time): SNS dt=0.1 ms → MN voltage →
  `stim2activation()` sigmoid → `mj_data.act[actuator_id]` → `mj_step` →
  `-gain*actuator_force` → Ib feedback (15 SNS inputs: 3 drive + 6 Ia + 6 Ib).
  Converted-model caveats: MJCF has `limited="false"` joints (the pelvis-pinned
  variant enables limits), keyframe qpos matches the FULL model (strip keyframe if
  you remove joints), and MJCF resolves `Geometry\` relative to the XML location.

## MATLAB workflow (follow exactly)

- Entry scripts are in `Code\Matlab\Mesh_Optimization\`: run `Opt_sanity.m` / `Opt_sanity_Ext.m`
  **before** the full `Opt_run.m` / `Opt_run_Ext.m`.
- Every run starts from a context struct: `ctx = buildKneeFlexorContext20mm()` or
  `buildKneeExtContext20mm()`; objectives/constraints are wrapped as anonymous functions of
  `(x, ctx)`. Flexor and extensor use different constraint machinery.
- Optimizers use `parpool(10)` + `surrogateopt` (≈7000 evals) then `patternsearch` (≈15000) —
  **hours, not minutes.** (parpool(10) fits easteregg2; on DESKTOP-5Q16KE9 cap at 6.) Run via
  the `matlab` skill (`matlab -batch`) in the background with
  output tee'd to a log; check the log tail, never babysit the run turn-by-turn.
- Path setup: self-locate the repo root, then `addpath(genpath(root/Code/Matlab))` with
  **Mesh_Optimization winning any shadowing contest** (duplicate filenames exist in Robot_Data,
  Knee_Torque_revision_3, Bone_Mesh_Plots). `Debug_RouteElim_Ext.m` shows the canonical setup.
- Results (.mat/GIF/CSV) go to `Mesh_Optimization\Results\`.
- **2brk drivers keep Ben's "%% Pick best solution (later, flexible)" section VERBATIM** (lowercase
  `pick`, commented `sol_actual` lines, hand-editable hardcoded k1/k2/k3, the two disp tables and
  two Mean fprintf lines) — only the evaluator call line may carry the extra evaluator arguments.
  Do not rework it into a PICK-strict or mean-of-all-tests version.
- Never casually re-run entry scripts; ask/state intent first — runs are long and results get
  overwritten.

## Current state (Sept 2026)

- Flexor: **solved and feasible** with 2 BPAs, but currently assumes double force at the same
  attachment points instead of mirrored routes (commit 0453a2f caveat — still true).
- Radius controls (commit fd5963d): `ctx.bpaRadiusMode` = "scalar" (geo.bpaRb / geo.bpaRs) or
  "bpaR" (candidate-dependent physical radii); `geo.bpaRbOffset` / `geo.bpaRsOffset` defaults 0.
  Three-radius scheme: `bpaRb` = 20 mm nominal pWrap standoff, `bpaRs` = 16 mm min tibia
  clearance, `wRap` = 25 mm Xi3 bend-length radius.
- Extensor route elimination **REWORKED AND VERIFIED (Sept 2026, commits through 416f08f)**:
  - Ben's rule in `buildDistalRingLocation20mm.m>candidateEliminationTest`: femur rows p2:p5
    tested in femur frame (vectors from previous active row, rotated +90°; remove when big more
    CCW); tibia rows p6:p8 in t1 frame (vectors from next active row, rotated −90°; remove when
    big more CW). ONE `atan2(cross,dot)` on the rotated pair — comparing two principal atan2
    values separately IS seam-bitten at ±180°. p3/p4 anchor to p7 while active; cascade guard
    forces release order p5,p6,p4,p3,p7; ONE removal per step (priority must init
    `bestMargin = -inf`, NOT 0 — hysteresis margins are negative).
  - Collision gates ON: chord must clear envelope by `geo.bypassTol` (0.5 mm), 2 mm endpoint
    trim, contraction-radius slack `geo.bypassRelaxFemur`/`bypassRelaxTibia` = 3 mm (envelopes
    use inflated 19.25 mm radius; routed BPA is contracted per `bpaR`), p7's gate = the two
    tibia cylinders only (wall band is virtual). `geo.marginTolD` = 1.5° release hysteresis —
    releases BEFORE the wrapped contact goes collinear and collapses the moment arm (the old
    "torque jog"). Tibia seeds CONSTRUCTED: p6 = +30° on upper clear circle, p8 = −30° on
    lower, p7 = (tibiaWallX, mid-y); knobs `geo.seedP6AngleD/seedP8AngleD/seedP7X`.
  - Verified x0 schedule: p5@−91.1°, p4@−50.4°, p6@−47.8°, p7@−20.2°, p3@+6.1°, one per step,
    min moment arm +45.2 mm, max torque step 0.53 N·m (jog eliminated).
  - `Debug_RouteElim_Ext.m` (tracked) = ~1 min verifier: elimination table, gate audit, p1/pEnd
    bound scans. Tiled route figures in Opt_sanity_Ext/Opt_run_Ext: 4 poses per row, legend in
    a 5th column spanning all rows.
- `MonoPamDataExplicit_balanceX3.m` has exactly ONE copy: `Robot_Data\` (`mif = BPAcount*Fmax`).
  The nested Knee_Torque_revision_3 shadow copy was deleted — do not resurrect it.
- Open design question: replace the 4 per-contact wrap terms (geometricBendMeasure) with p2's
  term kept + ONE unified R×(total polyline turn) term for p3–p8, plateau-calibrated so the
  Xi3·bend product stays consistent; prototype behind `ctx.bendModel` flag. Ben hasn't decided.

## Readmes

- `Code\Matlab\Mesh_Optimization\Mesh_Optimization_Readme.md`, `Code\Matlab\Functions\Functions_Readme.md`,
  `Code\Matlab\Robot_Data\Robot_Data_Readme.md` — agent-readable copies; **the .docx versions are
  canonical for humans.** If you edit a .md copy, remind Ben to update the .docx.
- `Code\Matlab\Mesh_Optimization\Knee_Torque_revision_3\Knee_Torque_revision_3\README_revision_3.md`
  — detailed flexor BPA/route model doc (modes, how to run tests).

- Known data concern RESOLVED (Ben, 2026-09-08): the angle-shifted test was the **flexor
  47 cm test** — its encoder read ~5.3° low. `minimizeFlxPin2brk.m` adds +5.3° to that
  test's experimental angles at build time (**Angle only — phiD is NOT shifted**; an
  earlier version mistakenly shifted both, fixed late 2026-09-08). Any other script that
  refits the 47 cm test must apply the same Angle-only shift; pre-2026-09-08 fit results
  predate the correction.

## Hazards — do NOT open these as source

- Known data concern RESOLVED (Ben, 2026-09-08): the angle-shifted test was the **flexor
  47 cm test** — its encoder read ~5.3° low. `minimizeFlxPin2brk.m` adds +5.3° to that
  test's experimental angles at build time (**Angle only — phiD is NOT shifted**; an
  earlier version mistakenly shifted both, fixed late 2026-09-08). Any other script that
  refits the 47 cm test must apply the same Angle-only shift; pre-2026-09-08 fit results
  predate the correction.
- `Solid_Models\Biomimetics_2022-Knee_Test\Point_cloud\Tibia_copy.txt` (7.1 MB point cloud);
  `Spine_Mesh_Points.txt` (172 KB, duplicated in 3 places); `HX711*sempio.txt` (1 MB);
  any `.mat` in `Previous Optimization Code\Trial Results\` (up to 95 MB).
- `.asv` files = stale MATLAB autosaves that mislead; the real source is the `.m` beside them.
- Logs: grep/tail them, never read whole. Point clouds/binary .mat: process with MATLAB, never
  ingest.

## Agent efficiency rules (credit/token hygiene)

- Delegate broad searches and audits to subagents; keep main-session context small.
- Fresh session per task; long MATLAB runs in background; plan mode before expensive execution.
- Ben uses **GitHub Desktop** for git and does not know git/Git Bash CLI. Read-only git
  (status/log/diff/fetch) via my own tools is fine, but **never commit, amend, merge, or
  push without his explicit go-ahead for that exact action** — even when he asks for git
  work, prefer preparing the changes + a draft commit message and letting HIM click
  Commit/Push in GitHub Desktop. Give GUI steps in GUI terms (Fetch/Pull origin, file
  checklist, commit message, Push origin) — never hand him shell git commands. (Rule set
  2026-09-08 after a two-machines-one-branch mix-up between parallel AI sessions.)
- Prefer COM/API automation (SolidWorks skill, Overleaf file edits, Zotero local HTTP) over
  screenshot-driven GUI automation.
- Model strategy: default GLM-5.3-Flash; escalate to GLM-5.3 only for hard debugging or
  architecture decisions. Peak billing is Mon–Fri 11 PM–3 AM Pacific — schedule expensive
  (GLM-5.3) work before 11 PM or on weekend nights; evenings/weekends are 50% off-peak.

## Extensor campaign state (Sept 2026, desktop session — merged here for the laptop session)

- **Extensor Xi0 is NEGATIVE**: driver bounds [-2,0] cm; transfer from flexor = `-g(1)`
  (minimizeExt10mmX3.m). Screens that pass positive flexor Xi0 into the extensor are wrong
  (bit once, Sept 9).
- `minimizeExt10mmX3.m` now loads its Xi1/Xi2 lock from
  `minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat` (corrected-encoder 2trans flexor
  front; pick pair 4.35e4/1.7e4) — NOT the legacy 20260730 mat.
- **Pbr variants recovered from git** (minimizeExtX3.m): rib midpoint [-3.84,-46.44,62.5]
  (current, ACTIVE), lower bolt hole [-6.26,-29.69,75.06] (the ORIGINAL "lower bolt hole";
  [-2.65,-54.71,75.06] later inherited that comment and is now "lower section" — do not
  confuse them), midpoint between two bolt holes [-7.325,-22.27,75.06] (EXCLUDED per Ben).
  minimizeExtX3.m has transMode (arg 6 / env EXTX3_TRANS; 1T≡2T for its y-symmetric
  K=[X2,X1,X2], verified). Rib midpoint vs lower bolt hole CVs (1trans, fronts captured)
  are essentially identical: pick Xi3 0.62 vs 0.60, Xi0 -10.1 vs -11.6 mm, bio-ext
  2.195/2.35/4.73 vs 2.204/2.369/4.727 → extensor Xi identification is robust to the
  bracket-point choice. Filter pass 72/90 vs 83/90.
- **minimizeExt.m (BIOMIMETIC extensor) stays at its committed state** — own centroid
  bracket [8.38,20.75,25.1], K=[X1,X2,X1] (Ben-set), no transMode; do not "extend" it.
  Extensor drivers have NO save statement — capture fronts by appending save(...) to the
  launch line (pattern: minimizeExt10mmX3_front_*.mat).
- **Results under the new pair**: extensor front collapses to 4 unique members, Xi0
  -6..-12 mm, Xi3 0.50-0.78 (all >> 0.02), pinned RMSE 0.96-1.64, bio-ext RMSE 1.99-2.26
  (FVU > 1 there — known gap). Best bio-ext member: Xi0 -6 mm, Xi3 0.504. Flipped-screen
  note: smaller |Xi0| (-1.6..-3.3 mm) beat -5..-12 mm on the pinned pool; Xi3 grid capped
  at 0.35 in screens — all leaders hit the cap.
- **New Dig scripts**: Dig_ExtPinX3_screen.m (flexor candidates × Xi3 grid through both
  extensors), Dig_ExtPin_frontBio.m (extensor front → bio-ext), Dig_FlxPin_frontScan.m
  (flexor front → bio-flexor), Dig_FlxPin_cornerCheck.m (direct Xi check vs pinned
  baseline), Dig_FlxBio_refine2.m (STAGED, not yet run — extends the flexor refinement
  grid past its Xi0/Xi1 edges).
- **Extensor front reproducibility (2026-09-09/10):** relaunching the identical
  noT3newXi-pair extensor CV with a fresh GA seed reproduced the front EXACTLY — same 4
  unique members (Xi0 -6.0/-8.6/-10.1/-11.9 mm, Xi3 0.504/0.661/0.621/0.783), same
  ranking. The family is stable, not a stochastic artifact.
- **noT3 name note:** the desktop's current noT3 pair (NEW bounds + corrected 47 cm,
  `minimizeFlxPin10_results_20260908_2brkt_{1trans,2trans}_noT3.mat`) supersedes the
  old-bounds campaign pair, which is archived in `Dig_out\old_bounds\`. Watch for name
  collisions if the laptop session also produced 20260908 noT3 mats.
- **Results tracking (Ben, 2026-09-10):** Ben will create his own pivot table to track Xi
  results across evaluators/configurations. Do NOT build unprompted tracking tooling for
  it. When he defines the format, future sessions maintain/populate it; until then,
  present results as simple tables with all 3 GoF and always name the source .mat.
- **Opt_run_Ext with the new Xi pair (2026-09-10, desktop):** completed feasible —
  ctx (Xi0 -10.1mm, Xi1 4.354e4, Xi2 1.701e4, Xi3 0.621 from the 20260910_noT3 front).
  Pick: Xi0 -12.3mm, Xi3 0.281 (locked pair 1.448e4/1.355e4 = pass-1 wide-search match).
  Max path length 0.832 m, contraction 0.69 KMAX, torque margin met (thinnest +0.063% at
  -18.9 deg), shortfall penalty 0, binding constraint -0.000442 (near-active at optimum).
  Results: Mesh_Optimization\Results\Vas_Pam_20mm_Result_20260910_{0446,0528}.mat
  (XiUsed recorded inside). Log: 2022_02_Festo\Dig_out\Opt_run_Ext_newXipair_20260910.log.
  Launch recipe needs cwd = 2022_02_Festo (OpenSim Vasti txt) AND Code\Matlab on path
  (Colors.m) AND Mesh_Optimization on path.

## Xi values for the dissertation text (SETTLED 2026-09-10, Ben-approved picks — verified by direct .mat loads)

A separate chat is updating the dissertation text with these. Sources (Testing_Data\2022_02_Festo\):

- **Flexor, pinned 2brk**: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat, pick 107
  Xi0 = +8.9 mm (0.0089 m), Xi1 = 5.62e4 N/m, Xi2 = 1.85e4 N/m
  Pinned per-test RMSE 1.98/1.36/2.22/1.47/1.23 (48/46/47/40cm-t/41cm), all FVU < 0.16
  1trans twin (minimizeFlxPin10_results_20260908_2brkt_1trans_noT3.mat): Xi0 +5.95 mm, Xi1 4.36e5, Xi2 2e4
- **Extensor, pinned**: minimizeExt10mmX3_results_20260910_noT3.mat, pick 1
  Xi0 = -10.1 mm (-0.0101 m), Xi3 = 0.621; Xi1/Xi2 LOCKED to the flexor pair 4.354e4/1.701e4 (never searched)
  Held-out RMSE 1.04-1.64 vs baselines 2.44-3.62 (72/90 filter pass)
- **Signs matter**: flexor Xi0 positive, extensor Xi0 negative. Extensor Pbr = rib midpoint
  [-3.84,-46.44,62.5] (lower bolt hole [-6.26,-29.69,75.06] tested: equivalent, pick Xi3 0.60 vs 0.62).
- Mesh optimizations with these values: Opt_run_Ext completed feasible
  (Mesh_Optimization\Results\Vas_Pam_20mm_Result_20260910_0528.mat, XiUsed inside);
  Opt_run launched with buildKneeFlexorContext20mm.m re-pointed to the 2trans_noT3 pick 107.
- **FVU caveat**: bio-ext 52cm FVU 2.35 > 1 at this pick — improved vs baseline 4.49 but do not overclaim.
- **Opt_run runtime reality check (Ben, 2026-09-10):** configured eval budgets (surrogateopt
  ~7000 + patternsearch ~15000) are MAXIMA — actual stages exit early on FunctionTolerance;
  Ben's observed patternsearch stages wrap in ~5 min, not 1.5-2.5 h. Extrapolate ETAs from
  measured rates and his historical runtimes, not from configured caps.
