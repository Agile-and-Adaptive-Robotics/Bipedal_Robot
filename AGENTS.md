# AGENTS.md — Standing context for AI coding sessions in Bipedal_Robot

Loaded automatically at session start. Keep it current; keep it lean.

## Machines

- **EB475WS4 (third machine, added 2026-09-10)** — repo at `D:\Github\Bipedal_Robot`,
  has a D: drive but **no `D:\Anaconda`** (the easteregg2 notes below don't apply here).
  Anaconda base = `C:\ProgramData\anaconda3` (conda 26.5.3; NOT writable — env creation
  fails under ProgramData and `D:\Users\...`). Spinal SNS env = **conda env `myo` at
  `C:\Users\Ben Bolen\.conda\envs\myo`** (py3.10.21; pip numpy **1.22.4** — scipy 1.9.3
  needs the 1.22 ABI, the easteregg2 pin 1.21.6 fails with "compiled against API version
  0xf"; + mujoco 2.3.7 + sns-toolbox 1.5.2 --no-deps + torch 2.14 cpu --no-deps +
  filelock/typing-extensions/sympy/networkx/jinja2/fsspec/graphviz/tqdm + scipy 1.9.3 +
  matplotlib 3.7.5). Call as `C:\Users\Ben Bolen\.conda\envs\myo\python.exe`. Second env
  `MyoSuite` (py3.9, mujoco 3.1.2, MyoSuite 2.8.5, NO sns_toolbox) at
  `D:\Users\Ben Bolen\.anaconda3\envs\MyoSuite`. No git CLI on PATH (Ben uses GitHub
  Desktop). **Inline `python -c "multi-line"` gets EATEN by this shell — write script
  files instead.** Hardware/installs (from the Sept-9 stash, still true): 12 logical
  cores, 32 GB RAM — **not a heavy-runs box** (cap `parpool` ~6–8; long optimization
  runs still belong on easteregg2). MATLAB **R2025a Update 1** at
  `D:\Program Files\MATLAB\R2025a` (on PATH); PSU network license: `license('test')`
  readings are unreliable here — trust a runtime test (Simulink/Simscape/Optimization/
  Parallel/CurveFitting/Statistics and gamultiobj/surrogateopt/patternsearch all verified
  RUNNING). **Simscape Multibody (`smimport`) is BLOCKED here** — the license server
  reserves all SimMechanics seats (error 101.2), so the URDF→smimport CAD route cannot
  run on this machine (it works on the laptop). SOLIDWORKS 2025 SP3 (33.3.0) at
  `D:\Program Files\SOLIDWORKS Corp` (contains a duplicate `SOLIDWORKS (2)` folder, same
  version); SW2URDF add-in registered (enable via Tools > Add-Ins); Simscape Multibody
  Link NOT installed.
  **2026-09-20 rebuild notes (after Ben replaced the C: drive; registry links to
  D:-installed software were wiped):** VS Code 1.138 at `D:\Program Files\Microsoft VS
  Code` — extensions (python/jupyter/matlab/latex-workshop/cpptools/xml) + settings
  rebuilt; project `.vscode\` (untracked) + Desktop shortcut "Bipedal_Robot (VS Code)";
  python interpreter = the myo env. **SW COM was broken** (LocalServer32 used 8.3 short
  paths; 8.3 generation is disabled on D:) — fixed per-user via
  `HKCU\Software\Classes\CLSID\{6AF263BB-EB9F-4176-89E9-4F892EB0CA3D}` → real path; API
  verified (33.3.0) through the solidworks skill with **pywin32 now pip-installed in the
  myo env**; `sw_session.py` gen_py fix (pip pywin32 caches in %TEMP%\gen_py — the old
  site-packages guess failed; REPO COPY EDITED, awaiting Ben's commit). MiKTeX 22.1
  (`D:\Program Files\MiKTeX`) was dead (its config died with the old registry) — removed;
  complete package set installed **user-scope at `D:\MiKTeX`** (pdflatex
  `D:\MiKTeX\miktex\bin\x64`, on user PATH; Ben declined the all-users/admin install for
  now — staged one-click: `D:\temp\miktex_setup\elevated_install2.ps1` + local repo
  `D:\temp\miktex-repo`). MATLAB R2025a U1 re-verified headless (Optimization + Simulink
  license OK).
- **DESKTOP-5Q16KE9 (laptop), the main workstation** — newest/premium
  MATLAB + SolidWorks, but modest hardware (6 cores, 16 GB RAM). MATLAB **R2025b** at
  `C:\Program Files\MATLAB\R2025b`; SOLIDWORKS **2025 SP4.1** (33.4.1) at
  `C:\Program Files\SOLIDWORKS Corp`. Its GitHub repos live under
  `C:\Users\Ben\Documents\GitHub\` — **no D: drive on this machine**. Conda =
  `C:\Users\Ben\.anaconda3` — env `myoconv` (py3.10.21, mujoco 2.3.7,
  sns-toolbox 1.5.2, numpy 1.25.2), call
  `C:\Users\Ben\.anaconda3\envs\myoconv\python.exe`; also `snsenv`.
  **MuJoCo↔Simulink bridge BUILT + PROVEN on R2025b here (2026-09-16; tests
  a/b/c bit-exact, E1/E2 reproduced — BRIDGE_REPORT.md "LAPTOP PORT" section;
  rebuild driver `mujoco_bridge\matlab\laptop_build_bridge.m`)**.
  Simscape Multibody license WORKS on the laptop for smimport AND hand-built
  models (E0 PASSES here; the block-ADD failure is EB475WS4-only) — Ben's
  cylinder-elbow plant belongs on THIS machine.
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
- Custom skills are version-controlled in the `ZCode_Skills` repo. **EB475WS4: the repo is
  at `C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills`** (NOT `D:\GitHub\` — corrected
  2026-09-16; easteregg2: `D:\GitHub\ZCode_Skills`; laptop: `Documents\GitHub\ZCode_Skills`).
  On EB475WS4 `C:\Users\Ben Bolen\.zcode\skills\` holds **directory junctions** into that
  repo (verified 2026-09-16 23:24 — the earlier "plain copies under `.agents\skills\`"
  note was wrong/stale; `.agents\skills` does not exist here). Live set, ALL junctions:
  animatlab, latex-overleaf, matlab, mujoco, neuroscience, opensim, python,
  sns-toolbox, solidworks. **NEW skills added 2026-09-16**: `neuroscience` (CPG circuit
  design: NaP half-centers/τh, laminated inhibition, commissural latch trap, curriculum
  tuning), `mujoco` (converted-model pitfalls: contact pairs ignore contype, ligament
  surrogates, muscle actuators), `python` (AARL workflow: inline-python trap, stdout
  wrappers, background chains, optuna hygiene) — and `sns-toolbox` is now IN the repo
  (updated with the backend-override trap + fixed-τh NaP measurements) and junctioned
  like the rest. `opensim` gained a .sto/synergy-analysis section. The repo still carries
  the retired `myoconverter/` folder (superseded by `opensim`). Edit the repo copy, then
  Ben commits via GitHub Desktop.

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
    **Blocker RESOLVED 2026-09-10 (EB475WS4 session; full chain in
    DESIGN.md):** pathpoint followers' keyframe-sampled range limits
    fought their equality couplings (knees jammed AND huge constraint
    forces → the NaNs; stripped + `limited="false"` + follower
    `armature="0.5"`), boundmass 0.01 was too light at human muscle
    forces (singular mass matrix → 0.1), predefined contact pairs ignore
    contype (patch the <contact> section instead), pelvis yaw had zero
    damping (rig now springs all 3 rotations), the standing solve was
    degenerate (solve on a ground-on model, joint rows only, preloads
    hip -40/knee -60/ankle +40), and the II-afferent BASELINE needed the
    stance gate (ungated i0_ii was a global ~0.2 co-contraction floor).
    **Knee convention (Ben-corrected, real-actuator-verified): the
    converter PRESERVED OpenSim's flexion-negative knee; knee_angle
    NEGATIVE = flexion, and the joint ships `limited="false"` (range
    inert). An earlier "knee range flip" was reverted the same evening.**
    audit_signs.py's moment-injection signs are inverted vs real
    activation (docstring caveat); use `_muscle_direction_test.py`
    (ctrl=1 per muscle) as ground truth. After the stack: deafferented
    suspended-air 22 s clean AND full ground stand→walk→stand 22 s with
    afferents, stayed up. **Night 2026-09-10/11 (Ben's staged plan):**
    (1) air-stepping lit (Ivanenko 2002): preferred 0.3 Hz cycles (3×
    SLOWER than walking, not 5× faster), E-duty 53%, sinusoidal joints;
    (2) deafferented air + `--no-interleg` WORKS: 21 s clean, knee
    −75..+15° true swing flexion (amplitude/frequency DECOUPLED: low
    DRIVE had shrunk all network amplitudes to sub-mV — doubled
    rg_to_pf/pf_to_mn/descend conductances + slow ADAP τ=1.9 s for
    rhythm), subtalar/mtp ligament-surrogate springs (OpenSim coordinate
    stops lost in conversion; foot flopped to 128°), follower armature
    1.0; (3) afferented ground + interleg at --drive 2.5: 21 s up, COM y
    ±2–3 cm, adduction ±8° (stance-gated BAL_LAT abductors working,
    mujoco −y = opensim +z). PELVIS-LIMBO DIAGNOSED: pelvis_tilt 36°
    despite 400 N·m/rad rig + new IMU trunk controller (BAL_TRK_EXT/FLX,
    torso up-vector PD → ercspn/obliques) — the IMU levels the TORSO
    (lumbar counter-tilts) but saturated hip extensors (glut_max 0.94,
    semimem 1.0) pitch the planted-leg pelvis backward; only balanced
    activation patterns fix it. **LIMBO FIXED at the source 2026-09-11
    night — the IK/NNLS chain RAN END TO END (details + traps in
    spinal\DESIGN.md night section):** `bsolve_ik.py` validates the
    converted model vs OpenSim along subject01_walk1_ik.mot (lengths
    median r 0.89, moment arms sign-exact after folding in MuJoCo's
    transmission minus) and back-solves per-timestep activations
    (SO-style ridge + [0,1] bounds via lsq_linear — plain NNLS explodes;
    measured GRF applied at the CoP; 6 Hz filtering like SO);
    `fit_pf.py` refits W_PF_MN/W_POSTURE from the back-solved group
    profiles → `fitted_walk_params.json`; hip_ext stance drive stack
    1.17 → 0.27 vs human peak 0.34. `optuna_walk.py` v3 (study
    ground_walk_v3; winner in best_walk_params.json with pf_gain) adds a
    global pf_gain dimension — back-solved weights are ~10× SMALLER than
    the old hand-tuned table the gain structure was tuned against (at
    gain 1 the sim barely moves) — plus a duty-0.6 reward (runner --eval
    reports duty). v3 winner score 1.521 vs v2's 1.091. FULL 22 s ground
    walk at rig 1.0 with `--fitted --best`: knee −22° REAL flexion, hip
    amp 38°, tilt 29°, stays up. Weaning ladder (`wean_rig.py`;
    `--rig-scale S` = rig stiffness×S, damping×√S): S=0.8 completes and
    keeps stepping (knee −23°, hip 37°) but leans to 37° — the support
    boundary is S≈0.8–1.0 until the pelvis-balance piece exists.
    **THREE SILENT TRAPS (probes diag_force/diag_frontal/diag_knee/
    diag_trans.py):** (1) pelvis slides load IDENTITY — the converter
    preserved OpenSim's coordinate values (keyframe
    qpos[pelvis_ty]=0.95 ↔ pelvis world z 0.95); the z-up remap lives in
    body frames; GRF VECTORS still remap os(x,y,z)→mj(x,−z,y). (2)
    Equalities MUST be disabled for INVERSE dynamics — at deep knee
    flexion the polyfit followers generate ~890 N·m of spurious
    knee-row wrench (forward sim still needs them). (3) data.act drives
    muscle force under mj_forward; data.ctrl is inert there
    (dyntype=muscle; ctrl only feeds act dynamics in mj_step).
    **Ben's night addendum (verified, see spinal\DESIGN.md):**
    MuJoCo's actuator_moment does NOT propagate through the eq couplers
    — at the knee use bsolve_ik.py's fd_moments (central-difference
    dl/dθ with followers re-projected = OpenSim's definition). ALL 36
    knee couplers are the vastii pathpoints — vastii→tibia actuation
    verified (mean arm −0.047 m vs OpenSim −0.045 m; the small AC
    anti-correlation = coupler-polyfit artifact, same family as the
    peronei/subtalar SO mismatches). "78 muscles" = 92 minus 14 with
    Fmax ≤ 5 N: the 6 intentional prunes PLUS 8 the converter shipped
    at 1 N — ercspn/intobl/extobl/ext_hal r/l — so the IMU trunk
    controller currently drives 1-newton muscles; fix = set those Fmax
    from stock gait2392 in patch_xml, then refit + v4 (needs Ben's go,
    invalidates current v3 tuning). plot_run.py load_benchmark now
    honors the .mot inDegrees header (subject01_walk1_ik.mot is
    DEGREES — the np.degrees() there was Ben's suspected double
    conversion; fig4 legend renamed "SNS sim mean"); ground figs +
    ground_walk.gif refreshed in Dissertation\CPG_airstepping_figs from
    the post-v3 S=1.0 run.
    **v4 KINEMATICS CAMPAIGN (2026-09-12, Ben: "fine-tune until the
    kinematics are similar to OpenSim"):** `kine_ref.py` is the
    acceptance metric — reference cycle from subject01_walk1_ik.mot
    phased by GRF onsets (1.23 s / 0.81 Hz, duty 0.61, knee −69.7°, hip
    43°, ankle 23°); runner --eval now reports `kine` + `kine_score`
    (cycle-normalized hip/knee/ankle shape RMSE + peak-knee + range +
    duty; 0 = perfect). optuna_walk v4b (study ground_walk_v4b_kine; a
    first v4 died to a broken −25 no-rhythm sentinel plateau —
    no-rhythm must score ≈ −65, what a frozen model really costs) ran
    60 trials with 2 new knobs (`desc_f` = DRIVE→RG-F,
    `e2_adapt` = PF_SHAPE["E2"] adapt; runner --best applies both):
    best kine_score −61.2 vs baseline −63.6 — converged. HONEST
    DIAGNOSIS (details + priority list in spinal\DESIGN.md 2026-09-12):
    the remaining gap is ARCHITECTURAL, not scalar — (1) E-duty 0.27 vs
    0.61 (half-center+adaptation tops out ~0.3 → needs sensory phase
    reset into the RG, which also fixes cycle-to-cycle phase jitter),
    (2) swing knee still extension-dominant (needs phase-specific quad
    SUPPRESSION, not more flexor drive), (3) ankle PF-dominant
    (POSTURE_OVERRIDE soleus/tib_post tone rides into gait), (4)
    cadence 1.18 vs 0.81 Hz (rg_adapt at range edge). Current best
    config = the v4b winner (runner --fitted --best reproduces); the v3
    winner is still better on the stability-shaped objective (its study
    remains in the db). **NEXT: ground duty 0.6 (E-duty still 0.16–0.34
    vs human 0.6), ankle balance, the pelvis-balance piece to wean
    below S=0.8, then vestibular/ocular (Ben) and cerebellum/BG
    layers.** New
    tools: `diag_stab.py`, `diag_phase.py` (adaptive thresholds),
    `_muscle_direction_test.py`, `draw_circuit.py` (Rybak-style schematic
    PNG), `neuro_scope.py`; runner flags `--no-ground`, `--no-interleg`,
    `--leg-damping X`, `--time N`, `--view` (3D window, full speed),
    `--scope` (live neural traces), `--realtime`; per-0.5 s
    forensics + NaN dof dump; `diag_nan.py` superseded. Lit grounding:
    Rybak/McCrea RG+PF, Bunz 2026 (reflex
    speed control), Ben's Zotero "Sensory Afferent Database" collection;
    read Di Russo/Ijspeert/Bouri 2023 JNE before any novelty claims.
    **2026-09-12/13 overnight (v5/v6 + a CRITICAL pre-existing fix —
    morning report at TOP of spinal\DESIGN.md):** (1) **`--best` pf_gain
    bug FIXED**: it tested `json["params"]` for pf_gain (pf_gain lives at
    the JSON TOP level) so the gain-scaling branch NEVER ran — every
    pre-fix `--fitted --best` reproduction since v3 silently evaluated a
    pf_gain=1.0 config (the recorded v4b/v3 full-run numbers — E-duty
    0.27, knee −14..+26, −63.3 — are gain-1.0 configs, NOT the true
    study winners). Fix scales ONLY fitted-file entries (trunk keys are
    not in the fitted json; scaling them was also wrong but inert — 1-N
    trunk muscles). Regression gate now BIT-EXACT:
    `runner --fitted --best --eval --drive 2.5868` = −61.1741850610298
    on current code; draw_circuit composites fixed identically
    (circuit_* weights numbers changed slightly — regenerated + recopied
    2026-09-13). (2) **v5 sensory phase-reset**: HIP_EXT_SIG/HIP_FLEX_SIG
    ports → PRESET_E/F INs → RG (ext: E↑F↓ prolongs stance; flex:
    F↑E↓ triggers swing), params.G phase_reset_e/f + PHASE_RESET dict,
    defaults 0. **CONDITIONAL TOPOLOGY: new neurons are built only when
    their gain > 0** — zero-g synapses alone change BLAS summation order
    and chaos turns the 1e-16 drift into real kine shifts (byte-identity
    at 0 is the regression contract). v5 study ground_walk_v5_phase
    (150 trials, seeded v4b): best −60.987 (trial 137, gains 0.22/0.22)
    — TPE AVOIDS the gains; tonic stance-gated position/velocity input
    (≤1 nA vs ~4 nA DRIVE) is effectively inert; phase_reset_e ~0.6-0.75
    + weak drive LATCHES RG-E on (stability constraint). (3) **v6
    swing-knee quad suppression — THE LEVER THAT WORKS**: PF_F1 → KINH
    inhibitory IN → knee_ext MNs (phase-gated by F1 = swing-only;
    params.G f1_kneext_inh). Study ground_walk_v6_kneext (110 trials,
    seeded v5): best −59.426 (trial 106, f1_kneext_inh 0.59), full-22s
    −60.59, duty 0.47, knee_min still +10, cadence rose to 1.4 Hz —
    next bottleneck = RG duty/cadence architecture. Reproduce: `runner
    --fitted --best6` (v5: --best5; old v4b: --best) — v5/v6 evals use
    repr(drive) so their jsons reproduce bit-exactly (4-dp rounding
    costs only ~0.005). Audits _phase_reset_audit.py / _kinh_audit.py
    PASS; hand sweep v5_sweep.csv; per-trial CSVs v5_results.csv /
    v6_results.csv; global-best full runs v5_best_trial*.npz (9) /
    v6_best_trial*.npz (10); runner flags --phase-reset E F,
    --kneext-inh X, --best5/--best6; RUNNER_DUMP_STATE=<file> env dumps
    all mutable params after arg-parsing for config diffing. (4)
    **Phase 1 deliverable**: `draw_circuit.py --vclasses` renders
    circuit_full with a Shevtsova-2026-eLife-RP107480 / Rybak-2015
    V-class strip + tags (our circuit: V0D/V0c+dI6 analog = commissural
    inhibition, V1/V2b = half-center+PF+Ia inhibition, V2a = RG→PF/PF→MN
    excitation, **NO V3 analog** — all cross-side connections
    inhibitory) → circuit_full_vclasses.{pdf,svg,png} in
    Dissertation\CPG_airstepping_figs; full class-by-class mapping table
    in spinal\DESIGN.md. **2026-09-13 (day, Ben's go received): trunk
    Fmax fix APPLIED** — `_fmax_audit.py` (every actuator gainprm[2] vs
    stock osim; the 84 active muscles match within the converter's
    RoM normalization ≤15%; ONLY ercspn/intobl/extobl/ext_hal r/l were
    1 N vs stock 2500/900/900/162) + patch_xml repair 2c (build-asserted
    8 tags). Post-fix the v6 winner IMPROVES (eval −59.241, tilt 27.4 —
    BAL_TRK finally drives real muscles); no retune needed; old
    --best/--best5 recorded numbers remain pre-fix values. Torque
    budget (`_torque_budget.py`, fd arms × stock Fmax vs bsolve ID
    demand): hip 3.4x, knee 2.8x, trunk 2.1x, ankle 1.4x (tightest;
    ankle demand 2.7 N·m/kg runs hot vs ~1.5 literature — known
    subtalar/CoP residuals). MuJoCo-muscle fidelity note (Ben Q):
    converted actuators = rigid-tendon simplified Thelen (activation +
    F-L-V, NO series elasticity); options = tendon stiffness
    (approximate) or custom Millard/Thelen plugin (full). v6 render
    ground_walk_v6_trial106.gif + v6_trial106_fig1..6.png in
    Dissertation\CPG_airstepping_figs (spinal_run.npz now holds the v6
    winner run). Dissertation-ready draft section:
    Dissertation\CPG_spinal_section_draft.tex. **2026-09-13 (Ben's
    circuit-figure critique): `draw_circuit.py --which deng` →
    circuit_dengstyle.\*** (figures/ + Dissertation folder) — the
    Deng/Nourse-Fig-6A-style layered schematic, STRUCTURE-DRIVEN: it
    builds the representative network, reads every connection from the
    compiled SNS object (`_net_edges.py` is the standalone inventory
    dump) and ASSERTS each (src,dst,sign) edge group is drawn — the
    figure cannot drift from the code. Deng Table A6 quoted in the
    footer. Honest accounting built in: Renshaw NOT implemented
    (ghost), RG mutual inhibition + Ia reciprocal are direct/lumped
    (Deng's Ext/Flx-IN and IaIN drawn as ghosts), Ib = autogenic-inh +
    stance-gated IBEXC reversal (Deng's IbIN→MN is EXCITATORY —
    opposite sign), ours adds II afferents, monosynaptic Ia→MN,
    commissurals, PRESET/KINH INs, and the real conversion-map
    formulas (a=clip(V/5mV,0,1); (L−Lmid)/Lhalf, L̇/0.6, F/Fmax →
    afferent currents). ADAP cells = burst termination, NOT fatigue.
    Full Deng/Nourse connection table in spinal\DESIGN.md 2026-09-13
    section. **2026-09-14: JOINT RoM LIMITS ON** (patch_xml repair 2f2:
    all 14 driver hinges limited="true" with stock ranges — knee
    [-120,+10] deg enforced; hyperextension artifact GONE, ground knee
    now -0.3..+11.4). **RENSHAW CELLS IN** (per-pool RC, MN→RC 1.0
    exc / RC→MN inh G["renshaw"] / RC↔RC once per pair, default 0,
    --renshaw X, ran 0.5). Post-limits+Renshaw: ground stays up, knee
    -0.3..+11.4, tilt 30.5, but gait re-timed (E-duty 0.15, 1.32 Hz —
    v7 retune needed); **AIR-STEPPING now 0.41 Hz, knee -97..+10, hip
    -21..+50** (Ivanenko regime; v6_rom_{ground,air}.npz +
    hindlimb_style_{ground,air}.png/pdf + {ground,air}_rom.gif in the
    Dissertation folder). **Toolbox-native diagram**: spinal_layers.py
    (Tutorial-4 Network subclasses RG/PF/motor) + _render_diagram.py →
    sns_diagram_spinal.png via the OFFICIAL sns_toolbox.renderer
    (graphviz BINARY installed into the myo env — <env>\Library\bin,
    NOT on PATH unless the env is activated). **Tutorials 1/2/4/8
    executed** (all 9 notebooks in spinal\sns_tutorials); **NEW SKILL:
    C:\Users\Ben Bolen\.zcode\skills\sns-toolbox\SKILL.md** (junction
    into ZCode_Skills — Ben commits; env, 1.5.2 API traps:
    add_population needs shape=[1], backends take net.compile() params
    not the Network, two class trees, connections-dict schema,
    add_network flattens; net-sizing: numpy backend fine at 8k neurons
    tested, practical wall ~25-30k dense-matrix memory — our
    410-neuron net nowhere near limits).
    **OPTIMIZER-RELEVANT INSIGHTS (2026-09-13, Simulink-realization
    session — read before interpreting v7+ trials):** (1) The network IS
    a self-sustaining limit cycle at CONSTANT DRIVE in numpy — 20 s clean
    alternation at bare DRIVE=2.5, period 1.21 s, amplitude constant
    (`spinal\check_selfsustain.py`). Tonic collapse in an eval means the
    CANDIDATE PARAMS killed the cycle (adaptation/drive balance), not
    that constant-drive testing is invalid. (2) **The network sits near a
    bifurcation**: an EXACT Simulink port (wiring verified to 4.2e-6 mV
    at t=0.3 s, same 2 ms Euler) settles to the TONIC fixed point where
    numpy oscillates — the reggate_v5_0 summation-order chaos tips
    marginal limit cycles either way. Winners that oscillate with only a
    small basin (drive/adaptation at the edge of the oscillatory range)
    will NOT transfer to Simulink/co-sim or survive BLAS-order changes;
    consider a basin-robustness gate on finalists: perturb params ±1 %
    (or shuffle a summation order) and require the cycle to persist 20 s;
    report the basin margin alongside the score. (3) Pointwise agreement
    between ANY two integrators dies at ~0.4 s (1e-9 @ 0.1 s → 1e-6 @
    0.3 s → O(1) @ 0.5 s) — eval-vs-eval differences late in a cycle are
    PHASE, not signal; only seeded identical binaries reproduce
    bit-exactly (existing contract). (4) Deng-style persistent-Na
    half-centers need tau_h FIXED — sns_toolbox's tau_h(V) quenches them
    (`spinal\deng_cpg_ode.py`: toolbox-tau flatlines, fixed-350 ms
    self-runs at 1.94 s; Simulink files `SNS_Simscape\demos\SNS_Deng_*.slx`
    reproduce). Same transfer risk if tuned candidates are ever ported to
    Animatlab/Simulink realizations. (5) **CORRECTION 2026-09-16 (Ben)**:
    sns_toolbox 1.5.2 DOES ship `NonSpikingNeuronWithPersistentSodiumChannel`
    (Tutorial 8, executed in `spinal\sns_tutorials` 2026-09-14) — earlier
    "toolbox can't express NaP / ADAP is the only burst-termination
    substitute" claims were WRONG. The spinal RG can be built as literal
    Deng persistent-Na HC neurons; FIRST test the real class's tau_h(V)
    with tau_max_h at depolarized V (the quenching result in (4) was
    measured on the hand-coded formula in deng_cpg_ode.py, not the class).
    **(6) Basin gate IMPLEMENTED (2026-09-16, laptop): `spinal\basin_gate.py`**
    (myoconv env) operationalizes this: loads a RUNNER_DUMP_STATE json
    (`state_v10_best.json` / `state_v10_study.json`), perturbs every scalar
    param ±1% (seeded), runs the network-only constant-DRIVE rhythm 20 s (as
    check_selfsustain.py), PASS = last-5 s swing ≥ max(0.1 mV, ½ baseline
    first window). RESULT: **v10-best (drive 2.929) is basin-robust — 12/12
    perturbations persist, final swing 0.98× baseline, ~6.5 mV, period ~2 s
    (5 peaks/10 s — don't impose a fixed peak floor, slow winners mislabel);
    `state_v10_study` does NOT self-sustain at constant drive (baseline decays
    to 0 — it depends on the runner input schedule, consistent with E2's
    tonic collapse)**. Caveat logged: an identical-seed trial EXPLODED (1e164)
    in a 2 s smoke run but was stable in the 20 s run — multithreaded-BLAS
    summation-order nondeterminism is real even between processes; treat
    single borderline trials as noise. Results: `spinal\basin_gate_results.json`
    + `basin_gate_full.log`.
    **2026-09-20 (EB475WS4) — CURRICULUM FIXED + FIRST GROUND GAIT; full
    chain in spinal\DESIGN.md "2026-09-20" section:** the 09-18 stage-2
    "winner" was a STATIC-POSE EXPLOIT (rises=0; score 36.906 = knee term
    alone; bit-exact repro) — the air objective now RHYTHM-GATES (rises<3
    → ≈−10), new gains enter at [0,0.5] (full-range sampling killed the
    stage-1 rhythm in all 25 trials), seeds are full-dicts (missing keys
    get SAMPLED — JSON-rule lesson), exploited studies archived
    (curriculum_exploit_archive_20260920.json) and rerun as curr_s2b/s3b.
    NEW KNOB `contact_onset` (params.G, runner-side edge-triggered heel/
    toe transients, default 0 = bit-identical, JSON-RULE loader updated).
    OUTCOMES: s2b 50 trials — seed (stage-1 winner, drive 3.15) still
    best 129.05; best afferented 128.3 ties it with all new gains ≈0
    (afferents tolerated, not additive). s3b 60 trials — 11/60 real
    gaits, winner trial 52 = −85.22: 20 cycles, duty 0.69 (ref 0.61),
    knee −53.5° cycle / −76..+10 raw, hip range 16.9° (43.3 ref), ankle
    34° OVER ref, tilt 29.7°, no NaN. Scalar search PLATEAUED — remaining
    duty/cadence/hip-range gaps are ARCHITECTURAL (sensory phase reset
    into RG). Also: runner npz q = named KEY_JOINTS already in degrees
    (the old "q[m,4] is a quat" bug is resolved; _diag_stage3's double
    np.degrees fixed); tuned-config constant-DRIVE self-sustain is LOST
    (schedule-sustained; `_debug_rhythm_tuned.py`). Run instructions for
    Ben's Spyder: `Code\MuJoCo_SNS\HowToRunCode.md`. Dissertation
    CPG_spinal_section_draft.tex \fillme slots FILLED (ZCODE 2026-09-20
    fenced blocks) + 3 new figs curr3b_ground_* — deconfliction note for
    the other dissertation chat is in CHATGPT_HANDOFF.md "DISSERTATION
    CPG SECTION — DECONFLICTION".
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
    **`minimizeFlxPin.m` is the evaluator of record (2026-09-11)**: Ben's original
    evaluator re-worked in place with the full 2brk mechanics. **Convention since
    2026-09-13: MIXED — insertion/tibia bracket ONE-TRANSFORM (pitch-only)
    `Tkbr = RpToTrans(RkbrZ, Pbri')` with K = [X1,X2,X1] and the 1-rotation
    pbrBnew reconstruction; origin/hip bracket TWO-ROTATION with
    K2 = [X1,X1,X2]** (supersedes the 2026-09-12 isotropic-K2 test and the
    all-2trans port; the y-symmetry proof makes the tibia 1trans switch
    numerically neutral — it is Ben's structural preference). Pbri =
    [-48.11,-107.81,13.8]/1000, Pbr2 = [-52.61, 0, 75.06]/1000, `persistent kf`
    data cache, +5.3° encoder offset on kf(3) (47cm; Ben moved it back from
    kf(4) on 2026-09-12 — see data-concern note below).
    **Test-5 pB shift (2026-09-13): env `FLX_T5_YMM`** = millimeters of permanent
    +y offset in the INSERTION-BRACKET frame for the 42 cm BPA only
    (`pbrBnew = [norm(pkbrB(1:2)), t5y, pkbrB(3)] + eB`). A pure knee-frame +5 mm
    y is UNREACHABLE by thetabrB rotation (point B sits nearly straight above the
    bracket, θ=91.8°, asin-arg 1.079); the bracket-y offset is equivalent to
    thetabrB +4.56° ("a few degrees more") and models a permanent ~5 mm bend of
    the bracket arm along the bending axis. Unset/0 = off. Commented 1trans
    alternatives for the Tkbr/Thbr frames AND the pbrBnew/pbrAnew reconstructions are
    kept in Lok (restored 2026-09-12 after Ben's pitch-only toggle experiment hit the
    missing lines) — toggle frame AND reconstruction together. Signature unchanged
    `(Xi0,Xi1,Xi2,idx_val)`; escape paths kept (X1=X2=kSpr=Inf → rigid; X1=X2=Inf +
    tendon → cable-only). Port verified vs minimizeFlxPin2brk 2trans: M_p identical on
    all 5 tests (finite Xi and escape), GoF identical on tests 1/2/5, differs only on
    3/4 via the offset move (smoke_FlxPin2brk_port_20260911.m in Dig_out).
    `minimizeFlxPin2brk.m` (the previous evaluator; 6th `transMode` arg / driver's
    `FLX2BRK_TRANS` env selects 1trans|2trans — the laptop session's separate-file 1trans
    evaluator was DELETED 2026-09-08 at Ben's ruling: transMode is the vehicle of record)
    is now STALE: pre-mixed-convention (2trans both brackets, pre-isotropic-K2
    [X1,X1,X2] 2trans / [X2,X1,X2] 1trans) — do not re-run it as-is. Mechanically it matches the
    port: insertion bracket two-rotation frame K = [X1,X2,X1] (Sept-8-morning arm-1
    runs used [X1,X2,X2]); origin bracket at `Pbr2`: 1trans pitch-only frame →
    K2 = [X2,X1,X2]; 2trans two-rotation frame → K2 = [X1,X1,X2] (Ben, late 2026-09-08).
    `USE_BRACKET2` flag).
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
  abduction DOF, LH_HipZ→LH_Hip rename pending. **LATCH PRIME SUSPECT (Simulink
  session, 2026-09-13, verified in `spinal\_tau_h_check.py`): the Na h-gate time
  constant. sns_toolbox's tau_h(V) formula collapses to ~0.1 ms at depolarized V
  (removes burst termination, quenches the oscillator; fixed τh=350 ms self-runs) —
  when opening the .aproj, FIRST check how LinearHill implements the Na h-gate tau
  (Animatlab treats tau_h.max as a fixed constant, which is why the Simulink Deng
  port works). Related trap: SNS_Library vs sns_toolbox use OPPOSITE synapse
    saturation conventions (ThrPre/Elo) — keep straight when porting values.
    **2026-09-16 update: sns_toolbox DOES ship
    `NonSpikingNeuronWithPersistentSodiumChannel` (Tutorial 8) — "toolbox
    can't express NaP" was wrong; if the MuJoCo-side RG moves to that
    class, verify its tau_max_h semantics against the quenching result
    before assuming cross-platform consistency (see
    SNS_Simscape\README_SNS_Simscape.md correction note).****
    **2026-09-14b: START POSE = normal.mot** (Ben supplied the gait2392
    "normal" Coordinates values — pelvis_tilt −1.87, knees −3.9/−8.2,
    hips R+24.6/L−16.6, ankle_l +9.8 dorsi, lumbar set, ty 0.96;
    `START_POSE_DEG` in runner.py, applied AS GIVEN — an auto-sign-flip
    heuristic was REMOVED after it mangled the canonical values; applied
    before standing solve + rig so springs/solve hold it; 21-22 muscles
    engage). **v8 (hand pose) BREAKTHROUGH: −47.9 eval in 60 trials**
    (v7 plateau was −62.5) — pose quality is a first-order lever; winner
    used quad suppression 1.0 + phase_reset_f 1.72. v8b on normal.mot
    FIRST LAUNCH COLLAPSED onto the −65 frozen sentinel (real walkers
    score −76..−86 on the pose) — sentinels recalibrated (−100 frozen /
    −110 NaN), eval schedule lengthened to 16 s (11 s walk window) so
    0.3–0.9 Hz gaits yield ≥3 countable cycles, seeded with the v8
    winner; fresh study `ground_walk_v8b_normal` (optuna_walk_v8b.py,
    --best8b). RULE: when plant/pose/objective change, sentinels must
    sit BELOW the worst genuine walker. **v9 (flat-foot height +
    amplitude objective): −65.4/60 trials → −62.24/120 trials**
    (reproduced bit-exact); **amplitude now MATCHES OpenSim** (hip
    46.4/43.3, knee 66.9/70.5, ankle 24.1/23.1, knee_min −78/−69.7);
    remaining = hip phase (inverted vs RG anchor), ankle −45° PF
    OFFSET (posture tone), duty 0.17, cadence 0.3 Hz — regenerated
    overlay in the Dissertation folder. **2026-09-14 evening: ankle
    trim + TRANSIENT reset** — `ankle_post_walk_trim` (POST bias of
    ankle_pf group scales toward 0 with drive; 1.0 = v9-identical);
    PRESET_E/F now have a FAST adaptation loop (PREA τ0.08, gain 1.5)
    = high-pass ONSET detector (tonic ≤1 nA proven inert; rectifying
    synapses pass only the onset pulse); **v10 study
    `ground_walk_v10_transient` running** (seeded v9, +trim searched,
    --best10). **ENV INCIDENT: the conda graphviz install clobbered
    myo\python.exe (repaired via --force-reinstall python=3.10.21;
    pip pins survived) — verify python.exe after any conda
    transaction in this env.** **v10 (transient reset + trim): −65.375
    (trial 56, 60 trials), reproduced bit-exact via `runner --fitted
    --best10`** after TWO catches: the --best10 flag branch was lost in
    successive same-anchor edits (chain silently ended at --best9 —
    state-dump diff + missing loader print caught it), and the trim
    loader branch was missing (JSON RULE again). Winner = knee −75,
    tilt 15.7, amplitudes hold (45/62/25 vs 43/70/23), trim winner
    0.097 (near-zero standing PF tone wanted in gait — confirms the
    set-point diagnosis). Duty/cadence/hip-phase remain
    architecture-level; study resumable. **2026-09-15: CURRICULUM +
    mechanosensory stance feedback + IaIN (Ben's methodology question
    became the plan)** — (a) heel/toe contact mechanosensors (per-foot
    MuJoCo contact normal forces → HEEL_c/TOE_c ports → HEEL_IN→RG-E
    S2W trigger + TOE_IN→RG-E late-stance prolongation; gains
    heel_rge/toe_rge default 0); (b) LBIN stance-Ib group IN (per
    Dominguez 2020: INs in the rhythm-generating layer) → RG-E
    prolonger (gain ib_rge default 0); (c) IaIN population replacing
    direct Ia→antagonist when G["ia_in"]>0 (PF_F1 phase gate +
    RC→IaIN disinhibition per Hultborn); (d) mutual Renshaw fix (was
    one-directional). ALL conditional topology, defaults 0 =
    v10-identical at renshaw=0. **Regression gate now measures the
    DELIBERATE mutual-RC change**: v10 winner re-scores −79.9 (was
    −65.4 on lopsided wiring) — retune needed, not a bug. Staged
    curriculum (_curriculum.py): s1 air-deaff (25 trials, best 104.6
    objective = 3×rises + swing depth), s2 air-aff (+heel/toe, 25
    trials, kine −68.1), s3 ground (+ib_rge/ia_in/trim, 30 trials,
    kine −72.5). All new pathways tuned nonzero (heel 0.91, toe 0.63,
    ib 0.68, ia_in 0.62, trim 0.63). Final v11 run: stayed up, knee
    −78.8..+12.7 INSIDE RoM (deep swing flexion), hip −24..+58, tilt
    −2..+19, ankle −90..−0.8. Gaps remain: duty 0.14, cadence 0.46 Hz,
    ankle PF bias. START POSE = normal.mot values (see
    START_POSE_DEG in runner.py). Zotero Web API access (AARL group
    735051) verified for full-text PDF retrieval; credentials in
    D:\Github\api_credentials_local.txt (outside repo). Musculoskeletal
    audit script: _muscle_force_compare.py (R² 0.87–0.98 matched
    activation). Literature audit report:
    LIT_CIRCUIT_AUDIT.md (P1a/P1b IMPLEMENTED; P2a/P2b/P3 open).
    **2026-09-14 (v7 retune under corrected physics): ANKLE DORSIFLEXION
    MECHANISM FOUND** — boosting swing DF drive does nothing (9 PF
    muscles ~10 kN vs 3 DF ~1.6 kN); `params.G["f1_anklepf_inh"]`
    (KINH→ankle_pf MNs, swing-gated, same pattern as knee) at 1.0 gives
    air-stepping ankle **−55..+7.1° true dorsiflexion** (tests
    `_ankle_test.py`/`_ankle_anh_test.py`). **v7 study**
    `ground_walk_v7_rom` (140 trials, Renshaw 0.5 fixed, anklepf gain
    searched): best **−62.523 (trial 100)**, plateaued; reproduces
    bit-exact via `runner --fitted --best7`. Winner = slow stiff shuffle
    (duty 0.18, knee pinned AT the +10 cap, 0.4 Hz, tilt 28.6) — scalar
    tuning CONVERGED; duty/cadence/flexion gaps are architecture-level
    (transient reset / FSA-analytic seeding next). **JSON RULE**: every
    new params.G knob a study uses MUST be saved in the json params AND
    have an `if key in best` branch in the --best loader — v7's renshaw
    omission cost a silent 0.14 kine delta, caught only by the
    reproduction check.
  **Toolchain + format traps:**
  `Biped_2xCPG_wSubs\tools\` (build_A2→build_B2→repair_dup_ids2→fix_gas_ia3→
  rebuild_pages_v3→move_links_by_page→set_fmax2; pipeline from pristine git 6437441
  backup — re-run the WHOLE chain, never patch a patched file). Traps: page CDATAs are
  NOT covered by whole-file XML validation (validate each page separately); tag
  depth-counters must skip CDATA and match `<Node attr>` opens; page rebuilds must emit
  `<Version>` + closing tags; cloned blocks need child-object GUIDs re-rolled; **NEVER
  hand Ben an .aproj without opening it in AnimatLab2.exe first** (zero Error dialogs =
  pass; each broken page throws one dialog) and never save from the GUI (taskkill /F).
  Full brief: `CHATGPT_HANDOFF.md` AnimatLab section. **RG oscillation root cause
  (2026-09-15; full map in `Biped_2xCPG_wSubs\tools\RG_contact_drive_notes_20260915.md`):**
  the working reference is `origin/AddingStepSensor_CoMorrow_stw:...Walker_2_Layer_CPG/
  Walker_2_Layer_CPG.aproj` — its RG is CONTACT-driven (spiking foot-contact neurons →
  own RG ext G=6 + contralateral RG flx G=6, SpikingChemical Equil −50; NO tonic anywhere;
  NonSpiking RG InitialThreshold +50 mV; connexion G 0.5 where ours is 5). Proven on our
  standalone .asim: a tonic-driven half-center CANNOT oscillate at any drive/G/threshold
  (latches; the Ca plateau holds the ON cell). The "-55 mV RG fix" was a class mix-up
  (−55 belongs to W2L's SPIKING contact neurons) — reverted same day. **The standalone
  .asim's Root (pelvis) ships `<Freeze>True</Freeze>`** (31 kg box at y=1.02) → the biped
  hangs, feet never touch ground, ContactCount=0 on all 41 bodies; Freeze→False + generic
  PhysicalToNode adapters (`SourceDataType=ContactCount`, source=toe body GUID,
  TargetDataType=ExternalCurrent, Gain C = amps/contact) makes contact drive WORK headless.
  Harness builder: `%TEMP%\build_contact.pl <gainC> <ve> <vf>`. AnimatLab traps: chart
  DataColumns need valid GUIDs (invalid → that chart silently writes 0 bytes);
  SimEndTime == chart EndTime kills the end-of-sim flush (keep sim end longer); an
  adapter's Gain child object needs its own fresh GUID. Next: standing tone →
  contact-asymmetry stepping → G rescale toward W2L values → port via tools pipeline.
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
  **2026-09-16: ACTIVE on the SECOND route too — Simscape Multibody Link add-in v7.4
  export of `09_BA_003.SLDASM` → `Knee assembly\09_BA_003.xml` (+`09_BA_003_error.txt`).
  Exporter DROPPED the 4 Hinge mates (Ben: make them revolute — replace each Hinge with
  Concentric+Coincident in SW, re-export) and the 6 BPA1/BPA2 mates to the assembly root
  (Ben: ignore; BPAs will be built NATIVELY in Multibody — expanding-diameter force law
  + N-segment sleeve with Spatial Contact Force vs the STL meshes — not imported from
  CAD); PathMate drops = patella, ignore. CORRECTION (Ben): the KB_R_003↔TI_R_006
  Concentric set is the tibial head BOLTED to the shank = rigid group, NOT the knee
  joint. Import inventory (which joints smimport made) still to run.
  Full brief for any assistant: `CHATGPT_HANDOFF.md` "ACTIVE WORK D".**
  sw2urdf v1.6.1 INSTALLED on the laptop (Ben, 2026-09-10; official build targets
  SW2021, so on SW2025 watch for the vanishing-dialog issue, issue #147); fallback:
  Simscape Multibody Link IS installed+registered (disabled — enable in SW Tools >
  Add-Ins). `import_simscape_when_ready.m` takes URDF or XML. Tendon parts
  (`Tendon_Extensor/Flexor.SLDPRT`) are built but NOT yet inserted into `09_BA_003`;
  GUI steps in the README. Gotchas: mask icon drawing commands take numbers only
  (no LineSpec/name-value; `color('black')` not RGB), and keep masked blocks square
  so circle icons stay round. **R2025a copies for the other machines (2026-09-10):**
  `export_slx_to_R2025a.m` wrote `{SNS_Library,KneeReflexDemo,sns_simscape_EI_test}_R2025a.slx`
  in this folder. The demo's 20 library links were re-pointed to `SNS_Library_R2025a`
  (exported slx files get renamed to match their filename), so those two `_R2025a`
  files must stay together in the same folder. Verified block-identical (67/67 and
  54/54 blocks, same type inventories) and 20/20 links resolve. The EI-test copy is
  best-effort only: Simscape-block downgrades are officially unsupported and it needs
  `simscape_sources\SNS_lib.slx` (not exported — if needed on R2025a, rebuild there
  from the `+SNS` sources via `sns_build_simscape_lib.m`). On R2025a machines, open
  the `_R2025a` copies, not the originals (R2025b format won't load).
  **MuJoCo↔Simulink bridge — INSTALLED & PROVEN (2026-09-12, same-day resume):**
  `mujoco_bridge\` = mathworks-robotics/mujoco-simulink-blockset on this machine,
  MJ_VER 2.3.7 (3.3.6 cannot parse `collision="predefined"`), mex-compiled with
  MinGW gcc 8.1 (README's 12.2+ requirement is NOT binding). All 3 prove-it
  tests PASS (`matlab\run_bridge_tests.m`, log `logs\test_abc.log`): fire-one-
  muscle (ctrl→force, knee extends 88°), two clock rates (0.001 s ctrl source →
  Rate Transition → 0.005 s plant), sensor readback BIT-EXACT vs Python myo-env
  ground truth. **One LOCAL PATCH to upstream `src\mj.cpp` is REQUIRED and
  applied**: upstream initData() zero-fills qpos (no mj_resetData) — our model
  starts 0.95 m underground and explodes; patch = `mj_resetDataKeyframe(m,d,0)`
  when nkey>0 (marked AARL LOCAL PATCH; rebuild via tools\setupBuild[MINghW]+build).
  Sensor-patched model `gait2392_simbody_cvt3_simbridge.xml` (cvt3.xml itself
  NEVER touched) generated by `matlab\make_simbridge_xml.py` — 2.3.7 sensor
  elements are `jointpos/jointvel/actuatorpos/actuatorvel/actuatorfrc` (the
  old report's "actuator_length" was wrong). Simulink gotchas banked in
  BRIDGE_REPORT.md (xmlFileRel needs an ABSOLUTE path; log the sensor bus via
  port DataLogging, not To Workspace; Rate Transition lives in Signal
  Attributes). Downloaded binaries are gitignored. Not yet wired into any
  pipeline — the SNS-side counterpart (Part 3) now exists, below.
  **LAPTOP PORT DONE (2026-09-16, R2025b):** blockset re-extracted from
  `mujoco_bridge\blockset.zip`, MJ_VER 2.3.7, AARL patch re-applied (correct
  signature is `mj_resetData(m, d)` — two args), MinGW 8.1 support package
  already configured, all 4 mexw64 built; tests a/b/c PASS **bit-exact** and
  E1/E2 reproduce (logs `logs\*_laptop.log`). E2 re-saved
  `results\SNS_SpinalNetwork.slx` in R2025b format — R2025a machines use the
  kept `SNS_SpinalNetwork.slx.r2025a` backup or regenerate. New Simulink/
  smimport gotchas (find_system option ORDER matters; Mechanism Configuration
  only reachable by direct path `sm_lib/Utilities/Mechanism Configuration`;
  gravity param is `GravityVector`; a Solver Configuration block IS required
  for Multibody nets; smimport's return value is not the model name — import
  a sanitized temp copy) all banked in BRIDGE_REPORT.md "LAPTOP PORT".
  **2026-09-20 SNS_Simscape session (laptop, all details in
  README_SNS_Simscape.md):** (1) ALL FOUR demos verified on R2025b —
  KneeReflex 42.72°, CPG 9.7-48°, Beer 5.13°/9.66°, Deng 1.938 s r −0.857;
  fixed sns_run_deng_demo.m corr() (no Statistics Toolbox on laptop →
  base-MATLAB Pearson). (2) **SNS_Library.slx REBUILT with Rybak-style
  synapse icons** (pass-through axon + E-triangle/I-ball terminal at the
  postsynaptic edge; restyle in sns_build_library.m>synapseIconCode) — now
  R2025b-native, demos' links still resolve, SNS_Library_R2025a.slx
  exported; R2025a machines regenerate via the two build scripts. Actuator
  validation re-PASSED (4.9e-10 N). (3) **Dissertation exports**:
  sns_export_diagram.m upgraded — per-block 12 pt (R2025b has NO model-level
  FontSize param), synapse names hidden, white canvas, styling in-memory
  only; EPS = MiKTeX pdfcrop + mgs.exe eps2write chain (print -depsc is
  REFUSED for Simulink systems; R2025b bundles NO ghostscript; pdfcrop
  cannot write to C:\ root). NEW `demos\sns_build_circuit_view.m` →
  `demos\KneeReflexCircuit.slx` = print-only neural-circuit view (blocks
  COPIED from the demo, mask values 1:1). snsfig/sns_draw_circuit fonts
  9.5-11 pt on a 16.5 cm canvas (= text width, no LaTeX downscale).
  (4) **OpenSim→Simscape COMPLETE on the laptop**: Gait2392_simbody_simscape.slx
  (osim_import\) imports (1446 blocks), COMPILES + sims after patching the 19
  File Solids to absolute STL paths — **mesh param is `ExtGeomFileName`**
  ('FileName' doesn't exist); `osim_import\Geometry` is a gitignored junction
  into the cvt3 Geometry; a junction created AFTER import does NOT fix an
  already-saved model. Joint inventory: 85 Prismatic (pathpoint slides) +
  20 Revolute + 44 Weld = 149. (5) **09_BA_003 Multibody-Link XML imported**
  (`mdl_knee_rig_xml_imported.slx`, runner dev\imports2_20260920.m): 96
  blocks, joints = 2 Cylindrical + 3 6-DOF, ZERO revolute — **knee DOF
  missing** (dropped Hinge mates → welds); Ben's SW re-export (Hinge →
  Concentric+Coincident) still required. The older
  mdl_knee_rig_import_tmp_imported.slx = 1-link sw2urdf stub. (6) CRASH
  HAZARD: a FAILED `SimulationCommand update` leaves Simscape's GUI tree
  poisoned — physmod_sm_gui_app_tree.dll access-violation at ANY later touch
  incl. process teardown; make the compile succeed or expect crash-at-exit
  (log survives). .gitignore gained `Code/Matlab/SNS_Simscape/osim_import/Geometry/`.
  **Part 3 — tuned spinal network as an EDITABLE Simulink model (2026-09-12,
  same session, all VERIFIED):** `spinal\export_network_json.py` (myo env)
  dumps the tuned `--fitted --best` network (v4b winner) to `spinal\
  spinal_net_export.json` (410 neurons / 1392 synapses / 376 inputs / 92
  MN→actuator outputs; the SNS Network object names every input port
  "Input" — real names come from SpinalNetwork.inputs order).
  `SNS_Simscape\sns_units_test_2n.m` proves the units mapping (PASS 6e-4 mV:
  C=tau µF → Cm=1000·tau nF, Gm=1 uS, Vrest=0, synapse e_lo/e_hi=(0,5 mV) →
  ThrPre=0/SlopePre=5, g µS, Esyn mV, currents nA). `sns_build_from_json.m`
  generates `results\SNS_SpinalNetwork.slx` (~2600 SNS_Library blocks; u
  [376×1] input currents nA, S [92×1] MN drives in MuJoCo ctrl order;
  REGENERATE after any retuning — values live in block masks).
  `sns_verify_from_json.m` = wiring proof: **PASS, 4.2e-6 mV max dev over
  ALL 410 neurons at t=0.3 s** vs the numpy 2 ms Euler reference
  (`export_verify_ref.py` → verify_ref.mat; scipy string lists load in
  MATLAB as padded char rows — strtrim). **The network is CHAOTIC: any two
  integrators agree only to ~0.4 s** (1e-9 @ 0.1 s → O(1) @ 0.5 s; the
  same "reggate_v5_0" summation-order effect), and the RG period itself
  shifts ~35% between 0.1 ms and 2 ms stepping — **the v4b tuning is a
  property of the 2 ms Euler semantics; reproduce it in Simulink with
  fixed-step ode1 @ 0.002 s** (bit-compatible with the runner's numpy
  stepping). Extra Simulink gotchas banked in README_SNS_Simscape.md
  (ExternalInput unreliable under -batch — inject a Constant into the
  demux; log buses via port DataLogging; subsystem ports are numeric
  blk/1, never inner port-block names).
  **Deng RG/PF layers as separate Simulink files (2026-09-13):** `SNS_Deng_Library.slx`
  (persistent-Na HCNeuron, Nourse 2023 Tables A4-A7, tau_h FIXED 350 ms) +
  `demos\SNS_Deng_RG.slx` + `SNS_Deng_PF.slx` (IN-laminated half-centers) +
  `SNS_Deng_CPGDemo.slx`/`sns_run_deng_demo.m`: ONE 10 nA/20 ms pulse ->
  continuous 1.938 s alternation (matches `spinal\deng_cpg_ode.py` 1.94 s).
  KEY: toolbox tau_h(V) QUENCHES this circuit; fixed tau_h bursts. Same
  suspicion for the Animatlab RG latch. See README_SNS_Simscape.md.
  **First closed-loop experiments (2026-09-13, scripts `mujoco_bridge\
  matlab\e{0,1,2}_*.m`, results in `exp_out\`, details in BRIDGE_REPORT.md):
  E0 — hand-built Simscape Multibody (primitives, no CAD) FAILS on this
  machine; the license blocks at block-ADD time, so Ben's cylinder-elbow
  idea must run on the laptop. E1 — first CLOSED-LOOP SNS↔MuJoCo sim
  (Ia/II/Ib from vas_med_r sensors → tuned synapses → MN → ctrl): reflex ON
  bounds the post-collapse transient (469 N vs 43 kN, knee ~80° straighter,
  S recruits 0.40→0.97). E2 — full 410-neuron net → all 92 muscles via a
  Model block, single-rate 2 ms (new `simbridge2.xml` = the runner's
  timestep/implicitfast rewrite), 5 s clean — but bare DRIVE=2.5 lands on a
  TONIC fixed point (no rhythm): the production rhythm needs the runner's
  full input schedule (POSTURE/BAL/afferent gating), not DRIVE alone.**
- Repo root: `CHATGPT_HANDOFF.md` (brief for other AI assistants when ZCode is unavailable)
  and `CHATGPT_REPORT.md` (their report back; 2026-09-08 edition covers Overleaf
  manuscript-status edits) — keep both current when work is handed off.
- `SADb_audit\` — Sensory Afferent Database reconciliation (Zotero personal / AARL group /
  Airtable "Sensory Feedback" base `appMQTnobUNRytIp7`) + the standing curation backlog.
  Spec + live state in `SADb_audit\README.md`, per-batch details in `curation_log.csv`
  (WORKFLOW rows = current state). **State 2026-09-16: Papers table = 943 records, 500 with
  empty Notes; the 383 rest-import campaign has batches 1–5 DONE + audit PASS (50/383;
  batch 6 = queue CSV rows 41–50) with Ben's 09-15 rulings applied; the Sept-15/16 "task5"
  auto-curation pass created ~390 more records and flagged 108 as insufficient
  (`task5_progress\task5_insufficient.json` — manual curation needed); PDFs attached to
  PDFs attached to 318/886 DOI-bearing records; THE definitive no-PDF hunt list =
  `author_fix\remaining_no_pdf.csv` (581 records, doi/record_id/title/hunt_status;
  rebuilt by `rebuild_no_pdf_list.py`; the Sept-15 intermediate lists were deleted
  2026-09-20); author
  normalization done (Primary Author = one surname + Secondary Authors multi-select);
  VOSviewer citation map built (`SADb_audit\vosviewer\`; full-corpus rebuilder =
  `vos_build2.py`). **2026-09-18: batch 6 DONE (8 curated + echo-verified, 2 Elsevier
  chapters logged no-text pending Ben's library pass; twins Fujiki 2018 / Ekeberg 2004
  Models + Côté 2018 / Prochazka and Ellaway 2012 Review) — 451/943 records have notes;
  task5 did NOT cover the rest-import queue (`reconcile_queue.py` proves it); next =
  batch 7 = rows 52–60 (pre-grounded in `batch6\ground_*.txt`). Airtable-independent
  stack: `export_corpus.py` → `export\sadb_export.{json,csv}` (943 records, Excel-ready),
  `app\build_app.py` → `app\sadb_app.html` single-file OFFLINE explorer (Table
  search/sort/filter, Pivot with drill-through, bubble map; Airtable GET 422s on a
  fields[] filter because two fields share the name "Models copy" — fetch full records).**
  SADb work is DELEGATED TO
  CHATGPT during GLM peak hours (Mon–Fri 23:00–03:00 Pacific) — its brief is the SADb
  section of `CHATGPT_HANDOFF.md`; keys live in `D:\Github\api_credentials_local.txt`
  (rotation to scoped keys pending; Zotero key is READ-ONLY per Ben).

## REPO SIZE REDUCTION (standing objective — Ben, 2026-09-15)

Ben wants the repository size reduced ACROSS THE BOARD (folder currently
~13.8 GB) and is wary of any operation that forces a re-clone. ANY AI
session working on this MUST present the safeguard checklist below and
get Ben's explicit go PER STEP before running anything destructive
(filter-repo, push --force, gc/prune, bulk deletions). Do not fold these
steps into unrelated work.

**Measured composition (2026-09-15, EB475WS4):**
- `SADb_audit\pdf_staging\` = **5.8 GB, UNTRACKED** (another session's
  PDF-staging area for its Zotero upload workflow). Disk-only: deleting
  or archiving it involves NO git operation — but it is that session's
  working data: confirm with Ben/that chat before removing. This is the
  single biggest disk item.
- `.git` = 6.17 GB: history carries large CAD binaries (Solid_Models
  SLDPRT/STL/STEP across many commits) + the ~500 MB npz commit
  a8746f9 (HEAD of branch `KneeTestSetup_BenBo_stw`, PUSHED — so it is
  on GitHub; branch is otherwise in sync with origin).
- Tracked worktree ≈ 1.8 GB (Solid_Models 1.06 GB dominates).

**Phased plan (nothing executed yet — Phase A/B safe, C/D need the
checklist + go):**
- Phase A (disk, no git): decide pdf_staging (archive vs delete) —
  frees 5.8 GB immediately once the owning session's upload is done.
- Phase B (prevention): extend .gitignore (npz, __pycache__,
  SADb_audit/pdf_staging/, large binaries by policy) — Ben commits.
- Phase C (history, branch-scoped): strip `*.npz` from
  KneeTestSetup_BenBo_stw (the npz commits are unique to this branch —
  verify with `git log --all -- '*.npz'` first), via amend of a8746f9
  (tip) if that covers all of it, else git-filter-repo on that ref
  only; then `push --force` that branch; then reflog expire + gc.
  Requires: GitHub Desktop git at
  `C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.5\resources\app\git\cmd\git.exe`,
  git-filter-repo (pip) for the multi-commit variant.
- Phase D (policy, biggest long-term lever): CAD-binary history
  (Solid_Models revisions, Unused_Parts duplicates) — decide LFS vs
  external storage vs pruning duplicate-geometry commits. Most
  invasive; separate explicit go.

**MANDATORY SAFEGUARDS (any session, any step of Phase C/D):**
1. Backup FIRST, outside the repo: full `git bundle create
   <path>\pre_rewrite.bundle --all` + a copy of the dirty working-tree
   files (list from `git status --short`; 117+ entries as of
   2026-09-15). Verify both exist before proceeding.
2. filter-repo/reset does HARD RESETS — uncommitted tracked-file work
   is destroyed unless backed up and restored afterward. Re-apply
   intentional working-tree deletions after restore.
3. Always PROMPT Ben with this checklist + the specific commands before
   executing; he is wary of re-clone-requiring operations (his words).
4. After any force push: EVERY other clone (laptop DESKTOP-5Q16KE9,
   easteregg2) must re-clone or `fetch + reset --hard` + expire reflogs,
   or it will re-push the old history. Tell Ben explicitly each time.
5. filter-repo removes the `origin` remote — re-add
   `https://github.com/Agile-and-Adaptive-Robotics/Bipedal_Robot.git`.
6. Verify afterward: `git log --all -- '*.npz'` empty (Phase C), .git
   size reduced, `git status` matches the pre-rewrite dirty list,
   unaffected branches/colleagues' refs untouched.
7. GitHub-side full shrink may need GitHub's gc (contact support or
   wait for their maintenance) — local + push results are immediate.

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
- **Result-mat saves are FULL-WORKSPACE (Ben directive, 2026-09-12)**: bare `save(file)` — he
  loads the mat and runs any driver section from it; figures inside the mat are accepted.
  Every `*FlxPin10_results_202609*.mat` (root + Dig_out archives) now carries lowercase
  `allBPA` + `numBPA` (patched 2026-09-12 via `Dig_out\patch_FlxPin10mats_allBPA_20260912.m`;
  noT3 vintages: allBPA = [1 2 4 5], the training pool of that era).
- **Plotting convention (Ben, 2026-09-12)** in `minimizeFlxPin10mm` / `minimizeFlxPin10mmX3`
  color-scheme sections: his 3-line override block
  (`% allBPA = allBPA;` / `allBPA = [1, 2, 3, 4, 5];` / `numBPA = numel(allBPA);`),
  auto tileLabels (A),(B),… drawn bold at each tile's top-left, and the legend rule:
  even numBPA → legend in tile (1,2); odd → first empty tile (`lg.Layout.Tile =
  2*ceil(numBPA/2)`, e.g. tile 6 of 3×2). `minimizeFlx10mm` (1×2, different series per
  tile) and `minimizeExt10mm` (separate 1×1 figures) keep per-tile legends — no empty
  tile exists there; letters are auto tileLabels. Verified renders in
  `Dig_out\plotcheck_*_20260912.png`.
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

## Figure standards — PROJECT-WIDE (Ben, 2026-09-21)

Every figure produced for this project (MATLAB, Python, Illustrator
exports) must follow these rules; reference examples: the
Plot_KneeFlxPin family + `Documentation\Reports and Papers\
Knee_Torque_Test\Figures\Figure components\FlxPin_group\FlxPin_group.fig`.

1. **Page/size**: sized for regular letter paper, usable area
   7.5 x 10 in (8.5x11 with margins). Nothing smaller.
2. **Type**: minimum 10 pt everywhere (axis labels, ticks, legends,
   annotations). **No italic text** (set font.style normal; avoid
   mathtext italics — use \mathrm or plain text).
3. **Font**: Arial only (freely available in MATLAB, Adobe Illustrator,
   and Python/matplotlib on all three machines). MATLAB:
   set(groot,'defaultAxesFontName','Arial'); Python:
   rcParams["font.family"]="Arial", "mathtext.fontset":"custom" with
   rm/it/bf all Arial.
4. **Colors**: use the accessible palette from `Code\Matlab\Colors.m`
   (Paul Tol 7: #FFD700 gold, #FFB14E orange, #FA8775 coral,
   #EA5F94 pink, #CD34B5 magenta, #9D02D7 magenta2, #0000FF indigo) as
   the series palette, in that order. Greys #B0B0B0 (context) and
   light lavender (inactive circuit context) for de-emphasized
   structure, Di-Russo-style.
5. **Accessibility**: colorblind-safe by construction (Tol palette +
   distinct line styles/markers/shapes so no information is carried by
   hue alone); every figure ships with **alt text** (one file per
   figure set, e.g. <figure>_alt.txt, one block per panel: panel
   letter, what is plotted, key takeaway).
6. **Synapse/element shapes in circuit diagrams** (Ben's SNS
   convention, 2026-09-09): open circle = neuron, open triangle =
   excitatory synapse, filled circle = inhibitory synapse,
   ellipse = muscle; afferent/interneuron classes distinguished by
   Colors.m hue + label.
7. **Terminology**: refer to afferent pathways as Ia / II / Ib
   (e.g. "stance-group Ib interneuron" — the code symbol LBIN means
   exactly that; do not write "LB"). KINH = swing-gated inhibitory IN
   (PF_F1 -> KINH -> knee_ext/ankle_pf MN suppression).

## Readmes

- `Code\Matlab\Mesh_Optimization\Mesh_Optimization_Readme.md`, `Code\Matlab\Functions\Functions_Readme.md`,
  `Code\Matlab\Robot_Data\Robot_Data_Readme.md` — agent-readable copies; **the .docx versions are
  canonical for humans.** If you edit a .md copy, remind Ben to update the .docx.
- `Code\Matlab\Mesh_Optimization\Knee_Torque_revision_3\Knee_Torque_revision_3\README_revision_3.md`
  — detailed flexor BPA/route model doc (modes, how to run tests).

- Known data concern FINAL (Ben, 2026-09-12): the angle-shifted test is **#3, the
  flexor 47 cm (kf(3))** — its encoder read ~5.3° low; correction is
  **Angle-only (phiD NOT shifted)**, applied at build time in `minimizeFlxPin.m`
  (Ben moved the shift back to kf(3) himself on 2026-09-12, after briefly assigning
  it to test 4 / 40cm-tendon on 2026-09-11). HISTORY: the 2026-09-08 2brk campaign,
  dissertation flexor pick 107, and the extensor Xi1/Xi2 lock were all fit with the
  kf(3) shift; the 2026-09-11/12 `offT4` / `Xi1gtXi2` / `K2allX1_2folds` mats were
  run with the kf(4) shift — remember which shift a mat carries when re-picking.
  `minimizeFlxPin2brk.m` happens to carry the shift on kf(3) but is STALE re: the
  2026-09-13 mixed convention (tibia 1trans + hip 2trans K2=[X1,X1,X2]).

### Mixed-convention campaign COMPLETE (2026-09-13)

- Config: tibia 1trans K=[X1,X2,X1] / hip 2trans K2=[X1,X1,X2]; folds holdout
  {1,5} and {3,4} (nothing left out); Xi1>Xi2 constraint live; encoder shift on
  kf(3). Runner `Dig_out\run_FlxPin10mm_mixT_20260913.m` (first launch lost run 1
  at the save line to the driver-`clear`-wipes-runner-vars trap — fixed, rerun
  deterministic). Readout: `Dig_out\Dig_FlxPin_mixT_readout_20260913.m`.
- **Run 1 `..._20260913_2brkt_mixT_2folds.mat` (no T5 shift)**: fold 1 (train
  {2,3,4}) best Xi0 1.279 mm / 2.551e4 / 2.556e4 ratio 0.998 (boundary); fold 2
  (train {1,2,5}) Xi0 2.465 mm / 3.389e4 / 1.748e4 ratio 1.94 INTERIOR; pooled
  106 rows, median ratio 1.53. pick=1 per-test RMSE 2.196/1.363/2.281/1.187/1.499,
  FVU ≤ 0.177 — 40cm-tendon best-ever 1.187/0.023.
- **Run 2 `..._mixT_T5y5mm_2folds.mat` (T5 +5 mm bracket-y)**: pick=1 Xi0 1.572 mm
  / 2.695e4 / 2.434e4 ratio 1.107; fold-2 winner Xi0 pushed to 7.1 mm. Per-test
  RMSE 2.181/1.315/2.253/1.168/**1.556**.
- **T5-shift verdict (fixed-Xi isolation, tests 1–4 bit-identical both ways): the
  +5 mm bracket-y offset HURTS test 5** (1.499→1.602 at run-1's pick; 1.440→1.556
  at run-2's). Sign sweep: −2 → 1.69, −5 → 2.17 — monotonic degradation in BOTH
  directions from zero. Test 5's residual is NOT explained by an insertion-point
  bracket-y offset; the un-shifted pB is optimal for it. Plastic-deformation
  allowance in that axis is unsupported by the torque data.

## Hazards — do NOT open these as source

- Known data concern FINAL (Ben, 2026-09-12): the angle-shifted test is **#3, the
  flexor 47 cm (kf(3))**; +5.3° Angle-only at build time in `minimizeFlxPin.m`.
  2026-09-11/12 offT4-family mats carry the kf(4) shift instead.
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

**CAVEAT (2026-09-11):** the flexor pick below (and everything downstream of it, incl. the
extensor Xi1/Xi2 lock) was fit with the +5.3° encoder offset on the WRONG test (47 cm
instead of 40 cm-tendon). Ben ordered a corrected re-run — `minimizeFlxPin10_results_
20260911_2brkt_2trans_offT4.mat` (legacy minimizeFlxPin10mm driver) — which supersedes
these for future work once Ben re-picks. Section kept as the record of what the
dissertation text used.
**offT4 campaign outcome (completed 2026-09-12 01:47; ran 3× deterministic — identical
tables):** pick=1 = Xi0 ≈ 0 (1.1e-6 m) / Xi1 2.42e4 / Xi2 2.11e4 N/m, 265/265 filter pass,
per-test RMSE 2.21/1.37/2.61/1.65/1.49 (48/46/47/40cm-t/41), all FVU ≤ 0.214. Fold 1 =
Ben's preferred split (train {3,4,5} / holdout {1,2}): 53 candidates, Xi0 0.01–0.95 mm,
Xi1 2.11–2.23e4, Xi2 2.41–2.95e4; best-distance candidate VALIDATES better than it trains
(raw held-out 48cm 2.28/0.19, 46cm 1.60/0.064). Readout: `Dig_out\Dig_FlxPin_offT4_readout_
20260912.m`. minimizeFlxPin.m restored to Ben's kf() naming (kfCache persistent wrapper),
baseline+pick verified identical post-restore.

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
- **Display-script rehydration (2026-09-10):** dated `Vas_Pam_20mm_Result_*.mat` (rework
  onward) originally saved xBest/XiUsed but NOT ctx. NOW FIXED BOTH WAYS: (1) both drivers
  save `ctx` in their dated captures (`Opt_run_Ext` adds ctx; `Opt_run` adds ctx + a new
  `XiUsed` record alongside its legacy `Xi3` field); (2) the two current mats were
  rehydrated IN PLACE with verified fresh-built ctx — `Vas_Pam_20mm_Result_20260910_0528.mat`
  and `Bifemsh_20mm_Result_20260910_1234.mat` (gate before each write: builder Xi ==
  XiUsed and fBest/constraints reproduced exactly; temp backups + rehydrate script in
  `%LOCALAPPDATA%\Temp\rehydrate_backup`). `Vas_Pam_20mm_Result_20260910_0446.mat` was
  deliberately NOT rehydrated (its pass-1 Xi pair -12.3mm/0.281 differs from the current
  builder; attach an Xi-overridden ctx only if it ever needs displaying).
  `Results\Knee_Extensor_20mm.m` still rebuilds ctx via buildKneeExtContext20mm() with a
  builder-Xi == XiUsed guard (could now load mat ctx instead, but rebuild+guard also
  catches builder drift). The old Sep-1 `Vas_Pam_20mm_Result.mat` still carries a
  pre-rework ctx.geo — loading it with the current builder is what produced the original
  `bypassTol` error. `AnimateKneeBoneMuscle` writes its GIF next to the script folder
  regardless of cwd.
- **FVU caveat**: bio-ext 52cm FVU 2.35 > 1 at this pick — improved vs baseline 4.49 but do not overclaim.
- **Opt_run runtime reality check (Ben, 2026-09-10):** configured eval budgets (surrogateopt
  ~7000 + patternsearch ~15000) are MAXIMA — actual stages exit early on FunctionTolerance;
  Ben's observed patternsearch stages wrap in ~5 min, not 1.5-2.5 h. Extrapolate ETAs from
  measured rates and his historical runtimes, not from configured caps.

## Extensor K=[X1,X2,X2] single-bracket campaign + flexor 2brkt (2026-09-14/16, EB475WS4)

- **minimizeExtX3.m K is now [X1, X2, X2]** (was [X2,X1,X2]; old line commented in fortz) —
  Ben-directed, LIVE for all future extensor runs. Single bracket (no useB2 in this evaluator;
  arg 6 = transMode). **minimizeExt10mmX3.m** gained env hooks `EXTX3_HOLD` ('1,8' etc. =
  ONE custom fold) and `EXTX3_ALLTESTS=1` (allBPA = all 9 tests), fold-row loop fixes
  (`size(list,1)` not `length(list)` — a 1xN row list re-ran folds; `ind = (1:size(x2,1)).'`
  — `1:length(x2)` breaks on 1x4 single-row fronts; the latent length() bug is STILL in
  minimizeFlxPin10mmX3.m), and reworked plot sections: 2 figures per metric (Training /
  Validation, tiles subtitled "Training"/"Validation"). The same latent length() bug was
  fixed in minimizeFlxPin10mmX3_2brkt.m; minimizeFlxPin10mmX3.m untouched.
- **Four single-fold lock runs saved** (EXTX3_PASS=2; mats minimizeExt10mmX3_results_20260916_*):
  pick1_h18 (lock 4.354e4/1.701e4, holdout {1,8}: Xi0 −1.60mm, Xi3 0.797, mean RMSE 1.448/FVU 0.650),
  pick107_h18 (5.624e4/1.854e4, {1,8}: −1.55mm, 0.796, 1.457/0.653),
  pick1_all_h3479 / pick107_all_h3479 (allBPA = all 9, holdout {3,4,7,9}: Xi3 0.446/0.480,
  mean RMSE 1.195/1.163). **Both lock pairs converge to near-identical picks** — the locked
  pair barely matters in this configuration. Also `_20260916_1translock.mat` (full 10-fold CV,
  1trans pair: Xi0 −0.91cm, Xi3 0.842, mean RMSE 1.638). Runner/log: Dig_out\run_ExtPinX3_x122_fourcases_20260916.m.
  **Lock-pair provenance (spelled out): the runner read the pairs from TWO different mats** —
  pick1 pair (4.354e4/1.701e4) from minimizeExt10mmX3_results_20260910_noT3.mat sol_actual,
  pick107 pair (5.624e4/1.854e4) from minimizeExt10mmX3_results_20260914_pk107lock.mat
  sol_actual. Both originate in the same flexor front minimizeFlxPin10_results_20260908_
  2brkt_2trans_noT3.mat (filtered_results rows 1 and 107); each extensor mat carries its pair
  unchanged because the extensor never searches Xi1/Xi2 (EXTX3_PASS=2 lock input).
  **Ben's note (2026-09-16): he will still probably USE the two-bracket `_noT3` flexor front
  results — minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat at pick 1 or pick 107
  (the dissertation-settled line). The K=[X1,X2,X2] extensor runs and the flexor x3u runs are
  comparisons, NOT a replacement of that line unless Ben says so.**
- **minimizeExt10mmX3_results_20260910_noT3.mat rehydrated in place** (labels/allBPA/numHold/
  baselineScores added — pick/plot sections now run from it; pre-fix backup in
  Dig_out\..._BACKUP_20260916.mat; gates verified: pick = −10.12mm/0.621 settled row).
- **minimizeExt10mmX3.m plot sections**: truncation at strain < −0.03 is INTENTIONAL (Tor NaN
  rule; stretch-side festo4 unvalidated — manufacturer 2-3% limit, never characterized) —
  do NOT relax it, fake Go_OfF values, or plot extrapolation (Ben ruled 2026-09-16).
- **Flexor 2brkt (companion)**: minimizeFlxPinX3.m = original restored bit-exact;
  minimizeFlxPinX3_2brkt.m = two-bracket evaluator (bracket 2 at Pbri2 = [30.5,−103.41,0]mm
  tibia frame, K2t = [X1,X2,X1], compliance-only; screw-head CLAMP at −3.5mm tibia-X with
  contact force Nc in bpa; Xi3 = UNITLESS wrap-loss delta_L = Xi3·15mm·theta_wrap·comp²;
  30mm-circle tangency check, FLXPX3_TANGENCY=DIAG = store-not-enforce; FLXPX3_NOSHIFT=1 =
  legacy unshifted). Driver minimizeFlxPin10mmX3_2brkt.m (FX3B_LOCK1/2 env lock). x3u lock-CV
  mats 20260915: L107 DIAG = Xi0 ≈ 0.1mm / Xi3 0.150 / mean RMSE 1.587 / FVU 0.077 (best
  pinned-flexor fit on record); L1trans = Xi3 0.291. Non-_x3u 20260915 mats = superseded
  series-stiffness Xi3. Gotchas: pool workers never see client setenv after spawn (env-gated
  evaluator branches silently run the wrong mode in parfor — set env BEFORE parpool); wrap
  the Xi-factor handoff details: CHATGPT_HANDOFF.md ACTIVE WORK F (2026-09-14/16).
