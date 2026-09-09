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
  Bifemsh/Vastus-adjusted variants used for muscle torque targets). Also `gait2327.osim`.
  MyoConverter outputs (`Gait2392_Robotbody\myosuite_gait2392_robotbody\`,
  `Solid_Models\OpenSim\myosuite_gait2392_simbody\`) are **regenerable** (~116 MB, Geometry =
  stock OpenSim STL copies) — gitignored 2026-09-08; rebuild with MyoConverter instead of
  restoring from git.
- `Documentation\Reports and Papers\` — papers; `Dissertation\` = the dissertation (see priority 1).
  `Documentation\Research Notes\` — agent-written lit reviews with V/L/U provenance tags
  (board-sports/gymnastics motor control; kickflip neuromuscular hypothesis).
- `Neuromechanical_Models\` — Animatlab CPG walker models. The `_Standalone.asim` exports
  are STALE (predate the subsystem reorg); the active phase-1 file is the instrumented
  working copy `Biped_2xCPG_wSubs\walk new new tester added 2 axis_phase1.asim` (2026-09-08).
  A fresh standalone export from the current .aproj is still pending; DataTool_*.txt are
  run byproducts. **2026-09-09 .aproj surgery (UNVERIFIED in GUI — Ben opens it next):**
  RH-side wiring completed (RH_RG half-center, RG→PF drive, hip Ia/II + ankle Ib/II
  coordination — mirrored from LH), LH↔RH RG commissural inhibition added via 4 new
  OffPages on the top page (+ 2 nA tonic kick on L RG ext), and 4 biarticular muscles per
  leg added (Gas/BFlh/Semimem/RF: physical bodies+receptors+attachments with TEMPLATE
  gains and placeholder attachment spots, plus full Deng-style neural chains driven from
  the existing PF layers). Also fixed a pre-existing bug: RH knee/ankle muscle-drive
  adapters were wired from Renshaw cells, now from MNs. Pre-surgery backup in
  `Biped_2xCPG_wSubs\AnimatLab_backups\`.
- `Code\Matlab\SNS_Simscape\` — SNS neuron block library (`SNS_Library.slx`: non-spiking RC
  neurons, E/I synapses, Ia/Ib afferents) + `KneeReflexDemo.slx` (antagonist BPA knee reflex
  demo; runs in plain Simulink). See its README for the two open blockers: the LAPTOP's
  MATLAB license lacks **Simscape Multibody**, and the "Simscape Multibody Link" SolidWorks
  add-in is not installed — `import_simscape_when_ready.m` finishes the CAD→Simscape import
  once both exist. Tendon parts (`Tendon_Extensor/Flexor.SLDPRT`) are built but NOT yet
  inserted into `09_BA_003`; GUI steps in the README.
- Repo root: `CHATGPT_HANDOFF.md` (brief for other AI assistants when ZCode is unavailable)
  and `CHATGPT_REPORT.md` (their report back; 2026-09-08 edition covers Overleaf
  manuscript-status edits) — keep both current when work is handed off.

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
