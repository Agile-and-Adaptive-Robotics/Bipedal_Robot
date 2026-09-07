# AGENTS.md — Standing context for AI coding sessions in Bipedal_Robot

Loaded automatically at session start. Keep it current; keep it lean.

## Project purpose (priority order)

1. **Dissertation first** (deadline: this week, Sept 2026). LaTeX source lives in
   `Documentation\Reports and Papers\Dissertation\` (untracked); the `upload\` mirror is the
   canonical Overleaf copy. Use the `latex-overleaf` skill for Overleaf work.
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
- `Code\Arduino\`, `Code\Festo\` — embedded/valve hardware code.
- `Testing_Data\` — **important.** Immediate subfolders `2022_02_Festo\` and `2026_06_Festo\` matter.
  - `2022_02_Festo\` holds the Xi-minimizer family: outer CV drivers (`minimizeFlxPin10mm.m`,
    `minimizeFlxPin10mmX3.m`, `minimizeExt10mmX3.m`) call inner evaluators (`minimizeFlxPin`,
    `minimizeFlxPinX3`, `minimizeExtX3`, `minimizeExt`, `minimizeFlx`) that each carry their own
    `computeForceVector`/`Lok`/`fortz`. New (Sept 2026): two-bracket flexor method —
    `minimizeFlxPin2brk.m` (evaluator; (d) single transform, (e) K=[X1,X2,X1], (f) origin-side
    2nd bracket `Pbr2` (current: [-52.61, 0, 75.06]/1000, Ben-set), `USE_BRACKET2` flag).
    **Pbr2 applies ONLY to the pinned-knee flexor configuration** — the extensor evaluators
    (minimizeExtX3, minimizeExt) have their own independent bracket offsets; do not port Pbr2.
    + driver `minimizeFlxPin10mm_2brk.m` (env `FLX2BRK_MODE`
    = smoke|full, `FLX2BRK_SOLVER` = gamultiobj|surrogateopt), harnesses `crossPredictFlx.m`
    (flexor→biomimetic/extensor cross-prediction), `sweepExtX3.m` (high-Xi1 hunt, pool
    {1,2,5,6,7,8}; tests 3/4/9 EXCLUDED per Ben), `abSolver.m`, `runBatch.m` (overnight
    sequence). `Robot_Data` must be on the MATLAB path for the biomimetic evaluators.
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
  Bifemsh/Vastus-adjusted variants used for muscle torque targets).
- `Documentation\Reports and Papers\` — papers; `Dissertation\` = the dissertation (see priority 1).
- `Neuromechanical_Models\` — Animatlab CPG walker models.

## MATLAB workflow (follow exactly)

- Entry scripts are in `Code\Matlab\Mesh_Optimization\`: run `Opt_sanity.m` / `Opt_sanity_Ext.m`
  **before** the full `Opt_run.m` / `Opt_run_Ext.m`.
- Every run starts from a context struct: `ctx = buildKneeFlexorContext20mm()` or
  `buildKneeExtContext20mm()`; objectives/constraints are wrapped as anonymous functions of
  `(x, ctx)`. Flexor and extensor use different constraint machinery.
- Optimizers use `parpool(10)` + `surrogateopt` (≈7000 evals) then `patternsearch` (≈15000) —
  **hours, not minutes.** Run via the `matlab` skill (`matlab -batch`) in the background with
  output tee'd to a log; check the log tail, never babysit the run turn-by-turn.
- Path setup: self-locate the repo root, then `addpath(genpath(root/Code/Matlab))` with
  **Mesh_Optimization winning any shadowing contest** (duplicate filenames exist in Robot_Data,
  Knee_Torque_revision_3, Bone_Mesh_Plots). `Debug_RouteElim_Ext.m` shows the canonical setup.
- Results (.mat/GIF/CSV) go to `Mesh_Optimization\Results\`.
- Never casually re-run entry scripts; ask/state intent first — runs are long and results get
  overwritten.

## Current state (Sept 2026)

- Flexor: **solved and feasible** with 2 BPAs, but currently assumes double force at the same
  attachment points instead of mirrored routes (commit 0453a2f caveat — still true).
- Radius controls (commit fd5963d): `ctx.bpaRadiusMode` = "scalar" (geo.bpaRb / geo.bpaRs) or
  "bpaR" (candidate-dependent physical radii); `geo.bpaRbOffset` / `geo.bpaRsOffset` defaults 0.
  Three-radius scheme: `bpaRb` = 20 mm nominal pWrap standoff, `bpaRs` = 16 mm min tibia
  clearance, `wRap` = 25 mm Xi3 bend-length radius.
- Extensor: hit-or-miss; `Debug_RouteElim_Ext.m` (untracked) is a **read-only** diagnostic for
  the route-elimination state machine in `buildDistalRingLocation20mm.m` (2244 lines — the
  core; read its docstring first: `x = [p1(1:3), pEnd(1:3), rest, tendon]`).
- Open WIP (uncommitted): `Debug_RouteElim_Ext.m` + logs, modified
  `buildDistalRingLocation20mm.m` / `buildKneeExtContext20mm.m`, staged deletion of
  `Knee_Torque_revision_3\...\MonoPamDataExplicit_balanceX3.m`.

## Readmes

- `Code\Matlab\Mesh_Optimization\Mesh_Optimization_Readme.md`, `Code\Matlab\Functions\Functions_Readme.md`,
  `Code\Matlab\Robot_Data\Robot_Data_Readme.md` — agent-readable copies; **the .docx versions are
  canonical for humans.** If you edit a .md copy, remind Ben to update the .docx.
- `Code\Matlab\Mesh_Optimization\Knee_Torque_revision_3\Knee_Torque_revision_3\README_revision_3.md`
  — detailed flexor BPA/route model doc (modes, how to run tests).

## Hazards — do NOT open these as source

- `Solid_Models\Biomimetics_2022-Knee_Test\Point_cloud\Tibia_copy.txt` (7.1 MB point cloud);
  `Spine_Mesh_Points.txt` (172 KB, duplicated in 3 places); `HX711*\esempio.txt` (1 MB);
  any `.mat` in `Previous Optimization Code\Trial Results\` (up to 95 MB).
- `.asv` files = stale MATLAB autosaves that mislead; the real source is the `.m` beside them.
- Logs: grep/tail them, never read whole. Point clouds/binary .mat: process with MATLAB, never
  ingest.

## Agent efficiency rules (credit/token hygiene)

- Delegate broad searches and audits to subagents; keep main-session context small.
- Fresh session per task; long MATLAB runs in background; plan mode before expensive execution.
- Prefer COM/API automation (SolidWorks skill, Overleaf file edits, Zotero local HTTP) over
  screenshot-driven GUI automation.
- Model strategy: default GLM-5.3-Flash; escalate to GLM-5.3 only for hard debugging or
  architecture decisions. Peak billing is Mon–Fri 11 PM–3 AM Pacific — schedule expensive
  (GLM-5.3) work before 11 PM or on weekend nights; evenings/weekends are 50% off-peak.
