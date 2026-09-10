# ChatGPT handoff — Knee minimizer / Xi-factor work (Bipedal_Robot, Sept 2026)

This file briefs any AI assistant (ChatGPT or otherwise) working on this project while the
primary assistant (ZCode) is unavailable. Read it fully before doing anything. The repo-wide
`AGENTS.md` next to this file has additional standing context — both are plain markdown.
**TWO ACTIVE WORK THREADS (both 2026-09-09): (A) AnimatLab .aproj arrow audit — section
below; (B) MuJoCo gait2392 spinal cord network in `Code\MuJoCo_SNS\spinal\` (read its
`DESIGN.md` first — rhythm layer verified; leg-DoF NaN blocker open with prioritized
suspects).** Also: laptop-session mining findings in `Testing_Data\2022_02_Festo\HANDOFF_laptop_20260908.md`.

## ACTIVE WORK B — MuJoCo gait2392 spinal network (started 2026-09-09, easteregg2)

Two-level RG+PF spinal network (SNS-Toolbox 1.5.2) for the converted
`gait2392_simbody` MJCF: 92 MN pools + Ia/II/Ib afferents each, per-leg
half-center RG, 4 phase-shifted PF groups/leg, stance-gated Ib load sharing,
descending DRIVE + balance inputs, static-opt standing posture injected as
per-MN bias. **Rhythm layer VERIFIED** (0.885 s period, antiphase −0.88);
closed-loop walk still NaNs in the leg joints.

- Run with `D:\Anaconda\envs\myo\python.exe` (Anaconda at `D:\Anaconda`,
  NOT on PATH; env `myo`, not `myoconv`), cwd `Code\MuJoCo_SNS\spinal`.
  Read `spinal\DESIGN.md` first: architecture, verified/unverified split,
  open problems. `audit_signs.py` (anatomical sign audit — must stay all-OK
  after ANY model change), `check_rhythm.py`, `runner.py`, `diag_nan.py`
  (first-NaN finder), `fit_synapses.py` (IK/SO → synergy NMF + per-phase
  NNLS back-solve scaffold).
- Model repairs are applied on the fly in `apply_harness()` (runner.py):
  hip_flexion/hip_adduction MJCF hinge axes WERE flipped vs OpenSim
  (negated, both sides); rect_fem lacked a patella wrap and pulled the knee
  into FLEXION (rerouted over the vastii's vas_med-P4 patella-tracking via
  point); quad_fem/gem/peri pruned (`PRUNE_MUSCLES`, Ben's list).
- **Do NOT weld the conditional pathpoints** (massless slide bodies with
  equality couplings): an earlier weld froze several muscle moment arms to
  zero. They stay coupled; massless bodies handled via compiler
  `boundmass`/`boundinertia`. The sign audit must run DYNAMICALLY (full
  activation → joint acceleration) — static `actuator_moment` misses the
  coupling paths and misdiagnoses.
- Remaining blocker (dynamics, not kinematics — audit is clean): NaNs at
  2 ms even under a rigid pelvis rig. Suspects in order: mesh foot contacts
  (try primitive collision), equality forces on the 0.01 kg pathpoint
  bodies, near-zero hinge damping, MN-command slew limiting.
- Speed-control design: descending DRIVE + presynaptic reflex-gain
  modulation (`params.MOD`); grounded in Bunz/Ijspeert/Schmitt 2026 Sci Rep
  and Ben's Zotero "Sensory Afferent Database" collection. Read Di Russo/
  Ijspeert/Bouri 2023 J Neural Eng before any novelty claims. MuJoCo gotcha:
  muscle Fmax = `actuator_gainprm[:,2]`.

# ACTIVE WORK A — AnimatLab .aproj neural surgery (started 2026-09-09, CONTINUES at 8am PT session)

Ben's asks, in progress: (1) complete the RH-side feedback/connections, (2) add
gastrocnemius, biceps femoris long head, rectus femoris, semimembranosus to each leg and
drive them from the existing pattern-formation layers, (3) connect LH and RH with a rhythm
generator layer. Ben also wants TFL eventually (needs a frontal-plane abduction DOF the Li
legs don't have yet) and the LH_HipZ/RH_HipZ labels renamed LH_Hip/RH_Hip (cosmetic, not
done). **Ben added some missing connections manually in the GUI on 2026-09-09 evening and
reports more are still missing — auditing drawn-arrow completeness against the standard
section is the FIRST task of the next session.**

## Current file state

- `Neuromechanical_Models/Biped_2xCPG_wSubs/Biped_2xCPG_wSubs.aproj` — GUI-verified clean
  open on 2026-09-09 (zero dialogs; earlier builds threw 1–15 dialogs per open).
- 20 muscles (12 original + 8 new biarticular: gas/bflh/semimem/rf per leg, Gait2392 Fmax
  applied: 2241/896/1288/1169 N), 122 neurons, 261 synapses, grids off everywhere.
- Completed: RH half-center + coordination wiring (mirrored from LH), RH knee/ankle
  muscle-drive adapters repointed off Renshaw cells onto MNs, LH↔RH commissural RG
  inhibition (4 synapses via OffPages on the top page), 2 nA tonic on L RG ext,
  full Deng-style chains for the 8 new muscles (MN, afferent Ia, Ia E, Ib, Muscle+SR
  nodes, 3 adapters each), Gait2392-derived Fmax. Gas Ia→Ia correctly retargeted to the
  ankle Ia F interneuron.

## Known remaining gaps (next session, in order)

1. **Missing drawn arrows** — Ben found more after the fix pass. Likely cause: links
   whose two endpoints never shared a page were left functional-but-undrawn. Audit each
   page's drawn arrows against the standard section; draw the missing ones via new
   OffPage instances on a page containing both endpoints.
2. **RG does NOT oscillate** — 2–10 nA tonic on L RG ext latches (L ext depolarized,
   everything else suppressed, even over 20 s). Sweep drive 20–40 nA; try the paper's
   10 nA 1 ms pulse; apply Deng Table A2 PF→MN conductances (hip 2.565/3.632, knee
   4.93/1.516, ankle 4.054/4.522 µS) instead of the 0.5 µS defaults.
3. **16 placeholder attachments** — new-muscle origin/insertion points were copied from
   anatomically similar sites; Ben must position them. Then set RestingLength = TSL + OFL
   (Gait2392: Gas 0.45 m, BFlh 0.435 m, Semimem 0.439 m, RF 0.424 m), LengthTension
   window = TSL+0.5·OFL to TSL+1.5·OFL, Kse ≈ Fmax/(0.033·TSL).
4. No Renshaw cells on new MNs; no II afferents on the new hip-spanning muscles.
5. TFL addition — needs a frontal-plane abduction DOF added to the Li legs first.
6. LH_HipZ/RH_HipZ → LH_Hip/RH_Hip label rename (cosmetic, pending Ben's OK).

## Toolchain (in `Neuromechanical_Models/Biped_2xCPG_wSubs/tools/`)

Full deterministic pipeline, in this order from a pristine baseline (git
`6437441:Neuromechanical_Models/Biped_2xCPG_wSubs/AnimatLab_backups/Biped_2xCPG_wSubs.aproj.bak_20260909_preRH`
is the pre-surgery original): `build_A2.pl` (RH wiring) → `build_B2.pl` (muscles) →
`repair_dup_ids2.pl` (re-GUID child IDs) → `fix_gas_ia3.pl` (Gas Ia retarget) →
`rebuild_pages_v3.pl` (regenerate all 15 page drawings from the standard section; drops
stale, dedups, recomputes all Org/Dst, grids off, exact counts) → `move_links_by_page.pl`
(place generated links into the fragment whose page draws them) → `set_fmax2.pl`
(Gait2392 Fmax). Verify with `verify_placement.pl` (link placement), `verify_final.pl`
(battery), `shape_diff.pl` (element-shape diff vs a reference file).
**If any input changes, re-run the WHOLE chain — never patch a patched file.**

## CRITICAL format gotchas (each cost hours — do not relearn)

- The .aproj has TWO layers: the standard section (functional model: neurons, links,
  muscles — full GUIDs) and per-page `<DiagramXml><![CDATA[...]]></DiagramXml>` AddFlow
  drawings (GUI only, objects referenced by `<Tag>`). The whole-file XML parse does NOT
  validate CDATA contents — validate every page CDATA individually as its own XML doc.
- Subsystem fragments are SCATTERED top-level blocks after `</NeuralModules>` (not
  nested); hierarchy is via SubSystemID references. A drawn link must live in the same
  fragment as the page that draws it, else the GUI shows nothing and clicking throws
  "No item with ID". A page's OffPage instances are that page's local references.
- Depth-counting `<Node>`/`</Node>` MUST (a) skip `<![CDATA[...]]>` spans (they contain
  drawing `</Node>` closers with no matching opens in that counting scheme) and (b) count
  both `<Node>` and `<Node attr...>` opens. Getting this wrong silently redirects
  insertions into page CDATAs → "GetAttribute(Int32 i) out of range" pop-ups (one per
  broken page).
- Page rebuilds must emit, in order: prologue (through `<AddFlow ...>`), **`<Version>`
  header** (inside the body, before the first node — dropping it throws the same
  GetAttribute error), nodes, links, tail (fillers + `</AddFlow></Diagram></Root>`).
  Missing tail = "An error occurred while deserializing the xml data."
- Cloning any block: re-GUID the block's own ID AND every child-object ID inside
  (StimulusTension, LengthTension, Gain, CaActivation, CaDeactivation, PID) or the C++
  sim throws "Attempted to add an object with the same ID twice". Strip inherited
  InLinks/OutLinks. Escape `&` as `&amp;` in Text. Patch `Value` AND `Actual`
  attributes together. Regex replacements that build XML MUST use the /e flag.
- The recurring silent-corruption warning: lvalue `substr($x,$p,length($seg)) = $seg`
  with a GROWN $seg eats trailing bytes (always use the original range length).
- **Verification protocol (MANDATORY before handing the file to Ben):** launch
  `D:\Program Files (x86)\NeuroRobotic Technologies\AnimatLab\bin\AnimatLab2.exe` with
  the .aproj path as argument, wait ~15 s, inspect the app state for an "Error" dialog
  window. Each broken page throws its own dialog at load, so zero dialogs ⇒ all 15 pages
  deserialized. Status bar must show "Load project complete". Then `taskkill /IM
  AnimatLab2.exe /F` — NEVER save from the GUI (it overwrites the file with stale memory).
  **Never deploy while Ben may have AnimatLab open** — his open session can't see file
  changes and a save from it clobbers the deployed fixes.
- Run headless sims via `bin\AnimatSimulator.exe <path>\file.asim` (.asim only, never
  .aproj). Chart files (e.g. "Rhythm Generator.txt") land in the sim file's folder; chart
  `<EndTime>` caps collection independently of SimEndTime.

## Background

The .aproj is Deng 2019 (Biomimetics 4(1):21) two-layer CPG (RG→PF→MN, Ia/Ib/II
afferents, Renshaw) ported onto the Li biped biomech; Ben+Connor's walker paper documents
the lineage. The old `_Standalone.asim` exports are stale; Ben exported a fresh
`Biped_2xCPG_wSubs_Standalone.asim` on 2026-09-09 (runs clean headless; contains a
Rhythm Generator chart with RH_RG columns + L Hip charts). Tuning ladder: kinematics →
virtual-walker ground walking → stable walking → transitions (details in AGENTS.md and
the animatlab skill's walker-cpg-architecture.md).



## Latest completed session — dissertation simulation figures, 2026-09-09

Ben requested the handoff and ended this session. The local dissertation edits
and figure review are complete. Read the dated **2026-09-09** section appended to
`CHATGPT_REPORT.md` for the exact scope and remaining work; the older report's
Overleaf compile counts do not apply to these new edits.

- Edited `Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/20-methods.tex`
  and `30-results.tex`. Preliminary simulation Methods now follows the historical
  AnimatLab walker section, after the actuator/joint methods.
- Added AnimatLab architecture/GUI, native MuJoCo model renders, and native
  Simulink library/reflex figures. Results contains a verified, explicitly
  preliminary bilateral hip/knee/ankle plot from the saved phase-1 run.
- Review: `Documentation/Reports and Papers/Dissertation/Notes/neuromechanical_figure_review.pdf`
  (seven pages). Detailed provenance and reproduction scripts are in the same
  `Notes` directory; figure assets are in `ProofFinal/figs/Preliminary/`.
- **Pending:** compare against the current Overleaf chapters, integrate the local
  edits/assets, then compile and inspect final float pagination. This session did
  not update Overleaf, the dissertation ZIP, or the full dissertation PDF.
- No new dynamics/optimization runs or scientific model changes. No commit or push.
  Preserve unrelated working-tree changes from other sessions.

## Latest completed session — Gait2392→MuJoCo conversion + SNS-Toolbox connection, 2026-09-08/09 (easteregg2)

This is the session ACTIVE WORK B builds on: it produced the converted MJCF that
`Code\MuJoCo_SNS\spinal\` drives, plus the `myo` env it runs in. Full environment
detail and Windows pitfalls live in `AGENTS.md` → "OpenSim / MyoConverter /
SNS-Toolbox on easteregg2". Highlights:

- **Stock `gait2392_simbody.osim` fully converted to MyoSuite/MuJoCo** — outputs in
  `Solid_Models\OpenSim\myosuite_gait2392_simbody\`: `gait2392_simbody_cvt3.xml`
  (vehicle of record), `gait2392_simbody.pdf` (per-muscle validation report),
  cvt1/cvt2 intermediates. Regenerate via `D:\GitHub\myoconverter\convert_gait2392.py`
  (~90 min). `gait2392_cvt3_pelvispinned.xml` (same folder) = pelvis joints +
  keyframe stripped, joint limits enabled — use for leg-swing demos.
- **`gait2392_robotbody.osim` is NOT convertible as-is**: right leg carries
  placeholder muscle paths (15 muscles with changed path-point counts, 19 muscles
  disabled) over stock Thelen parameters → path lengths off 7–480%; OpenSim
  `computeInitialFiberEquilibrium` fails (tfl_r) and force maps would be wrong
  (Ben: "wacky answers"). Conversion needs retuned fiber/tendon lengths first —
  Ben's modeling decision. `gait2327.osim` (GENERATED by
  `D:\GitHub\myoconverter\build_gait2327.py`) = stock params + robotbody paths,
  all 92 enabled, for path-only studies; robotbody's left leg is byte-identical
  to stock (it is a left-leg reference model).
- **SNS-Toolbox 1.5.2 installed in `myo`** (needs `--no-deps` + CPU torch +
  graphviz; the wheel mixes two class trees — use legacy `sns_toolbox.neurons`/
  `connections` classes with `Network.compile()`; details in AGENTS.md).
  `D:\GitHub\myoconverter\sns_gait2392_cosim.py` proves the co-sim loop
  (Nourse 2023 two-layer CPG → sigmoid → `mj_data.act` → `mj_step` → tension
  Ib feedback; ~4x real time, MNs antiphase −0.97; outputs
  `sns_gait2392_cosim.npz/.png`).
- **Three Windows pitfalls fixed** (documented in AGENTS.md, do not relearn):
  conda-MKL numpy native crash 0xc06d007e in np.dot (use pip OpenBLAS numpy
  1.21.6); `if __name__ == "__main__":` required in every MyoConverter driver
  (step-3 multiprocessing spawn re-imports the main script); zombie child
  processes hold the output log open (WinError 32) — kill before rerunning.
- No commits/pushes; new files live outside the repo (`D:\GitHub\myoconverter`,
  `D:\GitHub\Two_layer_CPG_SNS_Toolbox`) or in gitignored output folders.

## The person you are helping

Ben Bolen, Mechanical Engineering PhD researcher (PSU, Agile and Adaptive Robotics Lab).
He uses GitHub Desktop (not git CLI), and thinks in terms of the physical test hardware.
Be concrete, show numbers, don't pad.

## The science, briefly

He identifies correction factors (Xi0–Xi3) that map an ideal rigid-body + pneumatic-muscle
(BPA) model onto measured knee-test data, then feeds those into placement-optimization
scripts (`Code\Matlab\Mesh_Optimization\Opt_run*.m`) to design legs whose two 20 mm BPAs
meet or exceed human monoarticular muscle torque. Data: pinned-knee flexor and extensor
tests (5 and 9 length-variants), plus two biomimetic-knee validations.

**Advisor's requirement:** one (Xi1, Xi2) pair consistent across all four configurations
(pinned flexor, pinned extensor, biomimetic flexor, biomimetic extensor). Ben believes a
higher Xi1 is needed for the mesh optimization to reach human-level torque.

**Critical nuance:** Xi1/Xi2 are *effective system-stiffness* parameters — they lump the
bracket, fixtures, AND the cable winch of the test mechanism. "Bending" is a simple
Hooke-law spring in N/m (no EI/L length dependence). Ben shorthand-calls them "bracket
axial/bending stiffness" but they are not literal beam stiffnesses. The bracket reference
points (Pbr, Pbri, Pbr2) and frame conventions are interpretation choices; moving a point
changes the identified stiffness dramatically (one 4 cm move swung flexor Xi1 by 16x).

## Current state (2026-09-08 late, after the encoder-corrected CV rerun + biomimetic hand-tune)

- **Encoder question RESOLVED:** the angle-shifted test is the pinned-flexor **47 cm test
  (BPA #3)**; its encoder read ~5.3° low. `minimizeFlxPin2brk.m` adds +5.3° to that test's
  experimental Angle column at build time — **phiD/model angles are NOT shifted** (an
  earlier version mistakenly shifted both; fixed late 2026-09-08). The current CV campaigns
  exclude BPA #3 from training entirely.
- **Frame-convention question closed mathematically:** for y-symmetric stiffness arrays
  (K_x = K_z, e.g. [X1,X2,X1]) the 1-transform and 2-transform conventions agree to ~1e-17 —
  the frame math is not what separates the arms. What still differs is Ben's per-convention
  ORIGIN-bracket stiffness ordering (below), from his space-frame-z buckling observation.
  Both "arms" are carried forward in parallel; the pick is made by cross-configuration
  consistency, not by rerunning the pinned-flexor CV (its likelihood valley is flat —
  gamultiobj stochasticity alone swings Xi1 5.4e4 → 2.1e6 at equal fit).
- **Biomimetic-flexor hand-tune chain** (`Dig_FlxBio_dubfilt` → `_handtune` → `_refine`;
  mats/logs in `Dig_out\`) converged on **Xi0 ≈ +12 mm, Xi1 ≈ 5e5, Xi2 ≈ 1e4** (refine
  winner, score 0.057); the whole top-12 sits in Xi0 8–12 mm, Xi1 2e5–5e5, Xi2 8.5e3–1e4.
  Xi2 is consistent (~1e4–1.6e4) across configurations; Xi1 remains the flat, poorly
  identified one. The pinned-flexor driver's bounds and initial population were re-centered
  on this region on 2026-09-08.
- **Extensor side unchanged:** `minimizeExt10mmX3_results_20260907_2trans.mat`, Pbr = rib
  midpoint. Xi3 remains the helpful dial (extensor pool RMSE 1.41 at Xi3=0 → ~1.01 at
  Xi3=0.35 at flexor stiffnesses; best-generalizing folds carry Xi3 0.29–0.42, Xi0 −10 to
  −16 mm). Raising Xi1 above the flexor value monotonically worsens the extensor fit
  (tested to 100x).

## Update (2026-09-09, desktop session — extensor-side deltas)

- **minimizeExt10mmX3.m is re-pointed**: its Xi1/Xi2 lock now loads
  `minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat` (corrected-encoder 2trans
  flexor front; pick pair 4.35e4 / 1.70e4), NOT the legacy 20260730 mat. The
  "Extensor side unchanged" bullet above is thereby outdated.
- **Extensor Xi0 sign rule (bit once):** extensor Xi0 = NEGATIVE of the flexor value
  (bounds [−2,0] cm; transfer `-g(1)`). Screens that pass positive flexor Xi0 into the
  extensor are wrong.
- **minimizeExtX3.m now has transMode** (arg 6 / env `EXTX3_TRANS`; 1trans = pitch-only
  Thbr + `[norm(xy),0,z]` origin form; 2trans = two-rotation + `[norm,0,0]`; K =
  [X2,X1,X2] in both).
- **Pbr variant test (lower bolt hole vs rib midpoint), 1trans, fronts captured:**
  essentially IDENTICAL — pick Xi3 0.603 vs 0.621, Xi0 −11.6 vs −10.1 mm, filter 83/90
  vs 72/90, bio-ext 2.204/2.369/4.727 vs 2.195/2.350/4.727. Extensor Xi identification
  is robust to the bracket-point choice. Git archaeology: the ORIGINAL "lower bolt hole"
  is `[-6.26, -29.69, 75.06]` (active through the "results perfect" era); `[-2.65,
  -54.71, 75.06]` later inherited that comment and is now labeled "lower section".
- **Extensor front collapse:** with Xi1/Xi2 locked, the 72-row front holds only 4 unique
  solutions — Xi0 −6..−12 mm, Xi3 0.50–0.78 (all >> 0.02), pinned RMSE 0.96–1.64, bio-ext
  RMSE 1.99–2.26. Best bio-ext member: Xi0 −6 mm, Xi3 0.504.
- **Xi3 > 0.35 note:** every screen leader hit my 0.35 grid cap; the extensor CV's own
  fit chose Xi3 = 0.62. The Xi0↔Xi3 seesaw (both shorten the path) means high Xi3 with
  large |Xi0| is double-counting — Ben is skeptical of high Xi3; judge on plots.
- **minimizeExt.m (biomimetic) restored to committed state** — do not add transMode or
  touch it; Pbr variants are pinned-side only.

## Active code/model choices (do not silently change)

- `minimizeFlxPin2brk.m` — two-bracket pinned-flexor evaluator; signature
  `(Xi0,Xi1,Xi2,idx_val,useB2,transMode)`.
  - Insertion (tibia) bracket: two-rotation frame, **K = [X1, X2, X1]** for both arms
    (the Sept-8-morning arm-1 runs predate this and used [X1,X2,X2]).
  - Origin bracket at **Pbr2 = [-52.61, 0, 75.06]/1000** (pinned-FLEXOR only — never port):
    **1trans → pitch-only frame, K2 = [X2, X1, X2]**; **2trans → two-rotation frame,
    K2 = [X1, X1, X2]** (Ben, 2026-09-08).
  - `useB2=false` gives the single-bracket ablation. `fortz` takes transMode. Per-BPA
    structs carry `eA2` (origin-bracket deflections) and `unitD_p` (deformed force
    direction). 47 cm test: +5.3° on experimental angles only.
- `minimizeFlxPin10mm_2brk.m` — CV driver; everything runs through env vars:
  `FLX2BRK_MODE` (smoke|full), `FLX2BRK_SOLVER` (gamultiobj default | surrogateopt),
  `FLX2BRK_TRANS` (1trans|2trans), `FLX2BRK_ALLBPA` (e.g. '1,2,4,5'), `FLX2BRK_TAG`.
  Current bounds: Xi0 ∈ [0, 2.0] cm; Xi1 ∈ [3e4, 1e6]; Xi2 ∈ [5e3, 2e4] (log10 space);
  initial population pinned to 0.5–1.5 cm / 5e4–5e5 / 7e3–1.5e4. Flexor labels:
  BPA 1–5 = 48cm, 46cm, 47cm, 40cm-tendon, 41cm.
- Results naming (Ben-mandated):
  `minimizeFlxPin10_results_<yyyymmdd>_2brkt_<1trans|2trans>[_<tag>].mat` and
  `minimizeExt10mmX3_results_<yyyymmdd>_2trans.mat`. Current vehicle of record: the four
  20260908 `noT3` / `noT3noT5` mats. Every front that included the 47 cm test in training
  (the 20260907 mat + full/noT5 variants) is archived in `Dig_out\old_T3_results\`.
  Many old scripts load results by exact filename — never rename historical files.
- `minimizeExtX3.m` — pinned-extensor evaluator, one bracket (two-rotation),
  `Pbr = [-3.84, -46.44, 62.5]/1000` (rib midpoint). Line ~236 passes `Xi0` (not `[]`)
  into Contraction — that fix matters, do not revert it.
- `Dig_*` analysis harnesses (cwd = `Testing_Data\2022_02_Festo`; Functions,
  Functions\ModernRobotics, Robot_Data on the path; outputs to `Dig_out\`):
  `Dig_crossPredict`, `Dig_CVpatterns`, `Dig_ExtPinX3_CV`, `Dig_ExtPinX3_Xi3map`,
  `Dig_FlxPin_2brkt`, plus `Dig_FlxPin_2brkt_picksScan` and
  `Dig_allbpaNumHoldScan` (both still UNTESTED end-to-end), and the biomimetic chain
  `Dig_FlxBio_dubfilt/_handtune/_refine/_refine2` (refine2 STAGED, not yet run). New
  (desktop, 2026-09-09): `Dig_ExtPinX3_screen` (flexor candidates × Xi3 grid through
  both extensors), `Dig_ExtPin_frontBio` (extensor front → bio-ext scan), 
  `Dig_FlxPin_frontScan` (flexor front → bio-flexor), `Dig_FlxPin_cornerCheck` (direct
  Xi check vs pinned baseline). `Dig_FlxPin_2brkt_plots` was DELETED 2026-09-09 — the
  driver's RESULTFILE resume path covers plotting (load mat, set PICK, run). Caveat: the
  Dig_FlxBio_* scripts hardcode
  `D:/GitHub/...` paths (written on easteregg2), and dubfilt pools the pre-archive
  8-front layout (full/noT5 mats have since moved to `old_T3_results\`).

## Rules of engagement

1. Do NOT commit or push anything. Ben reviews via GitHub Desktop.
2. Do NOT launch long optimizer runs (full CVs take 45 min–2 h+ each) without saying so
   explicitly in your reply.
3. Two machines: **easteregg2** (`D:\GitHub\Bipedal_Robot`, MATLAB R2025a, 10 cores) and
   the laptop **DESKTOP-5Q16KE9** (`C:\Users\Ben\Documents\GitHub\Bipedal_Robot`,
   MATLAB R2025b, 6 cores — cap parpool at 6 there). Scripts run with
   cwd = `Testing_Data\2022_02_Festo`, with `Code\Matlab\Functions`,
   `Code\Matlab\Functions\ModernRobotics`, and `Code\Matlab\Robot_Data` on the path.
4. Never open: point-cloud `.txt` files, `*.mat` in `Previous Optimization Code`, `.asv`
   files (stale autosaves), or any log file beyond tail/grep.
5. If you change a model constant (bracket point, K order, bounds, angle correction),
   say so loudly and record the OLD and NEW values in your report.
6. Plot-quality bar: journal-publication ready if you make plots; otherwise make none.
7. **2brk drivers keep Ben's "%% Pick best solution (later, flexible)" section VERBATIM**
   (lowercase `pick`, commented `sol_actual` lines, hand-editable hardcoded k1/k2/k3, the
   two disp tables and two Mean fprintf lines) — only the evaluator call line may carry
   the extra evaluator arguments. Do not rework it into a strict or mean-of-all-tests
   version. The driver's `RESULTFILE` line loads a saved mat and skips the CV (set '' to
   run fresh); `TRANSMODE` must match the loaded file.

## Open questions you may be asked to work on

- Which arm (1trans vs 2trans K2 ordering) generalizes across configurations? Decide by
  cross-configuration consistency and Lm_p vs Lm_h length-match plots — not by rerunning
  the flat pinned-flexor CV.
- Can one (Xi1, Xi2) satisfy the advisor across all four configurations? Xi2 looks
  agreeable (~1e4–1.6e4 everywhere); Xi1 is the problem child (flat in the pinned-flexor
  fit at 5.4e4–2.1e6, biomimetic hand-tune prefers 2e5–5e5, extensor fit worsens with
  high Xi1).
- Xi0 bounds: lb is still 0, but the biomimetic evidence now favors clearly positive Xi0
  (+8 to +12 mm); earlier fronts pinned at lb=0. Whether to allow lb < 0 is open.
- surrogateopt vs gamultiobj for the Mesh_Optimization iii solver (the CV drivers settled
  on gamultiobj; the Opt_run question is open).
- Physical plausibility: is Xi1 = 2.1e6 N/m credible for an Onyx FDM bracket?

## When you are done — leave this for ZCode

ZCode (the primary assistant) reads this repo on Ben's main machine and will integrate your
work. Before finishing, either tell Ben to paste this to ZCode, or (preferred) write a file
`CHATGPT_REPORT.md` at the repo root containing:

1. **Files created or modified** — exact paths, and for each: what and why.
2. **Model changes** — any bracket point, transform, K order, bound, or flag that changed,
   with OLD and NEW values.
3. **Runs performed** — script names, configs, result file names (following the naming
   convention), and wall-clock time.
4. **Numbers and conclusions** — the Xi values found, fit metrics vs baseline, and what
   you recommend.
5. **Unfinished business** — anything attempted but not completed, with the error if any.

Do not commit; Ben handles git. ZCode will verify, re-test, and fold the work in.
(The current `CHATGPT_REPORT.md` is the 2026-09-08 Overleaf dissertation session — append
a new dated section rather than overwriting it.)
- **Front reproducibility (2026-09-09/10):** the noT3newXi-pair extensor CV rerun with a
  fresh GA seed reproduced the front EXACTLY (same 4 unique members, same ranking) — the
  collapsed 4-member family is stable. Captured front:
  `minimizeExt10mmX3_results_20260910_noT3.mat` (scan it with `Dig_ExtPin_frontBio.m`).
- **noT3 naming:** the desktop's current noT3 pair (new bounds + corrected 47 cm) uses
  plain `noT3` tags; the old-bounds campaign pair is archived in `Dig_out\old_bounds\`.
  If the laptop session also produced 20260908 noT3 mats, disambiguate before trusting
  either.
- **Results tracking:** Ben will create his own pivot table for tracking Xi results across
  evaluators/configurations (stated 2026-09-10). Don't build unprompted tooling; when he
  defines the format, maintain/populate it. Until then: present results as simple tables
  with all 3 GoF (RMSE, FVU, MaxResidual) and always name the source .mat file.
