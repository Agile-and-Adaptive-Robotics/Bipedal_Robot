# ChatGPT handoff — Knee minimizer / Xi-factor work (Bipedal_Robot, Sept 2026)

This file briefs any AI assistant (ChatGPT or otherwise) working on this project while the
primary assistant (ZCode) is unavailable. Read it fully before doing anything. The repo-wide
`AGENTS.md` next to this file has additional standing context — both are plain markdown.

## REPO SIZE REDUCTION OBJECTIVE (Ben, 2026-09-15) — READ BEFORE ANY GIT OPERATION
Ben wants the repo shrunk across the board (~13.8 GB on disk) and is WARY of
re-clone-requiring operations. If your task involves ANY git history operation
(filter-repo, amend, push --force, gc/prune, bulk deletion): read the
"REPO SIZE REDUCTION" section in `AGENTS.md` FIRST, present its MANDATORY
SAFEGUARD CHECKLIST to Ben, and get his explicit go PER STEP. Never run
destructive git operations unprompted; do not assume the ~500 MB npz commit
(a8746f9, branch KneeTestSetup_BenBo_stw, pushed) or the 5.8 GB UNTRACKED
`SADb_audit\pdf_staging\` folder may simply be deleted — the npz strip needs
backup + force push + all-other-clones-re-clone coordination, and pdf_staging
belongs to another session's upload workflow.

**SEVEN ACTIVE WORK THREADS: (0) AnimatLab .aproj wiring repair — PAUSED MID-FIX, read
`Neuromechanical_Models\Biped_2xCPG_wSubs\tools\CONTINUE_HERE.md` FIRST for exact state,
fix spec, and the GUI-verification protocol; (A) AnimatLab arrow audit — section below;
(B) MuJoCo gait2392 spinal cord network in `Code\MuJoCo_SNS\spinal\` (read its
`DESIGN.md` first — NaP/laminated architecture, air walk verified 2026-09-17;
NEXT TASK = the bilateral figure repair, see ACTIVE WORK B below);
(G) NEW 2026-09-18: convert the SCALED subject01_simbody.osim on EASTEREGG2
(Ben GO) — runbook inside
`Solid_Models\OpenSim\Gait2392_Robotbody\convert_subject01.py` (10-20 min;
then `_validate_subject_model.py` must PASS before any tuning touches it);
(C) Sensory Afferent Database curation in `SADb_audit\` — section at the END of this
file, and the thread Ben wants worked while ZCode is on peak billing (Mon–Fri
23:00–03:00 Pacific).** **NEW (2026-09-16, laptop): (D) SolidWorks → Simscape Multibody
knee-rig import + native BPA muscles — see the ACTIVE WORK D section below;
(E) AnimatLab bilateral-RG walker in `Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\`
(read its `tools\SESSION_NOTES_20260916.md` first) and (F) Knee Xi-factor
extensor/flexor evaluator work — both 2026-09-14/16, sections below.** Also:
laptop-session mining findings in `Testing_Data\2022_02_Festo\HANDOFF_laptop_20260908.md`.
(For the spinal/SNS thread, `AGENTS.md` is MORE CURRENT than this file — read it too.)

## ACTIVE WORK D — SolidWorks → Simscape Multibody: knee-rig import + NATIVE BPA muscles (2026-09-16, laptop; Ben: "add to writeup for chatgpt so it can work on this")

**Laptop-only thread** (DESKTOP-5Q16KE9): repo `C:\Users\Ben\Documents\GitHub\Bipedal_Robot`,
MATLAB R2025b. Simscape Multibody license WORKS here — smimport AND hand-built
models (`e0_multibody_license.m` PASSES); the same thing is license-blocked on EB475WS4.

**Where it stands.** Ben exported the knee assembly from SolidWorks 2025 SP4.1 with the
Simscape Multibody Link add-in v7.4: `Solid_Models\Biomimetics_2022-Knee_Test\Knee
assembly\09_BA_003.xml` (730 lines) + `09_BA_003_error.txt`. The exporter flagged 18
constraints; the XML "is valid but may not reflect the original assembly". Exactly what
dropped:
- **Hinge1/2/5/6** — "constraint not supported, ignored". These are the BPA LINKAGE
  hinges: Hinge1 `04_05_BL_001-1`↔`04_02_KB_R_003-1`, Hinge2 `04_06_FL_001-1`↔KB_R_003,
  Hinge6 `04_05_BL_001-1`↔`04_01_KT_R_003-1`, Hinge5 `04_06_FL_002-1`↔KB_R_003.
- **Coincident10/13/18/19 + Parallel1/2** — "constrained components are not resolved";
  all involve the `BPA1-1`/`BPA2-1` parts mated to the assembly ROOT.
- **PathMate1/2/3/5** — patella template path mates, exported as unknown constraint.
  Expected; ignore.
- Survived: 13 parts and 24 primitive constraint pairs (Concentric/Coincident) in the
  XML `<Constraints>` section.

**Ben's rulings (2026-09-16, follow them):**
1. **Do NOT chase the BPA Coincident/Parallel drops** ("don't worry about the BPA
   constraints") and **do NOT chase the PathMates** ("the path mates are for a patella").
2. **The 4 Hinge mates should become REVOLUTE joints.** Cleanest fix: in SolidWorks,
   replace each Hinge mate with its two primitives — Concentric (pin/hole) + Coincident
   (shoulder face) — and re-export; the Multibody Link exporter translates those
   primitives into revolute joints. Alternative (if Ben prefers not to touch mates):
   hand-edit the XML to add revolute joints — requires each hinge's axis (point +
   direction) in BOTH components' frames, readable from the SW mate definition.
3. **CRITICAL correction (Ben):** the Concentric1/2 + Coincident1 set between
   `04_02_KB_R_003-1` (ground knee bracket) and `05_01_TI_R_006-1` is the **tibial HEAD
   bolted to the tibia SHANK — a RIGID group, NOT the knee joint.** Do not mislabel it
   as the knee revolute. Whether the actual knee DOF survived anywhere in the export is
   UNKNOWN until the import inventory runs.

**Next steps, in order:**
1. Import + inventory: run `import_simscape_when_ready.m` in `Code\Matlab\SNS_Simscape\`
   (inputPath already points at `09_BA_003.xml`; it imports a sanitized temp copy in the
   same folder and saves `mdl_knee_rig_import_tmp_imported.slx`), then inventory the
   model: joint blocks (ReferenceBlock under jointlib), composite bodies, and which
   frame pairs they connect. Question to answer: does ANY knee revolute exist, or is the
   tibia welded to everything?
2. Apply the hinge fix (Ben GUI, or XML surgery), re-export, re-import, re-inventory.
3. **The real goal (Ben's direction): build the BPAs NATIVELY in Simscape Multibody —
   NOT imported from CAD.** Requirements: muscles that **expand in diameter as length
   contracts** (constant-volume braid kinematics: L → braid angle θ → D(L) →
   F = (P−P_atm)·πD²/4·(3cos²θ−1), pressure as input signal) and that **contact
   complicated geometry** (bone/bracket meshes). Recipe: Force Element block or custom
   Simscape component between two attachment frames for the force law; model the
   expanding sheath as N≈10–20 short rigid segments spanning the endpoints, each with
   diameter D(L) and its own **Spatial Contact Force** block against file-solid STL
   geometry (the repo has the tibia/bracket STLs). Suppress/exclude the BPA solids from
   the SW export — their mates are the broken ones and the muscles live in the
   Multibody layer now.
4. Prototype (when Ben says go): single-segment expanding BPA — two frames, force law,
   contact against one tibia STL — as an E0b-style demo script next to
   `mujoco_bridge\matlab\e0_multibody_license.m`.

**Toolchain gotchas already paid for on 2026-09-16 (do not rediscover):** Simscape frame
ports are LConn/RConn PORT HANDLES, not 'blk/1'; find_system search-option pairs
('LookUnderMasks','FollowLinks') must precede 'Type','Block' or the search silently
returns 0; Mechanism Configuration is reachable ONLY by direct path
`sm_lib/Utilities/Mechanism Configuration` (find_system cannot enumerate it); its
gravity param is `GravityVector`; a Simscape Multibody network REQUIRES a Solver
Configuration block (`nesl_utility/Solver Configuration`) even though smimport's own
output carries none; smimport's return value is NOT the model name in R2025b (import a
sanitized temp copy; save_system renames the model to the file base). Driver scripts:
`Code\Matlab\SNS_Simscape\mujoco_bridge\matlab\laptop_step1_license_import.m` (E0 + MEX
check + import), full list in `mujoco_bridge\BRIDGE_REPORT.md` "LAPTOP PORT". The older
sw2urdf-route plan (`...\Knee assembly\09_BA_003_URDF_export_plan.md`) is superseded by
this Multibody Link route — keep for reference only.


# ACTIVE WORK E — Walker_2_Layer_CPG_BilateralRG (built 2026-09-16, easteregg2)

Ben's 8-step list: copy Connor's Walker_2_Layer_CPG, verify it runs, build bilateral RGs,
asymmetric start poses, Shinohara commissural wiring, ground test, afferents to the
half-centers, heel/toe contact drive, RG/PF/MN subnetwork pages. **Steps 1–4, 6–8 DONE
and verified; step 5 (subnetwork pages) REMAINS.** Read
`Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\tools\SESSION_NOTES_20260916.md`
FIRST — full state, layout map, complete gotcha list.

**Standing rules for this folder:** the original `..\Walker_2_Layer_CPG` is UNTOUCHABLE
(Ben's rule). Every neural change must land in BOTH `Walker_2_Layer_CPG_BilateralRG.aproj`
(Ben's GUI-visible vehicle of record) and `Walker_2_Layer_CPG_BilateralRG_Standalone_modern.asim`
(the headless test vehicle), in lockstep, sharing deterministic GUIDs (prefix `cafe…`).
Never deploy while Ben has AnimatLab open; never save from the GUI.

## Current state (all GUI-verified, zero error dialogs)

1. **Bilateral RG**: R RG ext/ext IN/flx/flx IN (clones of L's) wired ipsilateral to R
   Hip/Knee PF; the 4 crossed L-RG→R-PF connexions removed. Stimulus_1 (10 nA, 0–10 ms)
   → L RG ext and Stimulus_2 → R RG flx = antiphase kickoff. Start pose: femur_L Z=+12°,
   femur_R Z=−12°, both tibiae Y=−28° (DEGREES in .aproj, RADIANS in .asim).
2. **Shinohara commissural** (Shinohara et al. 2025, bioRxiv 2025.11.11.687930,
   Aoi/Rybak/Danner — Ben's "like the Shinohara paper"): S RG flx → `S c1` → (inh) contra
   RG flx; S RG ext → `S V3` → (exc, WEAK) contra RG ext IN. Types: "c1 Commissural
   Inhibit" (Equil −70, SynAmp 2.749), "V3 Commissural Excite" (Equil −40, SynAmp 0.1 —
   with RG-Excite strength 2.749 the V3 path LATCHES both half-centers).
3. **Afferents**: flexor Ia + NEW II chains (stretch receptor → PhysicalToNodeAdapter
   `DataTypeID=II` → "L/R Hip|Knee flx II" relay) → exc flexor HCs of PF and RG;
   extensor Ib → exc extensor HCs (24 links). Type "Afferent HC Excite" (Equil −40,
   **SynAmp 0.01** = air operating point; 0.05 slows to ~0.8 s, 0.1 latches).
4. **Contact drive**: `L/R heel contact` (foot_*_contact body) + `L/R toe contact`
   (toe_*_contact body) neurons via PhysicalToNodeAdapters (SourceDataType=ContactCount,
   Target ExternalCurrent, Gain C=20) exciting ipsilateral extensor RG/Hip PF/Knee PF/
   Hip MN/Knee MN (20 links). Inert in air, active on ground.
5. **Ground variant**: `Walker_2_Layer_CPG_BilateralRG_Ground_Standalone.asim` — Root
   unfrozen (y=−0.15) + pelvis shelf box (top y=−0.28) = partial support; plus a
   `Contact.txt` diagnostic chart. Regenerate with `tools\ground_test.pl`.
6. **Rhythm** (5.1 s headless): air 0.466 s, 11/11 bursts, clean antiphase with all
   feedback live; ground rhythm survives contact (toe duty L 16%/R 38%; walker leans
   right on the shelf — no balance layer exists, expected, out of scope until Ben asks).

## Hard-won format facts (each cost hours — do not relearn)

- **Effective synapse strength in the ASIM = the SynapseType's SynAmp; the per-connexion
  `<G>` is IGNORED by AnimatSimulator** (proved: G=1e-4 vs 0.15 → bit-identical run). In
  the APROJ the mirrored fields are the type's `<MaxSynapticConductance>` and per-Link
  `<SynapticConductance>`. Tune through the TYPE, keep all three consistent.
- **AddFlow page drawings: a drawn `<Link Org="n" Dst="m">`'s Org/Dst = 0-based index of
  the endpoint NODE entry in the CDATA's interleaved (nodes+links) file order.** All 208
  original synapse drawings match this rule exactly. Consequences: (a) appending shapes
  at the file END never shifts existing indexes (safe); (b) deleting a MID-FILE drawing
  shifts everything after it and silently re-docks all later arrows (this caused both the
  "Unable to cast MyLink to Lassalle.Flow.Node" GUI error and a 20-arrow drift);
  (c) cloning a drawing without recomputing Org/Dst renders the new arrow INVISIBLE on
  top of the template's arrow — this is the "links in tree but not on the page" bug Ben
  found 2026-09-16 evening. `tools\fix_drawings.pl` is the vehicle of record: deletes
  orphans, rebuilds all appended drawings with computed endpoints, normalizes every
  drawing, re-places nodes with collision checks; `tools\verify_handles.pl` re-audits any
  edit (final state: 296/296 endpoint-exact).
- Physical bodies (muscles, stretch receptors, foot/toe contact boxes) ARE drawn as page
  nodes; adapter links dock to them like any node.
- APROJ vs ASIM schema: `<SynapticTypeID>`/`<Text>`/`Value/Scale/Actual` attribute
  triplets/degrees vs `<SynapseTypeID>`/`<Name>`/plain values/radians. Cloned blocks must
  re-roll EVERY child `<ID>` (CaActivation/CaDeactivation/Gain) or the sim throws
  "same ID twice". Perl `s///` never writes `"".expr.""` (literal text — corrupted the
  CDATA once). Validate each page CDATA as its own XML document.
- GUI check loop: `powershell -File tools\gui_err_text.ps1 <aproj>` launches AnimatLab2,
  dumps any Error-window text, kills the app. WinForms error windows are NOT #32770
  dialogs — match by title. Headless: `bin\AnimatSimulator.exe <asim>`; charts land in
  the asim's folder. Analyzers (copies in %TEMP%, specs in SESSION_NOTES): rganalyze.pl
  (bursts/period/antiphase vs −60 mV), vmean.pl, connmap.pl (full connexion map).
- Backups per milestone: `tools\backup_v0` (pristine copy) → v1_bilateralRG →
  v2_commissural → v3_preafferent → v4_precontact → v5_predrawingfix.

## REMAINING work (next session)

1. **RG/PF/MN subnetwork pages** (Ben's step 5 — "no connection changes, just how they
   are represented"). Recipe discovered from Biped_2xCPG_wSubs.aproj (15 pages): each
   page belongs to a child subsystem `<Node>` inside `<NervousSystem>` that owns its own
   `<Links>` collection + `<DiagramXml>` CDATA page; the neurons stay in ONE flat
   `<Nodes>` list (not nested). W2L currently has a single NeuralModule Node (flat
   neurons + links + one page). Port = create 3 child subsystem nodes ("RG Layer", "PF
   Layer", "MN Layer"), partition the flat neuron Nodes + Links among them, build each
   page's CDATA from the existing drawing entries, then GUI-verify + headless re-run.
2. **Chart columns**: Rhythm_Generator.aform (+ the asim's RG chart) has RG + c1/V3
   traces only; add II relay and contact neuron columns (clone an existing DataColumn).
3. **Tuning** knobs: "Afferent HC Excite" type SynAmp (air point 0.01), contact Gain C
   (=20), commissural V3 SynAmp (0.1). Ground walking/balance is a later architecture
   task — do not start it unprompted.

# ACTIVE WORK F — Knee Xi-factor: extensor single-bracket K=[X1,X2,X2] runs + flexor two-bracket evaluator (2026-09-14/16, EB475WS4)

**Evaluator change (Ben-directed, LIVE — all future extensor runs use it):**
`minimizeExtX3.m` insertion-bracket stiffness is now `K = [X1, X2, X2]`
(was `[X2, X1, X2]`, old line kept commented in `fortz`). Single bracket — the
evaluator has NO useB2; its 6th arg is transMode (default 2trans via
`EXTX3_TRANS`). Xi3 here is the UNITLESS wrap-loss factor in [0,1] (linear in
the driver's x(4), NOT log10; the flexor X3 family's log10 series-stiffness Xi3
does not apply).

**Driver changes (`minimizeExt10mmX3.m`):** new env hooks in the Solver
section — `EXTX3_HOLD='1,8'` (comma list) = ONE custom fold with those holdouts
(skips nchoosek); `EXTX3_ALLTESTS=1` = allBPA becomes all nine tests. Empty =
original 10-fold CV. Latent loop bugs fixed (same `length()`-on-a-row family as
the flexor driver): CV loop and compile loop now iterate `size(list,1)` fold
ROWS (the old `length(list)` re-ran folds when list was a 1xN row — a
single-fold `[1 8]` list ran the fold twice and crashed), and the front index is
`ind = (1:size(x2,1)).';` (the old `1:length(x2)` + transpose broke on 1x4
single-row fronts). Plot sections reworked: each of the four sections (torque,
length, moment arm, strain) makes TWO figures — TRAINING (allBPA minus fold
holdouts) and VALIDATION (the fold holdouts), tiles subtitled simply
"Training"/"Validation". Ben's pick section (L199-216) untouched.

**Four completed runs** (all `EXTX3_PASS=2`: Xi1/Xi2 LOCKED to the pair, GA
solves Xi0 and Xi3; single bracket; full-workspace mats in
`2022_02_Festo\`; runners + log in `Dig_out\`):
- `minimizeExt10mmX3_results_20260916_pick1_h18.mat` — lock = 20260910 front
  pick-1 pair (4.354e4/1.701e4), one fold holdout {1,8}: pick Xi0 −1.60 mm,
  Xi3 0.797, mean RMSE 1.448 / FVU 0.650, filtered 18/18.
- `minimizeExt10mmX3_results_20260916_pick107_h18.mat` — lock = pick-107 pair
  (5.624e4/1.854e4), holdout {1,8}: Xi0 −1.55 mm, Xi3 0.796, 1.457/0.653.
- `minimizeExt10mmX3_results_20260916_pick1_all_h3479.mat` — pick-1 lock,
  allBPA = all 9 tests, one fold holdout {3,4,7,9}: Xi0 −1.20 mm, Xi3 0.446,
  1.195/0.533.
- `minimizeExt10mmX3_results_20260916_pick107_all_h3479.mat` — pick-107 lock,
  allBPA = all 9, holdout {3,4,7,9}: Xi0 −1.08 mm, Xi3 0.480, 1.163/0.484.
**Lock-pair provenance, spelled out:** the runner read each pair from a
DIFFERENT mat — the pick-1 pair (4.354e4/1.701e4) from
`minimizeExt10mmX3_results_20260910_noT3.mat` sol_actual, and the pick-107 pair
(5.624e4/1.854e4) from `minimizeExt10mmX3_results_20260914_pk107lock.mat`
sol_actual. Both pairs originate in the same flexor front
`minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat` (rows 1 and 107 of
filtered_results); each extensor mat preserves its pair unchanged because the
extensor never searches Xi1/Xi2 (they are its lock input, `EXTX3_PASS=2`).
**Finding:** within each case the two lock pairs converge to nearly identical
picks — with one bracket and K=[X1,X2,X2] the fitted Xi0/Xi3 barely depend on
the locked pair. all-tests cases fit better with lower Xi3 (0.45-0.48 vs 0.80).
**Ben's note (2026-09-16): despite all of the above, he will still probably USE
the two-bracket `_noT3` flexor front results —
`minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat` at pick 1 or pick 107
(the dissertation-settled line, "like we talked about"). The K=[X1,X2,X2]
single-bracket extensor runs and the flexor x3u runs are comparisons; they do
NOT replace that line unless Ben says so.**
Plots: `Dig_out\plot_x122_fourcases.m` (+ `plotblock_ext_x122.m`) re-creates
all 32 figures from the mats. Also this session:
`minimizeExt10mmX3_results_20260916_1translock.mat` (FULL 10-fold CV with
Xi1/Xi2 locked to the 1trans pair — pick Xi0 −0.91 cm, Xi3 0.842, mean RMSE
1.638/FVU 0.746), and `minimizeExt10mmX3_results_20260910_noT3.mat` was
REHYDRATED in place (labels/allBPA/numHold/baselineScores added so the pick +
plot sections run from it; pre-fix backup in `Dig_out\..._BACKUP_20260916.mat`).

**Flexor two-bracket evaluator (same sessions, companion work):**
`minimizeFlxPinX3.m` is the ORIGINAL single-bracket flexor X3 evaluator,
restored bit-exact (baseline verified). `minimizeFlxPinX3_2brkt.m` is the
two-bracket variant: bracket 2 in the tibia frame at
Pbri2 = [30.5, −103.41, 0] mm from the knee ICR with K2t = [X1, X2, X1]
(shared Xi1/Xi2), compliance-only (no path row); screw-head CLAMP (tibia-frame
X deflection pinned at −3.5 mm by contact normal force Nc, chain equilibrium
re-solved — bpa.screw2_N/screw2_hit); Xi3 = UNITLESS wrap-loss
delta_L = Xi3 · 15 mm · theta_wrap · comp² (theta_wrap = pi minus the angle
between the class force direction and the bracket2→insertion line — straight
through = no wrap); 30 mm-circle tangency check with env
`FLXPX3_TANGENCY=DIAG` to store-but-not-enforce (with today's 2-row paths the
check rejects everything — needs the wrapped ≥3-row path from Ben's CAD);
+5.3° 47 cm encoder shift ON by default (`FLXPX3_NOSHIFT=1` = legacy).
Driver variant `minimizeFlxPin10mmX3_2brkt.m` (calls the 2brkt evaluator;
`FX3B_LOCK1/2` env pins Xi1/Xi2 — GA solves Xi0/Xi3; its compile-loop
`ind = 1:length(x2)` bug FIXED — the SAME bug is still LATENT in
`minimizeFlxPin10mmX3.m`, untouched). Flexor x3u lock-CV mats
(20260915, DIAG): `minimizeFlxPin10mmX3_2brkt_results_20260915_L107_DIAG_x3u.mat`
= pick Xi0 ≈ 0.1 mm / Xi3 0.150 / mean RMSE 1.587 / FVU 0.077 (best
pinned-flexor fit on record) and `_L1trans_DIAG_x3u.mat` (Xi3 0.291). The
same-named 20260915 mats WITHOUT `_x3u` used the old series-stiffness Xi3 and
are superseded.

**Gotchas hit (do not relearn):** (1) pool workers never see the client's
`setenv` after spawn — env-gated evaluator branches (FLXPX3_TANGENCY) silently
run the wrong mode inside parfor; set env BEFORE parpool or evaluate
client-side. (2) `for k = 1:length(list)` on a 1xN row list re-runs folds;
iterate `size(list,1)`. (3) A 1x4 single-row GA front breaks
`ind = 1:length(x2)` row assembly — use `(1:size(x2,1)).'` (fixed in the two
drivers named above; still latent in `minimizeFlxPin10mmX3.m` and
`minimizeExt10mmX3.m`'s siblings). (4) Do not relax the strain < −0.03 NaN
cutoff or fake its values in Go_OfF — the stretch-side festo4 continuation is
unvalidated (manufacturer 2-3% limit; stretch characterization never done —
antagonist pairs + fatigue); curve truncation past the measured domain is
INTENTIONAL (Ben ruled 2026-09-16).

## ACTIVE WORK B — MuJoCo gait2392 spinal network (started 2026-09-09; state below verified 2026-09-17, EB475WS4)

Two-level RG+PF spinal network (SNS-Toolbox 1.5.2) driving all 92 muscles of
the converted `gait2392_simbody` MJCF, now on the persistent-Na (NaP) RG +
laminated-PF architecture with curriculum tuning, mechanosensory heel/toe +
stance-gated Ib/II feedback, IaIN + Renshaw, conditional KINH swing
suppression. **Air walk VERIFIED** (`_smoke_nap.py`, full MuJoCo–SNS loop at
2 ms: 5 RG bursts/leg in the 5--15 s window, period 2.11 s, E-duty 0.67,
knee −106..+12°, hip −20..+51°; `spinal_run.npz/.png` = that run). Ground
walk awaits a valid NaP stage-3 curriculum winner.

**NEXT SESSION'S FIRST TASK — the bilateral architecture figure repair.**
Ben's annotated screenshot found four presentation defects in
`circuit_literature` panel A: right E/F columns not mirrored about the
midline, PF-box connections that look unattached, PF interneurons that look
self-exciting, and sensory-feedback lines that look unanchored. The full
defect list and the 7-step SAFE implementation/verification sequence are at
the END of `CHATGPT_REPORT.md` (section "2026-09-17 — Bilateral architecture
visual-audit stop point"). Render a temporary review PNG first; do NOT touch
the promoted figure or its dissertation copies until the review render
passes all four annotated regions.

**Changed-design facts every future edit must respect (Ben-corrected
2026-09-17):**
- Commissurals are FOUR directional cells: RG-E_l -> V3_l-to-r -> InE_r
  (both synapses excitatory, mirrored r-to-l) and RG-F_l -> C1_l-to-r -|
  RG-F_r (mirrored). Never collapse them into one shared V3 or C1 — the
  compiled network carries CIN_E_r/l and CIN_F_r/l. The representative
  figure builds `rg_weak_exc=0`: no direct ipsilateral RG-E<->RG-F
  excitation.
- The E1/E2/F1/F2 PF cells are PHASE CHANNELS (effective temporal rank ~2;
  corr(E1,E2)=0.998), NOT validated synergies. NMF synergies are labeled
  S1–S6 only. Six is the smallest shared bilateral PF count clearing 90%
  held-out VAF (0.924 R / 0.912 L) on the converted-model activations — an
  engineering target for this ~1.6-cycle record, not a biology claim; one
  right sixth component is spatially unstable.
- FSA backsolve (`fsa_backsolve.py` -> `fsa_results\`): the
  excitatory-only dynamic PF->MN fit leaves ~16–17% of MN samples requiring
  net negative current, so closure needs phase-specific inhibition, not
  just membrane leak.
- Activation target of record = `bsolve_out.npz['acts']` (converted-MuJoCo
  ridge/NNLS from `bsolve_ik.py`). The surviving
  `ResultsBSolve\zz_bsolve_*_activation.sto` is RULED OUT (RMSE 0.39, corr
  −0.01). The `normal.mot` OpenSim SO experiment is QUARANTINED in
  `ResultsNormalSO\` (pelvis residuals carry ~739 N vertical — not a
  walking target; never use it to tune PF counts or synapses).
- Ankle dorsiflexors are NOT underpowered: summed DF torque 1.023× OpenSim
  in the IK gait range (`compare_ankle_df.py`, `ankle_df_results\`). Any DF
  deficit is recruitment timing/magnitude, co-contraction, or
  force-velocity — NOT Fmax or static moment arms. Use equality-aware FD
  tendon moments; raw `actuator_moment` is exactly 0 for the converted
  ankle (and knee) paths.
- KINH is a conditional swing-phase knee-extensor inhibitory IN
  (`f1_kneext_inh`, default 0 = absent). Always label it conditional; it is
  not a V-class or literature-named population.
- Figures: extensor = blue, flexor = vermillion/orange (literature
  convention). `draw_literature_circuit.py` asserts every displayed solid
  neural edge against a freshly compiled network; `draw_circuit.py`'s
  compiled-edge contract is two-way (missing AND extra edge groups fail).

**Env + tools:** on EB475WS4 run
`C:\Users\Ben Bolen\.conda\envs\myo\python.exe` (py3.10.21, numpy 1.22.4,
scipy 1.9.3, mujoco 2.3.7, SNS-Toolbox 1.5.2), cwd `Code\MuJoCo_SNS\spinal`.
(`D:\Anaconda\envs\myo` is the easteregg2 twin — use it only there.) Gates
that must stay green after ANY change: `_fix_check.py`, `_panels_check.py`,
`audit_signs.py` (DYNAMIC sign audit), `_live_plant_audit.py`. Analysis
chain: `bsolve_ik.py`, `fsa_backsolve.py`, `_fsa_rank_robustness.py`,
`_pf_basis_audit.py`, `_synergy_count_audit.py`,
`_activation_provenance_audit.py`, `_normal_so_audit.py`,
`compare_ankle_df.py`, `gait_phase.py`, `plot_gait_joint_angles.py`,
`_smoke_nap.py`. `runner.py` summaries now resolve channels by
`NEURO_NAMES`, plot degrees, and restrict rhythm metrics to the 5--15 s walk
window; `_figure_hindlimb_style.py` renders beside the target and replaces
atomically (Windows preview locks).

**Open after the figure repair:** (1) re-run curriculum stages 2–3 — their
current files hold −100 no-countable-cycles sentinels (decide retain vs
purge of those studies first; stage-1 winner 137.254 / trial 78 is valid);
(2) the supported-ground figure and `opensim_overlay_gait_cycles.png` wait
for a valid NaP stage-3 winner; (3) the six-synergy PF-layer decomposition
(hip/knee/ankle × ext/flex half-centers, convergent drive for biarticular
MN pools) is DESIGNED but not implemented — see the report's
"Phase-normalized synergies, joint kinematics, and ankle capacity" section.

**Standing model facts (unchanged — do not rediscover):** repairs live in
`apply_harness()`: hip_flexion/hip_adduction hinge axes negated vs OpenSim;
rect_fem rerouted over the vastii vas_med-P4 patella point;
quad_fem/gem/peri pruned; 8 trunk-Fmax repairs to stock, verified live
(ercspn 2500 N, intobl/extobl 900 N, ext_hal 162 N; 78/86 actuators within
15% of stock, worst gluteal ~31–34% — the draft now says exactly that).
**Do NOT weld the conditional pathpoints** (equality couplings +
`boundmass` 0.1; welding froze moment arms). The sign audit must run
DYNAMICALLY — static `actuator_moment` misses the coupler paths. MuJoCo
muscle Fmax = `actuator_gainprm[:,2]`. Full detail: `spinal\DESIGN.md` top
(2026-09-17 sections are authoritative) and the five 2026-09-17 sections of
`CHATGPT_REPORT.md`.

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
   **Update 2026-09-16 (ZCode session)**: prime suspect remains the Na h-gate tau
   (AnimatLab LinearHill treats tau_h.max as a FIXED constant, which is what makes
   the verified Simulink Deng port oscillate; a voltage-dependent tau_h that
   collapses at depolarized V quenches this circuit — see
   `Code\Matlab\SNS_Simscape\README_SNS_Simscape.md`, Deng section). Also note the
   python side's toolbox DOES ship a persistent-Na neuron class
   (`NonSpikingNeuronWithPersistentSodiumChannel`, Tutorial 8) — if the MuJoCo-side
   spinal RG is rebuilt on it, all three platforms can share one Deng-style
   burst-termination mechanism (mind the ThrPre/Elo opposite-saturation trap when
   porting values).
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

## Xi values for the dissertation text (SETTLED 2026-09-10, Ben-approved picks — verified by direct .mat loads)

- Flexor, pinned 2brk: minimizeFlxPin10_results_20260908_2brkt_2trans_noT3.mat, pick 107
  (Xi0 = +8.9 mm, Xi1 = 5.62e4 N/m, Xi2 = 1.85e4 N/m; per-test RMSE 1.98/1.36/2.22/1.47/1.23,
  all FVU < 0.16; 1trans twin = ..._1trans_noT3.mat, Xi0 +5.95 mm, Xi1 4.36e5, Xi2 2e4)
- Extensor, pinned: minimizeExt10mmX3_results_20260910_noT3.mat, pick 1
  (Xi0 = -10.1 mm, Xi3 = 0.621; Xi1/Xi2 locked to the flexor pair 4.354e4/1.701e4, never searched;
  held-out RMSE 1.04-1.64 vs baselines 2.44-3.62; bio-ext 52cm 2.195/2.350/4.727 vs baseline
  3.032/4.486/5.587 — FVU 2.35 > 1 there, do not overclaim)
- Signs matter: flexor Xi0 positive, extensor Xi0 negative. Lower-bolt-hole Pbr
  [-6.26,-29.69,75.06] tested equivalent (Xi3 0.60 vs 0.62). Mesh optimizations
  Opt_run/Opt_run_Ext ran with these values (Opt_run_Ext results:
  Vas_Pam_20mm_Result_20260910_0528.mat). Everything above verified by direct .mat loads.

## Rules for the dissertation-text chat (2026-09-10)

1. The dissertation lives in `Documentation\Reports and Papers\Dissertation\ProofFinal\`
   (LaTeX; main.tex + chapters). Compile via its own setup if asked; do not restructure.
2. Xi values in the text must match the SETTLED block above exactly, including signs
   (flexor Xi0 positive, extensor Xi0 negative) and units as written in the text.
3. If Ben asks about fit quality, the honest summary is: flexor all-5 FVU < 0.16;
   extensor held-out ~2x better than baseline; bio-ext better than baseline but FVU 2.35 > 1.
4. Do not run optimizers. Values are settled; text updates only.

## DISSERTATION CPG SECTION — DECONFLICTION (2026-09-20, ZCode on EB475WS4). READ FIRST if editing the CPG/spinal material

Another session (ZCode, EB475WS4) filled the curriculum result slots in
`Documentation\Reports and Papers\Dissertation\CPG_spinal_section_draft.tex`
on 2026-09-20. If you are editing that file (or copying its blocks into
ProofFinal chapters), FIRST `git pull` (Ben pushes from EB475WS4 via GitHub
Desktop), then grep the file for `ZCODE 2026-09-20` — every block added or
rewritten by that session is fenced with `% >>> [ZCODE 2026-09-20]` ...
`% <<< [ZCODE 2026-09-20]` (9 pairs). Rebase your edits around those
fences; do not overwrite them blind.

What changed there: (1) Methods/architecture: conditional-population count
clarified (510 neurons / 382 inputs with stage-2/3 sensory populations
active); (2) Methods/tuning: rewritten — full-dict seeding, real trial
counts, the rhythm-gate and small-gain-entry objective rules, and the new
contact-event drive (heel-strike/toe-off transients on the heel/toe ports,
Conway 1987 reset-to-extension); (3) Results/curriculum: the old "corrected
campaigns remain to be run" paragraph replaced by the static-pose-exploit
diagnosis + stage-2b (seed wins, afferents tolerated not additive) and
stage-3b (first supported-ground gaits, 11/60 trials) outcomes; (4)
Results/ground-walking: the \fillme stub FILLED with the stage-3b winner
(trial 52, score −85.22): 20 counted cycles, duty 0.69 vs 0.61, knee
−53.5° cycle / −76..+10° raw, hip excursion 16.9°, ankle 34°, tilt 29.7°,
COM 0.80, 0.12 m progress; (5) Results overlay figure repointed to
`curr3b_ground_cycles.png` + two new figures (`curr3b_ground_traces.png`,
`curr3b_ground_limbs.png` — generated from the winner run by
`Code\MuJoCo_SNS\spinal\plot_run.py`; old v9 overlay file kept untouched);
(6) Discussion: new subsection "Scalar Tuning Converges; the Remaining
Gaps Are Architectural" + the Limitations \fillme filled. No .bib changes
(all cites use existing keys). PART D bib block untouched.

Provenance if numbers are questioned: full session log =
`Code\MuJoCo_SNS\spinal\DESIGN.md` section "2026-09-20"; winner params =
`curriculum_stage3.json` (trial 52) and `curriculum_stage2.json`; rerun any
figure via `_diag_stage3.py` then `plot_run.py` (env + steps:
`Code\MuJoCo_SNS\HowToRunCode.md`). Files this session created/changed
beyond the tex: the three `curr3b_ground_*.png` figures, `HowToRunCode.md`,
`_curriculum.py`/`runner.py`/`params.py` (rhythm gate, contact_onset),
DESIGN.md 2026-09-20 section, and this note.

**AMENDMENT (later the same evening, ZCode):** Ben's figure audit showed the
stage-3b "winner" was a ONE-LEGGED gait (left foot planted 100% of frames —
the v1 objective scored the right leg only, removed joint means, and took
duty from the neural rhythm). The kinematic objective was REBUILT
(kine_ref v2: both legs, mean/amplitude/phase/periodicity vs each leg's own
OpenSim reference, cycles from per-foot contact onsets, frozen-leg guards;
runner now logs per-foot contact/foot height every ground step and takes
`AARL_PELVIS_TY`). A fresh 80-trial ground study (curr_s3c) ran under it;
corrected re-ranking shows its ENTIRE front is the frozen-left family
(best −181.6, trial 54). The tex fences were REWRITTEN to this corrected
story (grep `ZCODE 2026-09-20`; the ground/results/discussion blocks now
carry the evening version), the curr3b_* figures were DELETED from
CPG_airstepping_figs, and `curr3c_gait_isb.png` / `curr3c_gait_contact.png`
(t54, OpenSim/ISB conventions) are the current figure set. If you already
copied the earlier (pre-correction) numbers anywhere, replace them from the
current fences only. Details: DESIGN.md "2026-09-20 EVENING" + "s3c
OUTCOME" blocks.

## Dissertation text session 2026-09-10 (ZCode, easteregg2) — DONE, local only, NOT uploaded

Synced Ben's 11:50 AM Overleaf zip into ProofFinal (advisor edits preserved: abstract wording +
removed "Planned extensions..." sentence; intro "Muscle Mutt" rename + gap-statement rewrites;
BolenFrontiers22.bib bolen_2026 -> Actuators entry; thesis.bib/BolenFrontiers22.bib restored at
ProofFinal root, byte-identical to HEAD). Then edited 8 chapter files (all under
Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\):

- 94-AppendixC: flexor + extensor test tables now carry ACTUAL l0/tendon/P/Fm read from
  FlxPinBPASet.mat / ExtPinBPASet.mat (extensor shorthand was wrong: "42cm"=0.415, "43cm"=0.436/
  0.432, "46cm"=0.457, "47cm"=0.465, "48cm"=0.480; ke(6) tendon 0.022); Fm filled for all rows;
  pool/excluded roles marked. Adopted-values table rebuilt with the SETTLED picks (flexor row 107:
  +8.9mm/5.62e4/1.85e4; extensor pick 1: -10.1mm/0.621 on locked pair 4.354e4/1.701e4 = flexor
  front row 1, same 20260908_2trans_noT3 mat). New bounds table (flexor X1 [3e4,1e6], X2 [5e3,2e4]).
  Frame construction + compliance chains converted from inline math to numbered equations.
  All provenance line numbers re-verified against the working tree (subagent audit, 59 citations):
  corrected 2brk (frames 252-298, K 387/396-411, UD 536-547, fzero 439-446, encoder corr 66-70,
  labels line 51), ExtX3 (frames 461-483, K 575-584, UD 687-700), ext seed route 380-389, flexor
  seed 178-190, contexts 156-174/158-165; crossPredictFlx.m -> Dig_crossPredict.m line 46,
  sweepExtX3.m -> Collect_ExtPinX3_sweep.m lines 59-60; minimizeExt.m line 352 K-array is now
  [X1,X2,X1] (changed per Ben 2026-09-10; table updated).
- 92-AppendixA: l_m formula -> line 254 (not 368-377); strain discard -> line 308; hypot ->
  predictKneeFlexor20mm.m line 199; MonoPam_mult coupled balance -> 490-564; F* zero-clamp
  correctly attributed to festo4.m:35; "Note for review" RESOLVED (two-bracket fit of record uses
  diag(X1,X2,X1) matching eq:bktstiffness; legacy [X1,X2,X2] marked superseded; new display eq).
- 93-AppendixB: stale TODO block removed; file-of-record updated (fresh _Standalone.asim export
  2026-09-09 runs headless; .aproj authoritative); SNS-Toolbox section now describes the built
  406-neuron spinal network + verified antiphase rhythm.
- 20-methods: new paragraph in sec:improved_model describing the two-bracket flexor extension +
  locked-pair extensor refit (points to App C tables + sec:ongoing). Three preliminary-sim
  subsections updated with current status (AnimatLab bilateral reorg + RG latch; MuJoCo conversion
  verified + 0.1ms coupling + spinal net + NaN blocker, foot meshes first suspect; Simulink base-
  Simulink operation + URDF route verified + placeholders still open).
- 30-results: "Ongoing Experiments" -> "Follow-Up Identification and Route Redesign" (label
  sec:ongoing kept): settled values + verified GoF; 20mm bio-flexor l0 41.5->42.0 cm (mat says
  0.420); bio-ext 52.0->51.8 cm everywhere; extensor redesign meets target at all angles (thinnest
  margin <0.1%); flexor re-run flagged as IN PROGRESS (no 20260910 flexor result mat exists).
- 40-discussion: 48->48.5 cm for the fit BPA; "42 mm" typo -> 41.5 cm test; run-on fixed.
- 50-futurework: pipeline status paragraph (assembly done, NaN blocker remains) + tonic-test
  status sentence in Verification Tests.
- 60-conclusion: parenthetical noting conversion/coupling/spinal-net/rhythm already in place.

DATA-CHECK DISCREPANCY (needs Ben's eyes): the SETTLED block quotes flexor pick-107 per-test RMSE
1.98/1.36/2.22/1.47/1.23 and "all FVU < 0.16". Direct evaluation with the current working-tree
minimizeFlxPin2brk.m at (0.0088574, 56238, 18542) 2trans gives RMSE 1.94/1.51/2.29/1.63/1.21 and
FVU 0.14/0.06/0.17/0.04/0.03 (test 3 = 0.1651, i.e. NOT <0.16). The DISSERTATION uses the
reproducible numbers. Xi values themselves confirmed identical to the settled block (mat row 107
= filtered_results(107,4:6)). Extensor pick + pool GoF + bio-ext 2.195/2.350/4.727 all confirmed
(own evaluation + Dig_ExtPin_frontBio_newXiPair_20260910.log).

Also deleted the tracked-but-junk _zipreview/ folder (unstaged deletions; include in next commit).
No commits made. Ben uploads the 8 changed chapters/ files to Overleaf (built from his 11:50 zip,
so advisor text edits are preserved).

# ACTIVE WORK C — Sensory Afferent Database (SADb), `SADb_audit\` — THE THREAD FOR PEAK-HOUR COVERAGE

Ben's request (2026-09-16): work this thread while the primary assistant (ZCode) is on
peak billing (Mon–Fri 23:00–03:00 Pacific). The corpus is Ben's Sensory Afferent
Database — ~940 locomotor sensory-feedback papers in the Airtable base "Sensory
Feedback", curated against his Zotero libraries. **Read `SADb_audit\README.md` (the
CURATION SPEC section is the detailed rulebook) and the WORKFLOW rows of
`SADb_audit\curation_log.csv` (the running state) before doing anything.**

## Measured state (2026-09-16, direct Airtable counts + logs)

- Papers table = **943 records**, **500 with empty Notes**. A Sept-15/16 "task5"
  auto-curation pass (`SADb_audit\task5_progress\`) created ~390 new records and
  auto-curated many; its `task5_insufficient.json` lists **108 Zotero keys** whose
  auto-curation was judged insufficient (manual curation needed). Before curating ANY
  paper, reconcile `task5_classified.json` / `task5_created.json` against the live
  table — task5 coverage overlaps the older batch queue in ways not yet fully mapped.
- The original 383-record rest-import campaign: batches 1–5 DONE (50 papers) +
  independent audit PASS; Ben's 2026-09-15 rulings applied (see the curation_log
  RULINGS row). Batches 6–7 DONE (2026-09-18/20, 16 more papers, 4 logged no-text
  chapters) — **batch 8 = queue CSV rows 62–70** (grounding already fetched in
  `batch6\ground_*.txt`; running total 66/383; audit due after batch 10).
- PDFs: 318/886 DOI-bearing records carry attachments (2026-09-20). THE definitive
  no-PDF hunt list is `author_fix\remaining_no_pdf.csv` (581 records, columns
  doi/record_id/title/hunt_status; rebuilt by `rebuild_no_pdf_list.py`; every DOI
  in the corpus has been hunted once — 488 have no OA copy anywhere, 88 have
  candidate URLs that fail %PDF verification). The Sept-15 intermediate lists were
  deleted 2026-09-20. Zotero-local PDF inventory: `pdf_inventory.csv` +
  `author_fix\pdfs_from_zotero.csv`; staged copies live at
  `D:\sadb_pdf_staging` (OUTSIDE the repo on purpose — never move them in).
- Author normalization DONE table-wide (2026-09-15): Primary Author = ONE surname
  (no "et al."), remaining surnames in the new `Secondary Authors` multi-select.
- VOSviewer citation bubble map BUILT: `SADb_audit\vosviewer\` (`sadb_map.txt` +
  `sadb_network.txt`; open at app.vosviewer.com; regenerate with `vos_build.py`).

## Work menu (Ben assigns; suggested order)

1. **Curation batches** (10 papers per batch) for whichever records are still bare.
   Follow the README spec EXACTLY: grounding ladder (Zotero → PubMed → Europe PMC →
   publisher/archived full text), field rules (Notes style, Animal vocabulary,
   Feedback link vocabulary), Review Papers / Models twin-record creation BEFORE the
   paper update, ONE batched update per 10 papers, verify by re-pulling those ids and
   diffing every field echo, append one row per paper to `curation_log.csv`. Audit
   every 5th batch (read-only re-verification of the last 50 records).
2. **Task-5 manual curation**: the 108 `task5_insufficient.json` keys. These failed
   auto-curation because grounding was thin — they need the full treatment and the
   no-grounding-no-note rule is MOST important here.
3. **PDF hunt** on `remaining_no_pdf.csv`: find OpenAlex/publisher OA candidates,
   verify each URL with a range-GET checking the `%PDF` magic BEFORE attaching (an
   earlier session's unverified URL attaches were silently dropped by Airtable), then
   URL-attach. Local files cannot be uploaded by token — they go in via Ben's
   drag-drop or the tunnel bridge (`author_fix\tunnel_attach.py`, Ben runs it).
4. **Viz/app builds** if Ben asks: VOSviewer map refresh after batches land; a
   single-file HTML pivot/search/bubble app fed by an Airtable export. New artifacts
   stay small and text-only under `SADb_audit\` — the repo is in a size-reduction
   campaign, so NO PDFs or binaries into the repo, ever.

## Credentials + rotation (read them, never copy them)

- Keys live in **`D:\Github\api_credentials_local.txt`**, outside the repo. Read
  them from that file; NEVER write a key into any repo file, a chat message, or
  CHATGPT_REPORT.md. Airtable REST: `https://api.airtable.com/v0/appMQTnobUNRytIp7/<Table>`
  with `Authorization: Bearer <PAT>`; working script patterns in
  `SADb_audit\author_fix\*.py` and `task5_progress\auto_curate.py` (reads env AT_PAT).
- **Rotation is PENDING**: both keys appeared in AI chat transcripts. The plan is a
  fine-grained Airtable PAT scoped to this base only (data.records read+write) and a
  read-only Zotero key, with the old keys revoked after all processes migrate. If Ben
  has rotated by the time you read this, use the new keys from the file.
- **Zotero is STRICTLY READ-ONLY** (Ben's standing note in the credentials file).
  The web-API DELETE bypasses the Zotero trash and is PERMANENT — one 107-attachment
  purge already happened, on Ben's explicit per-item criterion. No Zotero writes,
  deletes, or file uploads without Ben's explicit per-action go.

## Hard rules (each prevents a failure that already cost real repair time)

- NEVER delete an Airtable record; never edit the 99 original papers' existing notes;
  never touch Recorders; never rename Feedback/Models records; no new Animal options;
  no proxied (`proxy.lib.pdx.edu`) URLs anywhere; the DOI is the canonical identifier.
- Grounding is MANDATORY. Notes are written FROM the paper's text (abstract minimum),
  never from model memory. **Output contract: every field value carries a verbatim
  quote + source locator, or the field stays empty and gets flagged. No text found →
  NO note, log row says `no-text`. Never guess; ambiguity → flag for Ben in
  `curation_log.csv`.**
- Grid views and interfaces are Airtable UI objects — NOT creatable via the API. When
  Ben asks about views (search/sort/filter/group-by are per-view, saved in the UI),
  give him click-paths; don't attempt API calls for it.
- Airtable attachment sizes populate asynchronously — never re-attach because an
  early size check read 0.
- Do NOT touch threads (0)/(A)/(B) (MATLAB, MuJoCo, AnimatLab) or any git operation.

## Machine traps that bit earlier sessions (details in AGENTS.md + README.md)

- This box's shell is cmd, not bash. PowerShell 5.1 misreads non-ASCII literals in
  UTF-8-WITHOUT-BOM .ps1 files (use `[char]0xNNNN` codepoint literals). Author JSON
  payloads from `-Raw` file reads, never console echoes (cp437 mojibake corrupted a
  CSV once). Inline multi-line `python -c "..."` gets eaten by this shell — write
  script files. Zotero LOCAL API (port 23119) is read-only and needs `&qmode=everything`
  for `q=` to match DOIs. Airtable REST writes are UTF-8-safe.

## When done

Append a dated section to `CHATGPT_REPORT.md` (never overwrite old sections): batches
done with record ids, flags awaiting Ben, counts before/after, files created. The real
state file is `SADb_audit\curation_log.csv` — keep it current and any later session
(ZCode or ChatGPT) picks up cleanly. Leave `AGENTS.md` and this handoff file to ZCode
sessions.
