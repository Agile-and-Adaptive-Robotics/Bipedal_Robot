# Spinal network for Gait2392 in MuJoCo (SNS-Toolbox) — design & status

Built 2026-09-09. Files in `Code\MuJoCo_SNS\spinal\`.

## Goal

Two-level spinal cord network (McCrea–Rybak RG + PF) with proprioceptors for
each of Gait2392's 92 muscles, driving the MyoConverter MJCF of
`gait2392_simbody` in MuJoCo: standing → walking → standing, speed via
descending drive AND reflex-gain modulation, with PF→MN weights back-solvable
from OpenSim IK/SO activation patterns.

## Architecture (per side; ~406 neurons total, 368 inputs)

| layer | cells | notes |
|---|---|---|
| descending | DRIVE, POSTURE, BAL_PF, BAL_DF | MLR/postural surrogates, external inputs |
| RG | RG-E, RG-F + ADAP-E, ADAP-F | half-center: mutual inhib. + slow (0.9 s) self-adaptation; frequency rises with DRIVE; stance-biased drive split (E stronger) sets duty |
| PF | PF-E1, PF-E2, PF-F1, PF-F2 + PFA-* | forced mode (no own rhythm); (tau, adapt) shape multipliers stagger their windows; reciprocal inhib. on conflicting pairs |
| MN+afferents | per muscle: MN, Ia, II, Ib | Ia ~ tendon velocity (pure-signal: no resting tone — baseline drove constant reciprocal inhibition and crushed flexors); II ~ length (small baseline); Ib ~ force |
| Ib load sharing | IBEXC per stance group | Ib afferents → group IN, gated by RG-E → excites that group's MNs (stance reflex reversal / load sharing) |
| balance | BAL_PF/BAL_DF | pelvis-COM PD → ankle PF/DF + hip flex/ext MNs (ankle+hip strategy), fades with DRIVE |

Standing uses a **solved posture pattern**: static optimization at the
keyframe (torque rows = gravity − contact − passive, plus co-contraction
preload; activations regularized toward a physiological prior) injected as
per-MN bias currents (`POST_<muscle>`). This is the same machinery as the
planned IK/SO back-solve, applied to standing.

Ia/II/Ib gains are speed-modulated presynaptically (`params.MOD`, scaled by
DRIVE) and phase-gated from RG potentials — the "reflex modulation = speed
control" knob.

## Literature grounding (verified; mostly in Ben's Zotero)

- **RG+PF two-level**: Rybak et al. 2006 J Physiol (PF = multiple
  populations exciting synergist MN pools, forced by RG, low persistent Na);
  McCrea & Rybak 2008; Rybak et al. 2015 eNeuro review. The 2006 model
  explicitly anticipates multiple phase-shifted PF units for bifunctional
  muscles → Ben's "similar but phase-shifted PF" observation.
- **Speed**: increased drive shortens stance, swing ~constant (Danner 2017,
  Rybak 2024 eLife RP98841 — both RG-only, no PF). **Reflex gain modulation
  as speed control**: Bunz, Ijspeert, Schmitt 2026 Sci Rep
  10.1038/s41598-026-48509-z (offline+online, walk↔run); Ross/Rybak-lab mouse
  spindle-gain work (2018, in Zotero).
- **Precedent for exactly this stack**: SNS-Toolbox paper (Nourse et al.
  2023, Biomimetic 8:247) MuJoCo case study = rat hindlimbs, two-layer RG+PF
  CPG (Deng-style), PF shapes joint commands, MNs adjusted by Ia/Ib, speed
  set by RG inhibition level. Our network is the human-92-muscle analogue.
- **Synergies**: Ivanenko et al. 2004 (5 modules, speed-invariant);
  Kibushi et al. 2018 (speed carried by timing/amplitude); Avaltroni 2024
  spinal maps. → back-solve design in `fit_synapses.py`.
- **Closest competitor**: Di Russo, Ijspeert, Bouri 2023 J Neural Eng
  (CPG+reflex human neuromechanical model) — read before claiming novelty.
- Also in Zotero Sensory Afferent Database: Markin et al. 2016
  (large-scale PF→MN neuromechanical cat model), load-receptor review 2000,
  stick-insect multisensory control, air-stepping by tonic drive (2009),
  human balance NN papers (2017/2000).

## Verified status (this session)

- `check_rhythm.py`: half-centers oscillate; DRIVE=4 nA → period 0.844 s,
  E duty 0.29-0.33, left–right antiphase corr −0.74; PF windows alternate
  (E cluster ~0.08, F cluster ~0.72 of the cycle).
- Network builds/compiles/steps; runner closes the loop on the full
  92-muscle model with logging/summary/NaN forensics.

## BREAKTHROUGH 2026-09-10 (machine EB475WS4 session): the NaN blocker is DEAD

Root causes, in the order they were found (all fixed in `runner.patch_xml`
+ `apply_harness`; every fix is asserted at build time):

1. **boundmass 0.01 is too small at human muscle forces.** With 92 muscles
   pulling on the equality-coupled pathpoint DoFs the mass matrix goes
   singular ("Inertia matrix too close to singular" → NaN) at as little as
   0.3 whole-body co-activation. → `boundmass=0.1` (survives 1.0
   co-contraction for 6 s in diag_stab.py).
2. **The equality-coupled "conditional pathpoint" slide joints carried
   range limits sampled only around the keyframe** (MyoConverter artifact).
   The equality drives them exactly along their polycoef curves, which
   leave those boxes within 10-30 deg of knee flexion (14 of 36 coupled
   dofs outside their range at knee=10 deg, 31 at 30 deg) → range-vs-
   equality fight = huge constraint forces + jammed joints. → strip the
   ranges, `limited="false"` (NOTE: default class sets limited="true"
   explicitly, so dropping `range` alone errors out at compile) and add
   `armature="0.5"` to the followers (0.1 was not enough: rect_fem P3_y
   went singular at t=8.9 s with the CPG co-contraction).
3. **knee_angle driver range: NO FLIP — Ben was right (2026-09-10 PM).**
   The converter PRESERVED OpenSim's flexion-negative convention: real-
   actuator tests (unjammed followers, `_muscle_direction_test.py`) show
   gravity buckling the standing knee NEGATIVE (flexion), semimem/
   bifemsh/med_gas driving NEGATIVE (to -77 deg), vas_lat/rect_fem
   driving POSITIVE (extension). The converted range [-2.094, 0.1745]
   (= -120 deg flexion / +10 deg extension) is CORRECT as shipped.
   Moreover knee_angle ships `limited="false"` — the range is INERT, so
   an earlier "range flip" (2g, applied and then REVERTED the same
   evening) changed nothing dynamically. CAUTION: a kinematics sweep that
   rotates knee_angle WITHOUT letting the coupled translation dofs follow
   (the tibia's rolling center) mirrors the apparent foot path and
   "proves" the wrong convention — don't repeat that. Related: the
   audit_signs moment-injection method reports moment-arm signs that are
   SIGN-INVERTED vs real activation (docstring caveat added); the hip
   "axis flips" from the previous session were physical no-ops (axis
   negation only relabels qpos signs).
4. **collision="predefined" + 19 explicit ground <pair> lines**: predefined
   pairs IGNORE contype/conaffinity, so muting the ground geom does
   nothing (the first "air" runs dragged on an invisible floor). → patch
   surgery on the <contact> section: remove it entirely for suspended
   tests; keep only the 6 foot pairs on the ground (OpenSim gait2392 has
   no self-collision either).
5. **The pelvis rig springs were landing** (they always had been in the
   original code — a mid-session refactor briefly broke them, caught by
   _rig_check). Rig now also writes its damping (was accepted but never
   written before: undamped rig spring) and the 3 pelvis ROTATION hinges
   get springs too (150/25 ground, 400/40 air): with zero yaw damping the
   pelvis freely spun about vertical at 40-68 rad/s on the ground
   (foot-friction spin-up from left-right standing-solve asymmetries).
6. **Standing solve fixes**: always solve on a GROUND-ON model copy (in
   air the solve degenerates: nothing carries gravity → all activations
   <=0.2, wrong muscles); fit only DoF rows with nonzero muscle moment
   (pelvis translation rows carry the full body weight but no muscle can
   act on them → the ridge smeared everything to mush); stronger
   co-contraction preloads (hip -40, knee -60, ankle +40 N·m), lam 0.25.

Results after the fix stack (myo env on EB475WS4):
- **Deafferented + suspended in air: 22 s, zero warnings, 11 alternating
  RG_F bursts, COM stable on the rig.** knee range still only ~-1..1 deg
  in air: quads (9.2 kN total Fmax) vs hamstrings (5.5 kN) co-contract —
  at PF-cell current levels the net knee torque never crosses zero. The
  tonic DRIVE->PF term was halved (drive_to_pf 0.4→0.2) and PF reciprocal
  inhibition raised (2.0→3.0) to narrow the windows; W_PF_MN retuned
  (F1 knee_flex 1.2, F2 knee_ext 0, E1 knee_ext 0.25/ankle_df 0.10).
- **Ground walk with afferents (default runner config): 22 s, stayed up,
  11 bursts, contacts 3-11.** After gating the II-afferent BASELINE with
  the stance gate (runner.py; the ungated i0_ii=1.0 nA was a global ~0.2
  co-contraction floor on every extensor MN - found via diag_phase.py),
  joint motion appears: knee range 0..+17 deg, ankle -87..0 deg. NOTE ON
  SIGNS: with the true flexion-NEGATIVE knee convention, that +17 deg is
  EXTENSION/hyperextension-side motion, not flexion - the swing knee
  still does not flex in the walk; the ankle sweep is dominated by
  plantarflexion (toe-pointing). Same observation, correct labels.
- NEXT (in order): (1) swing-knee FLEXION (negative excursion) - the
  flexors have full range now; the co-contraction standoff is the
  remaining foe (phase table: quads 0.25 vs hams 0.27 in swing); consider
  stance-gating the quad drive harder + the PFA clamp 1.5->2.5;
  (2) ankle plantarflexion overshoot (-61 deg after the rig fix, was
  -87 = toe-pointing; raise the gated II floor / tib_ant swing drive or
  lower E2 ankle_pf);
  (3) COM bounce/drift - BAL gains (kx=150 saturates max_current=6 at
  com_x -0.3) and stance duty 0.29 vs human 0.6; (4) wean rig springs
  (tx tether 1500 N/m, pelvis rotations 150 N·m/rad, NEW 2026-09-10 PM
  after Ben watched the viewer: lumbar_extension/bending/rotation
  150 N·m/rad + hip_rotation 50 N·m/rad - the torso was flipping upside
  down over the free lumbar hinge and the legs free-spinning on their
  long axes; bounce min-COM improved 0.62->0.71) toward true balance;
  (5) then vestibular/ocular (Ben) and cerebellum/BG layers. Tools:
  diag_stab.py (plant A/B), diag_phase.py (phase-aligned activations),
  _muscle_direction_test.py (per-muscle real-activation direction test).
  runner.py --view now plays in REAL TIME with an initial side-view
  camera (orbit: left-drag, zoom: scroll, pan: right-drag, double-click
  a body to track it).

## 2026-09-10/11 night session (Ben staged tuning plan)
1. Literature (Ivanenko 2002 JN air-stepping): preferred cadence at 100 percent BWS = 36+-8 steps/min ~ 0.3 Hz cycles - air stepping is ~3x SLOWER than walking, not 5x faster. Extensor duty 53+-3 percent, joints near-sinusoidal.
2. Deafferented air, no interleg (--no-afferents --no-ground --no-interleg): WORKING. 21 s clean, 0.58 Hz, knee -75..+15 deg true swing flexion, hip -10..+60. Keys: decouple amplitude from frequency (low DRIVE shrank all amplitudes to sub-mV - doubled rg_to_pf/pf_to_mn/descend conductances, slow ADAP tau 1.9 s); subtalar/mtp ligament surrogate springs 10 N-m/rad (OpenSim stops lost in conversion - foot flopped to 128 deg); follower armature 1.0. Remaining: ankle PF-heavy, E-duty 0.16 vs 0.53.
3. Afferented ground + interleg (--drive 2.5): 21 s, stayed up, 6 bursts, COM y +-2-3 cm, adduction +-8 deg (BAL_LAT abductors engaging under load). PELVIS LIMBO: pelvis_tilt still 36 deg with 400 N-m/rad rig + IMU trunk controller - the IMU levels the TORSO (lumbar counter-tilts) but pelvis pitch is driven by SATURATED hip extensors (glut_max 0.94, semimem 1.0) pitching the planted-leg pelvis backward. Trunk muscles cannot fight hip torques - needs the IK back-solve.
4. New pieces: interleg toggle (--no-interleg); BAL_LAT_R/L stance-gated abductor strategy (mujoco -y = opensim +z); BAL_TRK_EXT/FLX IMU vestibular surrogate (torso up-vector PD -> ercspn/obliques). W_PF_MN: E1 knee_ext 0.10, E2 ankle_pf 0.35, F1 knee_flex 1.80, F2 ankle_df 0.45.
5. NEXT - IK/SO back-solve (Ben plan): OpenSim IK walking kinematics -> map to MuJoCo joints -> per-timestep NNLS back-solve along trajectory -> fit W_PF_MN, extract human hip/pelvis torque balance (pelvis_tilt +-5 deg), EMG ordering. Then ground duty 0.6, ankle balance, wean rig springs.

## 2026-09-11 optimization (Ben asked for fast + efficient tuning on ground)
SNS-Toolbox ships no optimizer; literature standard for expensive gait sims = Bayesian optimization (Calandra 2014 biped BO, Antonova 2016, Ryu/Geyer 2021 CPG optimality). Built optuna_walk.py: Optuna TPE, sqlite persistence (optuna_walk.db, resumable), runner --eval mode (12 s schedule, quiet, metrics dict: nan/dx/kz/tilt_max/knee_min/hip_amp/burst_r).
v1 (60 trials, 8 params incl. walk_drive/rg_adapt/desc_e/rg_to_pf/e2_pf/f1_df/f1_kf/post_kneext/kx): best 1.071 - STABLE upright shuffle: 26 s up at tether 1500/300/80 N/m (identical), 0.82 Hz, COM y +-5 mm, adduction +-9 deg, but knees 0..+13 ext-side (stiff-legged). Lesson: stability penalties dominate motion rewards - shuffle wins.
v2 (40 trials, added post_hipext param + motion-demanding objective: knee reward saturating at 35 deg): the v1 winner ENQUEUED as seed stayed best (1.091) - deep-flexion attempts all scored worse on ground. Honest result: at current fidelity, ground stability costs swing flexion; the missing piece is balanced activations (IK back-solve), not more parameter search. v2 objective shaping is ready for a rerun after that.
runner --best loads best_walk_params.json (drive included). runner npz now saves qfull (full qpos) + cfg air/ground; render_frames/render_video are npz-driven (render exactly the recorded run, config auto-detected, mj_fwdPosition so muscles follow bones). Ground media: ground_walk.gif, ground_walk_phases.png, ground_fig*.png in Dissertation folder. SELF-SUPPORTED: not yet - vertical rig still carries weight; BAL cannot hold sagittal/lateral alone; next lever = IK back-solve then wean rig.

## Notes for next session (Ben, 2026-09-11 evening)
1. fig4_gait_cycles.png OpenSim overlay: Ben suspects DOUBLE degree-conversion (our load_benchmark does np.degrees() on the .mot values - verify the freshly generated subject01_walk1_ik.mot units first: hip_flexion_r ~ +-0.3 means radians (keep conversion), ~ +-20 means already degrees (drop np.degrees)).
2. Rename the sim legend entry in fig4 to SNS sim mean (distinct from OpenSim IK).
3. ground_fig4_gait_cycles.png (dissertation folder) should ALSO carry the OpenSim overlay (plot_run was run before the benchmark file existed for that variant).
4. OpenSim 4.3 CLI works locally (Scale+IK in 90 s, artifacts committed in Gait2392_Robotbody); stay on 4.3 unless Python-API scripting is needed (then 4.6 into a py3.11 env, side-by-side, nothing on PATH).

RESOLVED same evening: (1) subject01_walk1_ik.mot says inDegrees=yes - values
are degrees, the np.degrees() in load_benchmark was the double conversion
(plot_run.py now honors the header flag); (2) legend renamed "SNS sim mean";
(3) rerun plot_run.py on a ground run and refresh the ground_fig* copies in
Dissertation\CPG_airstepping_figs (NOT yet done - do it on the post-v3 run).

## 2026-09-11 night session (EB475WS4): the IK/NNLS back-solve chain

Ben's staged chain EXECUTED: MuJoCo-vs-OpenSim muscle validation ->
per-timestep NNLS activation back-solve along the IK trajectory ->
W_PF_MN/W_POSTURE refit (limbo fix) -> v3 optimizer rerun with the
balanced pattern -> rig weaning.

New tools (spinal/): `bsolve_ik.py` (validation + back-solve, writes
bsolve_out.npz/.png/report), `fit_pf.py` (refit -> fitted_walk_params.json),
`wean_rig.py` (--rig-scale ladder), `optuna_walk.py` v3 (study
ground_walk_v3), probes `diag_force/diag_frontal/diag_knee/diag_trans.py`
(keep - they are the regression tests for the traps below). runner.py
gained --fitted (load fitted_walk_params.json), --rig-scale S (all RIG
stiffness*S, damping*sqrt(S); ligament surrogates untouched), --best
pf_gain support, and an eval `duty` metric.

VALIDATION (subject01_walk1_ik.mot: 121 frames 0.5-2.5 s, cycle 1.23 s,
duty 0.70; OpenSim reference via opensim-cmd 4.3 AnalyzeTool on
subject01_simbody.osim - NOTE 4.x has NO standalone StaticOptimization
tool, SO runs as an ANALYSIS inside Analyze, and the lengths file is
`*_MuscleAnalysis_Length.sto`):
- Muscle lengths: median r 0.89, all 78 muscles r>=0.81, 45/78 r>0.9.
  Large RMSE on thigh biarticulars (~5.8 cm at r=0.96) is the
  subject01-vs-generic-gait2392 SCALE offset, not shape error - r is the
  metric. Moment arms r 0.74-0.95 with 100% sign agreement after folding
  in MuJoCo's transmission sign (qfrc = -F dl/dtheta). Conversion sound.
- NNLS activations vs OpenSim SO: glut max/med 0.7-0.93, tib_ant 0.70,
  lat_gas 0.73; peronei/tib_post negative (subtalar geometry + scale);
  median r 0.22, mean torque residual 0.33. Group profiles:
  bsolve_groups.png.

THREE SILENT TRAPS (each verified by a probe script; do not relearn):
1. PELVIS SLIDES LOAD IDENTITY. The converter preserved OpenSim's
   coordinate VALUES (keyframe qpos[pelvis_ty]=0.95 <-> pelvis world z
   0.95; qpos[pelvis_tz]=+0.1 moves the pelvis to world y=-0.1). The z-up
   remap lives in the BODY FRAMES/AXES. A y/z swap in the loading code
   puts the pelvis 2 cm above ground and 1 m lateral (diag_trans.py).
   GRF VECTORS still remap os(x,y,z) -> mj(x,-z,y).
2. EQUALITY COUPLERS POISON INVERSE DYNAMICS. At deep knee flexion the
   polyfit pathpoint followers (fitted near the straight keyframe)
   generate ~890 N*m of spurious constraint wrench at the knee row
   (diag_knee.py). The bsolve ID model zeroes follower armature AND
   disables equalities (mjDSBL_EQUALITY): clean tree ID along the
   measured trajectory (follower dofs are excluded rows anyway; their
   boundmass inertia is virtual). FORWARD SIM STILL NEEDS the
   equalities - ID only.
3. ACTUATOR FORCE CHANNEL: the actuators are dyntype=muscle -
   data.act IS the activation and drives actuator_force under
   mj_forward; data.ctrl is the excitation target and is INERT under
   mj_forward (only feeds act dynamics in mj_step). ctrl=1 force 0,
   act=1 force -2655 N on soleus_r (diag_force.py).

Back-solve formulation: per IK frame, tau = mj_inverse(q, v, a) with
measured GRF subtracted via mj_jac at the CoP on calcn_r/l (forces
remapped, see trap 1); solve lsq_linear with the SO-style ridge
(lam = 0.05 * median column norm) and bounds [0,1]. Plain NNLS is
FORBIDDEN here - it dumps a~1000 into muscles whose FL~0 (near-zero
columns soak residual). 6 Hz zero-phase filtering of the kinematics
matches SO's lowpass_cutoff_frequency_for_coordinates=6. Frontal
signs (adduction/subtalar/list/rotation) are UNSEEABLE by muscle-length
matching - bsolve_ik sweeps them against ID+GRF consistency (frontal
residual 627 -> 30 N*m; found hip_adduction_l -1, subtalar_l -1,
hip_rotation_r -1, hip_rotation_l +1; signs stored in bsolve_out.npz).

REFIT (fit_pf.py): per functional group, NNLS of the back-solved
activation profile over the four recorded PF-cell phase windows (from
spinal_run.npz; only right-side PF cells are logged - windows are
identical per side). THE LIMBO FIX, at the source: hip_ext stance drive
stack 1.17 -> 0.27 vs human peak 0.34 (the old W_POSTURE hip_ext 0.22
alone was ~2/3 of the human PEAK, and E1+E2 stacked on top -> glut_max
0.94/semimem 1.0 saturation -> pelvis limbo). Output:
fitted_walk_params.json (full W_PF_MN/W_POSTURE + pf_gain=1.0).

V3 OPTIMIZER (optuna_walk.py, study ground_walk_v3): new pf_gain
dimension (log 0.5-8) - the back-solved weights carry honest human
amplitudes ~10x SMALLER than the hand-tuned table the network's gain
structure was tuned against (old F1 knee_flex 1.80 vs fitted 0.148); at
gain 1 the sim barely moves (hip_amp 0.85 deg). Plus the 5 phase
multipliers [0.5,1.8], drive/rg/kx ranges, and a duty-0.60 reward.
CAUTION: do NOT pass --fitted inside the optuna eval call - the reload
wipes the trial's mutations (verified: identical metrics with and
without gain). The final config reproduces with `runner --fitted --best`
(--best applies pf_gain to the whole table; the 5 knob values in the
json are already effective and are written AFTER the gain to avoid
double-apply).

WEANING (wean_rig.py): runs the ladder rig-scale 1.0 -> 0.6 -> 0.4 ->
0.25 -> 0.15 with --fitted --best; a stage passes when the full schedule
completes, COM height > 0.62, pelvis tilt < 35 deg; stops at first fail,
keeps wean_stage*_S*.npz per stage.

## 2026-09-11 night RESULTS (v3 winner + weaning verdict)

- v3 study ground_walk_v3, 40 trials: best score 1.521 (trial 37) vs
  v2's 1.091. Winner: pf_gain 0.52 (the optimizer went DOWN from the
  seed - human-shaped patterns need LESS brute-force gain, the opposite
  of the hand-tuned table's direction), drive 1.49, rg_adapt 0.83,
  desc_e 1.33, rg_to_pf 1.98, kx 220. Saved in best_walk_params.json
  (effective values + pf_gain + multipliers; reproduce with
  `runner --fitted --best`).
- v3 winner 12 s eval: no NaN, kz 0.877, tilt 27.1 deg, dx 0.151 m,
  9 bursts - but knee_min -0.6/hip_amp 0.9 deg on the SHORT window.
- FULL 22 s ground run at rig-scale 1.0 (wean_stage0_S1.npz): PASS -
  knee -22.4 deg REAL swing flexion, hip amp 37.7 deg, tilt 29 deg,
  kz 0.87, dx 0.15 m. The v2-era "knees 0..+13 ext-side stiff shuffle"
  is GONE at full rig.
- Weaning ladder: S=0.8 completes 22 s and keeps stepping (knee -22.7,
  hip 37.4) but tilt 36.7 deg = progressive lean; S=0.6 tilt 46.3 deg.
  VERDICT: support boundary S~0.8-1.0. The pelvis rotation assist is
  what the weaned rig misses - the BAL_TRK IMU levels the TORSO, the
  (now human-scaled) hip extensors hold the legs, but nothing yet holds
  the PELVIS pitch in the frontal-sagittal sense. Next lever per Ben's
  plan: ground duty 0.6 + ankle balance, then the pelvis-balance piece,
  THEN wean below 0.8.
- Figures: plot_run.py rerun on wean_stage0_S1.npz (spinal_run.npz
  currently holds that run) with the inDegrees fix - fig4 now carries
  the OpenSim IK overlay with the "SNS sim mean" legend; ground_fig1/3/4
  + ground_walk.gif refreshed in Dissertation\CPG_airstepping_figs.

## 2026-09-11 late night: publication circuit figures (draw_circuit.py rewrite)

draw_circuit.py is now a figure SUITE (Ben asked for paper figures of the
circuit): `python draw_circuit.py [--which core|full|weights|all]
[--source params|fitted|best] [--fmt pdf,svg,png]` -> figures/.
  - circuit_core: one side, DRIVE/POSTURE/POST_i -> RG+ADAP -> PF(+PFA)
    -> MN ellipses; reflex block; each MN pool carries ONE label (its
    dominant W_PF_MN entry); edges drawn for w >= 0.03; full table =
    the weights figure.
  - circuit_full: both sides (left ghost) + interleg commissurals as
    arcs over the top + BAL family; rig marked external.
  - circuit_weights: W_PF_MN x groups heatmap + W_POSTURE column.
CRITICAL PROPERTY: numbers are never hardcoded - --source composites
params.py -> fitted_walk_params.json -> (best) pf_gain x knobs exactly
like `runner --fitted --best` (verified by spot-check vs the jsons), and
each figure prints a gray source note. RERUN THIS AFTER ANY RETUNE
(e.g. the v4 trunk-Fmax fix) so figures never drift from the sim.
Gotcha recorded the hard way: matplotlib arc3 with a +x chord bulges
DOWN for positive rad - over-the-top arcs need negative rad. Old
spinal_circuit.png is superseded (its weights were the stale pre-refit
hand-tuned table). Copies of all 9 files (pdf/svg/png x 3) live in
Dissertation\CPG_airstepping_figs\ as circuit_*.

## 2026-09-11 night, addendum: Ben's knee/patella + muscle-count notes

- MOMENT-ARM DEFINITION (Ben): OpenSim's arm = dl/dtheta with the whole
  geometry following. MuJoCo's data.actuator_moment is the RAW Jacobian
  and does NOT propagate through the eq couplers - at the knee it misses
  the 36 vastii pathpoint followers, at hip_flexion 5 more. bsolve_ik.py
  now computes moment arms by CENTRAL DIFFERENCE dl/dtheta with the
  followers re-projected at every perturbed pose (fd_moments), which is
  the matching quantity AND the right B matrix for the back-solve.
  Keep this in mind anywhere knee moments matter.
- VASTII -> TIBIA VERIFIED (Ben: "the opensim model has no patella"):
  all 36 knee eq couplers drive the vas_med/vas_int/vas_lat pathpoint
  bodies (P3-P5) on both sides - the vastii actuate the knee/tibia
  through them (rect_fem rides vas_med's P4 since the patella reroute).
  Quantitatively (fd_moments vs OpenSim MuscleAnalysis knee arms):
  vas_med/lat/int mean arm -0.047 m vs OpenSim -0.045 m, same sign, but
  the small AC fluctuations anti-correlate (r -0.67..-0.81) -> the knee
  coupler polys deviate from OpenSim's true pathpoint trajectories away
  from the keyframe (polyfit artifact, same family as the range-limit
  and ID-wrench issues). Same likely explains the residual peronei/
  tib_post SO mismatches (subtalar couplers) and knee median arm r 0.31.
- "ALL 78 MUSCLES" RESOLVED: the comparison table = 92 actuators MINUS
  14 with Fmax (gainprm[2]) <= 5 N. Only 6 are the intentional prunes
  (quad_fem/gem/peri r/l). THE OTHER 8 SHIP AT Fmax = 1 N FROM THE
  CONVERTER: ercspn r/l, intobl r/l, extobl r/l, ext_hal r/l. That
  means the IMU trunk controller (BAL_TRK_FLX/EXT -> ercspn/obliques)
  has been driving muscles with 1 newton of capacity - the torso is
  held by the rig springs alone. diag_fmax.py prints the audit.
  RECOMMENDED FIX (needs Ben's go - it invalidates the current v3
  tuning): set the 8 Fmax values from stock gait2392_thelen2003 in
  patch_xml (like the prune patch), re-run the standing solve, refit,
  v4. Until then, trunk-control claims about ercspn/obliques are
  placeholders.
