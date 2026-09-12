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
