<!-- SESSION_1_SIMULATION_RECOVERY_20261027_START -->
## Session 1 simulation recovery update

Reporting window: 26–27 October 2026. Updated 27 October 2026 on **easteregg2**.
Repository root: **`D:\GitHub\Bipedal_Robot`**. Python:
**`D:\Anaconda\envs\myo\python.exe`**. Paths below are relative to that repository
root unless given in full. This block records Session 1's simulation work;
Session 2's entries and the existing advisor/dissertation record are preserved.

### Verified W2L recovery and speed changes

- Corrected gain registration at the evaluator boundary, where the previous
  fresh-process verification discarded per-target settings. Unknown settings
  now raise an error rather than disappearing silently.
- Recovered trial 458 with all effective settings delivered. Its 30 s trajectory
  has 28 left and 21 right contact events, median phase 0.496063, forward
  displacement 2.246231 m, and fitted speed 0.091786 m/s. Four independent
  replays match exactly. All 12 fresh-process perturbations of the saved numeric
  settings by ±1%, excluding the fixed horizon, pass the supported-stepping gate.
- Audited 56 effective conductance classes and three sensor dimensions. Two
  saved settings have no corresponding edges. These are conductance classes,
  not 56 independently tuned anatomical synapses.
- The support ladder reaches nominal `wsup=0`, but the PD harness still carries
  about 16.76% of body weight. No harness-free walking claim is supported.
- Replaced the unsafe sparse-equivalence claim with an exact backend that
  preserves reduction order. Full 30 s traces match the reference exactly;
  isolated neural-step speed improves 1.425×. Independent processes use one
  BLAS thread each. The old FastTrial physics loop remains unsuitable.
- Re-ran perturbation panels with effective gains. Commissural deletion fails
  bilateral stepping; removing all feedback fails the contact-phase gate while
  cycling persists. These are supported W2L adaptation experiments, not an
  original Deng replication.

Evidence: `Code\MuJoCo_SNS\spinal\reports_20261027_codex\REPORT.md`,
`verification.json`, `full_basin_results.json`, `speed_parity_exact_fast.json`,
`perturbation_results.json`, and `w2l_verified_config.json` in that folder.

### Di Russo replication

Found the original author repository, `DiRussoAndrea/Spinal_controller`, and
custom `simgait/SCONE` fork. The previously cited neumovelab repository is a
comparative reimplementation. The original model is the reduced planar
`gait9dof18musc_Thelen_20170320.osim`, with nine active coordinates and 18 Thelen
muscles. It initializes locally in OpenSim 4.4.1. The optimized parameter file
contains 256 parameters, now inventoried.

Corrected multiple equation, sensor, and initialization errors in the existing
Di Russo-inspired adaptation. Its compiled neural implementation matches its
corrected Python equations, but that does not establish equivalence to the
author's controller. A 960-trial duration/speed search selected a somersault;
the candidate is rejected as walking. Further optimization of that surrogate
was stopped. Published walking remains unreplicated, and the original custom
SegmentalController runtime is not available locally.

Original sources are staged at `tmp\codex_dirusso_original\` and
`tmp\codex_simgait_scone\`. Evidence and rejected-result plots are in
`Code\MuJoCo_SNS\spinal\reports_20261027_codex\`.

### Actual robot geometry and corrected conventions

- Refreshed the open `10_AH_001.SLDASM` after Ben adjusted the left knee.
  Preserved both four-bar knees and both spherical hips. Current joint anchors
  meet within numerical precision. Source CAD was not saved by Session 1.
- Corrected Session 1's earlier frame mistake. The current robot pipeline uses
  **x anterior, y superior, and z right**, including gravity, ground plane,
  geometry, inertia, and muscle geometry. A pelvis-part datum transform removes
  the assembly's global repositioning while retaining the knee adjustment.
- Ben superseded the lumbar hinge with a **ball-and-socket joint**. The spine
  remains a rigid 0.25 m bar. Payload center is 100 mm anterior and 75 mm above
  the joint. The current model splits the balancing payload equally between
  torso and a rigid pelvis/ilium load.
- Loaded `Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design,
  Newest (7-17-20)\Foot Ankle.SLDASM`: all 22 components resolved with no open
  errors or warnings. Exported the geometry, mirrored it for the opposite side,
  and welded all foot internal joints and the ankle attachment to each shank.
- Total mass remains 75 kg: 5.356664 kg CAD lower body with feet, 0.12 kg spine,
  and 69.523336 kg balancing payload. Each payload portion is 34.761668 kg.
  Neutral torso extension demand is 34.1012 N·m. Printed parts use solid Onyx;
  hardware densities and the 20 mm square spine cross-section remain engineering
  assumptions rather than measured hardware properties.

The current model is `Code\MuJoCo_SNS\robot_torso\ah001_neural.xml`.
Earlier z-up routing files, including `opensim_routes.json` and
`ah001_with_opensim_routes.xml`, are superseded and must not be used as the
current robot model. Legacy converted human/W2L models were not globally
rotated; their verified baselines remain separate.

### Robotbody routes and neural control

The source is now **`Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim`**,
using its right-side routes and explicit left mirrors. The existing
`Code\Matlab\Mesh_Optimization\Routing27_v2\actuator_map_27.json` selects
**27 actuators per side, 54 total**. No new hip-muscle consolidation was performed.
The full source inventory retains other route definitions for auditing; they
are not additional robot actuators.

Moving path points use exact piecewise cubic functions exported from OpenSim,
including derivative and endpoint checks. Conditional membership updates at
runtime. Complete path Jacobians include moving-guide terms in force and
velocity calculations. Checks cover 492 source function evaluations, 25 native
conditional checks, and 75 constrained poses. Function errors are below
1.4 × 10⁻¹⁷ m, and full path-derivative errors are below 1.4 × 10⁻⁹ m.

The updated audit reports **zero reversed neutral sagittal actions**. Mirrored
attachment corrections for adductor magnus 2, gracilis, and semimembranosus are
recorded explicitly; maximum point shifts are 15.6, 21.0, and 18.1 mm. The
previous TFL reversal clears with the refreshed CAD pose. Conditional switching
of unselected semitendinosus routes has small length jumps after retargeting;
those routes are not actuators in the 27-per-side model and remain flagged.

Constructed a **262-neuron, 230-input SNS controller** driving the 54 actuator
controls, with Ia/II/Ib and contact inputs, postural force allocation, and
RG-driven stepping references. Reference tensions have a 30 ms activation
filter and the existing 27-actuator map's MIF caps, including 2500 N for each
erector. These are simulated reference capacities, not a calibrated BPA force
model. The latest test mode uses explicitly compensated MN target currents;
it is a hybrid controller, not a faithful Deng circuit replication.

**Current control status: a 20-second supported crouched-stance test is verified;
locomotion remains unverified.** Initial tests collapsed or lost posture. After
correcting force allocation, motor-neuron tracking, and the explicit support rig,
the final five seconds of the successful stance run had knee excursions of
1.35° right and 0.95° left, maximum pelvis tilt of 3.92°, and maximum lumbar
rotation of 0.58°. Mean foot loads were 17.63 N right and 21.61 N left. Actual
vertical rig support averaged **94.68% of body weight** in that interval, with
angular stabilization also active. This establishes supported postural control,
not independent balance or harness-free standing. Disabling all muscle forces
under the same support settings caused lumbar failure after 0.202 seconds.

The reference and exact accelerated neural backends produced **bit-identical
complete 20-second stance trajectories**, including positions, velocities,
activation, force, contact, route, and neural logs. Evidence is in
`Code\MuJoCo_SNS\robot_torso\controller_assessment.json` and
`runs\stand_20s_reference.*` / `runs\stand_20s_exact.*`. A replay image and
animation are `robot_torso\neural_stance.png` and `neural_stance.gif`.

The original stepping drive latched rather than oscillated. A separate rhythm
probe found bilateral oscillation at drive 3.0, and the corresponding supported
robot run completed 20 seconds without a numerical or posture failure. Its
whole-run forward displacement was 0.446 m, and mean actual vertical support
was 93.74% of body weight. These totals include settling and do **not** establish
alternating gait or locomotion; contact, phase, and displacement after stepping
onset still require assessment. Reducing vertical rig stiffness from 800 to
300 N/m failed the pelvis-height gate after 4.352 seconds, before stepping
onset. These runs are `runs\step_drive3.*` and
`runs\step_drive3_lower_rig.*`. Task 3's `s3k`, `w2lvar`, and `syn6` campaign
has not been launched.

Current evidence: `Code\MuJoCo_SNS\robot_torso\robotbody_routes.json`,
`robotbody_source_inventory.json`, `route_sign_corrections.json`,
`dynamic_route_verification.json`, `run_neural_robot.py`, and `runs\`.
Fresh CAD/foot inventories are in
`Code\MuJoCo_SNS\spinal\reports_20261027_robot_control\`.

### Scope and next steps

Continue stepping diagnostics and inspect contact, phase, and support forces
before classifying locomotion. Supported stance and complete-trajectory backend
parity are now verified within the assistance conditions stated above. Updated
ZCode prompts name easteregg2 and provide full paths, environment, and runnable
commands in
`Code\MuJoCo_SNS\spinal\reports_20261027_robot_control\ZCODE_PROMPTS_SESSION1.md`.
No
commits, pushes, outgoing ZCode messages, or dissertation source edits were
made in this simulation-recovery work. This Session 1 block is additive; it does
not replace Session 2's work or resolve advisor comment threads.
<!-- SESSION_1_SIMULATION_RECOVERY_20261027_END -->
