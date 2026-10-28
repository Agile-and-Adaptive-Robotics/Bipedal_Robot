<!-- SESSION_1_SIMULATION_RECOVERY_20261027_START -->
## Session 1 simulation recovery update

**Current task-5 disposition: failed crossed-leg diagnostic.** A subsequent
pelvis-frame audit confirms reversed ankle order in 91.8% of post-ramp samples.
The reduced-support numerical replay is reproducible but is not acceptable
physical standing/stepping. Ideal reference tensions are still used; the
optimized BPA hardware has not been imported. Details below supersede any
earlier candidate terminology.

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

### Supported-stepping continuation

Replayed the recorded robot simulation for Ben and added a stepping assessment
that separates neural rhythm, sustained contact loss, actual mesh clearance,
and progression. The drive-3 run has nine RG bursts per side and a left/right
correlation of -0.720, but maximum foot clearances are only 0.19/0.87 mm.
It fails bilateral swing and lifting gates and is not verified stepping.

Three gradual support-cap trials failed the pelvis-height gate: 90% actual
vertical-force cap at 17.890 s, 85% at 15.708 s, and 80% at 14.882 s.
Increasing posture gains alone also failed; these are retained as failed
diagnostics, not successful weaning. A torque-capacity check found substantial
available knee torque, while the existing joint gravity compensation omitted
explicit stance-load feedforward. The next bounded tests add the unsupported
body load through muscle torque allocation and retain unchanged external
rig-force accounting. Controller settings and outcomes are saved per run.

### Stance-load and activation-floor results

Stance-load compensation completed 25 s with actual vertical support capped
at 85%. Assessment now starts after the ramp ends (15 s), so pre-wean foot
lifting is excluded from low-support results. The direct-RG pattern still
failed repeated bilateral foot-lifting gates in that interval.

Ben asked about a 5% coactivation floor. The previous solver allowed zero.
Added an adjustable minimum allocated activation, delivered through SNS motor
outputs and the existing activation dynamics. This floor is not a measured
coactivation index. Compared zero and 5% without assuming a benefit. Both
larger-motion trials completed 25 s at the 85% cap; the 5% case had peak
post-ramp clearances of 9.6/18.2 mm, but failed repeated bilateral lifting.

An RG-event-anchored pattern separates early knee lift from hip advancement.
At the 85% cap, its zero-floor 25 s trial has five right and nine left
clearance-qualified swing bouts after the ramp, with peak clearances of
12.7/101.5 mm. Actual post-ramp support averages 84.99% body weight. Neural
rhythm and bilateral lifting gates pass, but backward drift of 0.153 m
remains; this is a supported stepping candidate, not forward walking.
At a 75% cap, the trial survived but repeated right-foot lifting failed.
The corresponding 5% floor at an 85% cap also failed right-foot lifting.
The default floor therefore remains zero while tuning continues; the 5%
setting remains available and is recorded in each relevant run.

### Further support reduction

Explicit swing-leg unloading plus a larger right-knee reference completed a
30 s trial with a 5% activation floor and an actual 75% vertical support cap
after 15 s. Its matched zero-floor case failed lumbar posture at 19.408 s.
The 5% floor therefore helped survival in this configuration, although its
right-foot clearance remained below the 5 mm repeated-lifting criterion.
The recorded 5%-floor activation minimum is 0.05 within floating-point error.

Added a foot-apex phase check: neural antiphase alone is insufficient when
the physical feet lift together or at inconsistent phases. The earlier
zero-floor 85% candidate fails that physical phase check and is not promoted
to verified stepping. Higher tracking gains, activation-feedback gain 3,
and initial inertia-aware posture trials failed and are retained as diagnostics.

The next bounded trials explicitly add RG-timed lateral pelvis assistance
and stronger horizontal centering while retaining hard 75%/65% caps on total
vertical rig force. These extra horizontal support forces are logged and
must be disclosed in any resulting supported-stepping claim.

### Verified 65% support checkpoint

The lateral-assistance configuration completed 40 s in fresh reference and
exact-accelerated processes, with bit-identical full trajectories and matching
source hashes. The final 25 s are capped at **65% actual vertical support**;
the activation floor is 5%. This reduces vertical assistance substantially
from the initial approximately 95% posture test, while horizontal and angular
support remain active. Mean absolute horizontal support forces are 13.20 N
anterior/posterior and 31.53 N lateral; mean absolute support torques are
13.93/3.50/17.24 N·m about x/y/z. The lateral support target moves 30 mm with
neural phase. Total model mass remains 75 kg.

Post-ramp maximum foot clearances are 55.2 mm right and 21.8 mm left, with six
clearance-qualified swing bouts per side. Net post-ramp forward displacement
is only +8.2 mm. Neural rhythm produces 17 bursts per side, but physical foot
cycles skip many neural cycles. A stricter phase-coverage gate therefore
classifies this as reduced-support motion with intermittent bilateral lifting,
not regular supported stepping or forward walking. Maximum lumbar rotation
in this interval reaches 17.1°, another reason to keep tuning.

Evidence: `robot_torso\supported65_backend_verification.json`,
`runs\verified65_reference.*`, `runs\verified65_exact_fast.*`, and
`runs\verified65_exact_fast_assessment.json`. The runner now snapshots source
at run start; the verified runs retain their exact runner copy and source
hashes. Construction and dynamic-route checks still pass with zero reversed
neutral actions. Slower intrinsic RG time constants are being tested to address
the mismatch between neural cadence and completed physical swings.

### Cadence and interpretation checkpoint

The intrinsic RG time-constant comparison completed 40 s at tau_h = 0.7 s,
with the same 65% post-ramp support and 5% activation floor. It produced ten
neural bursts per side in the assessment interval, but physical foot-phase
coverage still failed. Tau_h = 0.525 and 1.05 s trials lost lumbar posture.
Slowing the rhythm alone therefore does not resolve the gait mismatch.

Ben observed apparent trip recovery, backing up, and duck-footing in the replay.
No trip perturbation was applied, so no trip-recovery claim is made. Hip yaw
tracking remains poor despite a straight-hip reference. No five-stage curriculum
was run on this CAD model. This work reused the existing layered SNS architecture
and starting parameters, then added bounded posture/stance-load allocation,
compensated MN target currents, and RG-event-anchored motion references.
The controller is a hybrid bootstrap with substantial rig assistance.

The next curriculum must include actual hip yaw and foot heading, repeated
clearance-qualified alternating swings, mechanical knee-range checks, trunk
posture, and measured support force. The current full-trial evidence is a
repeatable reduced-support checkpoint, not a completed gait curriculum.
The current replay is `robot_torso\robot_support65.gif`; the machine-guarded
reproduction preset is `robot_torso\run_support65.ps1`.

### Exact controller architecture clarification

The replay used the base `spinal\build_network.py` builder, the s3k-lineage
phase-PF architecture with `joint_pf=0`, `syn6=0`, `full_rules=0`, `ia_in=0`,
and `renshaw=0`. It did not use w2lvar, syn6, separate joint PF half-centers,
or dedicated biarticular PF layers. Each side has RG-E/RG-F persistent-Na
half-centers and four shared PF phase cells E1/E2/F1/F2. Each of the 54 muscles
has an MN and Ia/II/Ib encoder. Functional antagonists have reciprocal Ia
inhibition, using direct afferent-to-antagonist-MN edges in this configuration.
Biarticular muscle PF weights combine the primary functional group with half
the secondary group's weight; the legacy map treats sartorius and gracilis
as single-group muscles and does not fully encode anatomical biarticularity.

Crucially, `motor_target_current` cancels the instantaneous net synaptic
current entering each MN to track the allocator's activation target. PF and
local reflex connections therefore exist in the network but do not have
their normal effective control authority over muscle activation in this
bootstrap. Effective motion comes from RG timing, engineered motion targets,
muscle torque allocation, and the support rig. Earlier broad descriptions of
PF/reflex involvement should be read with this explicit limitation. A proper
PF-driven neural gait and its five-stage curriculum remain to be implemented
and verified on the CAD robot.

### Routing-commit review and functional-design prompt

Reviewed the local and GitHub branch history since 26 October. Commit
34a6431d adds the survey3 originals, capability-v3 comparison, and VasMed
minimum-BPA results. The survey2 driver explicitly fixes master's route points
and changes BPA count, tackle, and resting/tendon lengths; these results are
not automatically updated CAD attachment routes. The current replay does not
import those optimized hardware parameters or use `bpa_muscle.py` forces.

Measured foot headings and ankle ordering in `foot_heading_audit.json`:
median headings +28.25° right/-32.04° left, minimum right-minus-left ankle
lateral spacing -0.328 m, and reversed order in 91.8% of 20 Hz samples after
15 s. This confirms the user's crossed-leg concern. The assessment now includes
an uncrossed-ankle-order gate. Mesh self-collision and clearance must also be
addressed; ankle ordering alone is not a sufficient collision test.

Prepared a full easteregg2-specific ZCode prompt at
`Code\Matlab\Mesh_Optimization\Functional_Routing_20261027\ZCODE_FUNCTIONAL_ROUTING_PROMPT.md`.
It proposes full signed torque-vector and attainable-set optimization, including
combined hip functions, antagonist coverage, physical BPA constraints, geometric
route variables, and staged validation. It separates functional actuator design
from two-/four-cell PF-module choices. No new design optimization was launched,
no message was sent to ZCode, and no existing route geometry was changed.

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
