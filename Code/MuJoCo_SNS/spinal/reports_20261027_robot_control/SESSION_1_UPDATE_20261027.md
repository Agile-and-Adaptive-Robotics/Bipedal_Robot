<!-- SESSION_1_SIMULATION_RECOVERY_20261027_START -->
## Session 1 simulation recovery update

**Current task-5 disposition: uncrossed supported hybrid-control diagnostic.**
The new v2 model completes 40 s with no reversed ankle ordering and 80% actual
vertical assistance after 15 s. Regular stepping, independent balance, and
calibrated BPA actuation are not established. The older support-65 trajectory
remains rejected for crossed-leg mechanics (91.8% reversed ankle ordering).

Reporting window: 26–28 October 2026. Updated 28 October 2026 on **easteregg2**.
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

### 28 October: abductor correction and defense assets in progress

Ben requested a new GIF with more abductor activation, an accompanying still,
and a dissertation-placement assessment. Measured abduction leverage in the
selected 27-per-side layout comes primarily from the gluteal seed route, TFL,
rectus femoris, and sartorius. Increasing every abducting route's floor to 10%
caused collapse; muscle-mediated separation feedback performed better.

Exploratory uncrossed runs were rejected as final evidence when the geometry
audit found up to 16.6 mm knee-loop closure error. Investigation found that
the allocator's full off-manifold constraint-Jacobian nullspace spuriously
removed hip authority. A corrected anatomical tangent keeps six hip rotations,
three lumbar rotations, and two planar knee motions independent. A separate
v2 diagnostic model tightens closure, enables opposite-leg CAD mesh collision,
and documents provisional passive knee stops (-70° to +3° relative to CAD).
A flexed-knee initial stance avoids the straight-knee startup transient.
The original model and older diagnostic XML are preserved.

Current trials, raw states, source snapshots, and activation comparisons remain
in `robot_torso\runs`. Final assets will be generated from accepted recorded
physics only. The latest dissertation already has a CAD integration paragraph
in 60-sns and a detailed subsection in 95-AppendixD (`sec:app_cad_supported_stance`).
The new still belongs with the scoped CAD-tool demonstration in the SNS chapter;
the measured comparison and assumptions belong in the existing Appendix D
subsection. Overall Discussion needs a cross-reference rather than duplicate
figures or a second results table. No dissertation chapter source has been
edited in this continuation.

### 28 October: accepted media, attachment sensitivity, and w2lvar probe

Final `abd_defense_exact` and fresh reference runs reproduce every saved array
bit for bit. In the 15–40 s interval, extra abductor control changes reversed
ankle ordering from 67% to 0%, minimum separation from -38.3 to 142.4 mm,
and right TFL mean activation from 22.3% to 81.8%. Vertical assistance is 80%
in both conditions. The control adds 160 N/m muscle-mediated separation
feedback and a 6% floor on currently abducting routes, against the 5% global
floor. Knee-loop maximum position error is 0.387 mm with control on, and
sampled opposite-leg mesh penetration is zero. Toe-out and very small
post-ramp knee excursions remain; gait acceptance gates fail.

Recorded media are in `Documentation\Reports and Papers\Dissertation\Defense\
slideshow\AH001_20261028`: 96-frame GIF (20–32 s; exact 12 s timing),
2400-by-1400 still (24 s), provenance, alt text, and placement recommendation.
The new defense copy is `Defense\AE_Revision_20261023\output\
Bolen_Defense_Assertion_Evidence_20261028_AH001_v8.pptx`: animation slide 29,
hidden still/comparison slide 39. All 37 original slide identities, text,
and notes are conserved, source v6 is unchanged, and the GIF is embedded
byte for byte. Both added slides were rendered in native PowerPoint and
visually inspected. Use the still in the scoped SNS CAD integration discussion,
with the matched table and assumptions in existing Appendix D subsection
`sec:app_cad_supported_stance`; Overall Discussion needs only a cross-reference.
No master dissertation chapter source was edited.

Ben identified a possible 60 mm tibia/foot interface mismatch. The current
foot weld uses the top of the ankle assembly and distal shank, not a verified
mating datum. Its mount-to-sole height is 103.7 mm. A separate hypothesis
model translates each entire foot 60 mm toward the shank's proximal origin;
it does not shorten tibial CAD, change tibial inertia, or validate an interface.
Total mass remains 75 kg. With the same s3k-lineage hybrid control, median
ankle spacing changes from 159.0 to 57.3 mm, and sampled interleg mesh
penetration reaches 1.316 mm despite no reversed ankle ordering. The
attachment assumption materially affects the result.

`probe_foot_mount_w2lvar.py` adapts the existing w2lvar builder to all 54
actuators in isolated processes, leaving production builders and the baseline
runner unchanged. The actual variant has 544 neurons and 230 inputs; its
548 named indices include four ankle aliases, not extra neurons. Each side
has HIP E/F and shared KNEE/ANKLE E/F PF pairs with inhibitory lamination.
Both original-attachment and 60 mm-shift trials complete 40 s at the same
80% post-ramp vertical support, without reversed ordering or sampled
opposite-leg contact. Median spacing is 63.8/99.0 mm respectively; maximum
knee-loop error is 0.516/0.784 mm. The shifted w2lvar timing trial passes
bilateral unloading and lifting counts but fails repeated alternating foot
apices and forward progression. It is not accepted as regular stepping.
Accelerated/reference neural stepping is bit-identical for 400 seeded-input
steps of this exact 544-neuron circuit; no full-trajectory w2lvar replay or
robustness claim is made.

Crucial limitation: the inherited MN current clamp cancels incoming PF/reflex
MN synaptic currents. This comparison therefore tests network timing and
central feedback within an engineering force allocator, not full w2lvar
muscle-pattern control. A proper PF-driven comparison must retain those
currents and separately bound the postural contribution.

Floating route points are body-fixed virtual attachments or guides, not
automatically snapped to the nearest CAD solid. They transmit force through
the named parent body's Jacobian, assuming ideal rigid offset brackets.
Bracket mass, compliance, packaging, and collision-aware wrapping are not
modeled. This is distinct from the validated source moving/conditional
path-point functions and their complete virtual-work derivatives.

The four-condition table, reproduction commands, and neural-comparison limits
are in `Code\MuJoCo_SNS\robot_torso\FOOT_MOUNT_W2LVAR_RESULTS.md`.
A separate 20–28 s GIF and 24 s still of the hypothetical shifted-foot w2lvar
timing probe are under the final media folder's `w2lvar_foot60_probe\` subfolder.
They are labeled exploratory and are not substituted into the defense baseline.

### 28 October: s3k/w2lvar/syn6 derivation audit and fast trials

Audit and evidence: `Code\MuJoCo_SNS\spinal\reports_20261028_family_audit\
AUDIT_REPORT.md`, `RESULTS.json`, six compiled edge inventories, and the
machine-specific `REPRODUCE.md`. Reviewed current builders, effective gain
loaders, source weight/basis files, corrected wiring rulings, and the primary
Rybak, Deng, Shinohara, Shevtsova, and Szczecinski sources. These families are
literature-informed engineering models, not three exact paper replications.

PF architecture per leg: s3k has E1/E2/F1/F2 phase cells; w2lvar has HIP E/F
and shared KNEE/ANKLE E/F antagonistic pairs; syn6 has six human IK-derived
NMF channels with a mixed S5 outside the default reciprocal inhibition.
The current CAD counts are respectively 262/230, 544/230, and 488/230
neurons/inputs. w2lvar's 548 named indices include four ankle aliases.
Its six trunk actuators have no PF input. syn6 covers all 54 CAD actuator
names but has no per-muscle central sensory projections at default gains.
The suspected w2lvar explicit-zero fallback problem does not occur in the
current weight file. The earlier zero-afferent claim about the original Deng
network is retracted; stale paragraphs in older audit files are identified.

The CAD mechanical check finds 20 of 54 routes without resolved free-joint
torque at five recorded poses, because distal joints are welded. Internal
route tension and sensory activity remain possible; full network connectivity
does not imply 54 useful control directions. Primary human-name groups also
omit some biarticular and hip rotation functions. These are substantive
limitations when transferring the controller to the robot's actual routes.

The human-model screen ran 13 configurations: three frozen baselines and
ten targeted central-feedback, PF-gating, or crossed-pathway changes. The
baseline scores reproduce -161.56754173676563 (s3k), -174.84835562693294
(w2lvar), and -200.37955834081458 (syn6). Every targeted change worsens its
own family's score. The s3k and syn6 baseline diagnostics show a frozen left
leg; the s3k anchor is a numerical regression reference, not accepted bilateral
walking. No winner or existing Optuna database was changed.

Added optional `spinal\w2l_cpg\preallocated_exact_eval.py`: sparse Numba
conductances, reused input/state/reduction buffers, and one-hot I/O maps,
while retaining the reference NumPy synaptic reduction order and fixed-tau
NaP equations. All six family/model combinations pass 400-step voltage and
h-gate parity checks. All three complete human baselines match both the old
exact-fast and reference backend arrays bit for bit. Warmed local neural-step
ratios relative to the previous fast evaluator: w2lvar 1.20x human and 1.41x
CAD; syn6 1.05x human and 1.19x CAD. Small s3k builds are 5–7% slower with
the extra preallocation, so the old exact-fast evaluator remains preferable
there. These are shared-machine neural timings, not whole-trial speed claims.
The new backend rejects non-fixed-tau gated classes rather than silently
changing a paper model's voltage-dependent channel kinetics.

CAD neural-authority experiments now expose the masking by the old MN clamp.
Full immediate release fails lumbar posture at 0.126/0.186/0.066 s for
s3k/w2lvar/syn6. A 5–15 s ramp to 10% uncanceled native synaptic-current
authority fails at 22.242 s for s3k and 13.830 s for syn6. w2lvar completes
the requested 25 s at the same 80% post-ramp vertical support, but still
fails physical stepping gates. A fresh preallocated replay matches every
saved array of that partial-neural run bit for bit.

Delayed PF-to-MN deletion at 15 s preserves the exact pre-intervention
history. In the partially exposed w2lvar controller it changes motor commands
by up to 0.06615 within the first 100 ms (6.615 percentage points). The
fully clamped negative control changes them by only 6.96e-11 in that window.
This verifies actual PF influence in the partial mode, unlike the fully
masked defense baseline. The remaining 90% engineering clamp, posture
allocator, reference tensions, rig, and geometry assumptions prevent any
claim of a pure neural, BPA-validated, or independent walking controller.

Production builders, `params.py`, `runner.py`, dissertation chapters, and the
defense deck were not edited by this audit. Syntax checks, numerical replay,
and scoped whitespace checks pass. No commits or pushes were made.

### Articulated Gait2392 functional PF experiments (2026-10-28)

Ben directed these experiments to the articulated human Gait2392 MuJoCo
model, separate from the CAD robot. Evidence and reproduction instructions:
`Code\MuJoCo_SNS\spinal\reports_20261028_functional_pf\REPORT.md`.
Optional implementation: `spinal\functional_gait2392\`.

The w2lvar extension adds separate four-cell modules (two PF neurons and
two inhibitory relays) for bilateral hip X/Y and foot inversion/eversion,
global lumbar X/Y/Z, separate sagittal ankle drive, and six biarticular
families per leg. Combined census: 980 neurons/400 inputs; temporal-waveform
variant: 980/424. Eighteen proximal/knee/ankle biarticular routes receive
dedicated family PF drive. Long toe muscles do not yet receive dedicated
multiarticular modules. Native SNS motor-neuron outputs drive the human
model without the CAD controller's motor-current cancellation.

An equality-aware tangent projection agrees with finite-difference torque
moment arms to 2.183e-11 m. Recruitment follows signed moment arms, updating
both cached and underlying synaptic matrices every 10 ms. Twenty-six
compartments switch hip-rotation action across the sampled hip-flexion
sweep. Finite perturbations confirm that positive coordinates named
hip_adduction in this conversion move both feet laterally (abduction).
The anatomical x-anterior/y-superior/z-right frame is explicit; legacy
world placement uses a documented adapter and unchanged baseline physics.

IK/movement PCA and activation NMF are kept separate. Rank-six held-out
activation reconstruction VAF is 0.760; IK-only activation prediction is
0.608. The paired record is only two seconds of quasi-static IK/GRF
back-solved activation, not EMG, and retains a subject/model scale mismatch.
These are mathematical engineering seeds, not identified biological modules.
Figures: `motor_primitives_vs_synergies.png` and `hip_rotation_routing.png`.

Bounded air tests include axis-only, biarticular-only, combined, gain/time-
constant variations, and RG-phase-driven temporal waveforms. A finite
unilateral initialization pulse produces antiphase RG activity, but no
candidate meets the repeated physical knee-alternation gate. Smaller knee
excursion alone is not accepted as improved gait. Added combined trials
attempt afferented air, supported standing, and ground locomotion; these are
staged diagnostics, not a completed five-stage optimization. Standing
completes eight seconds with foot-normal-force/body-weight approximately
0.575 and inferred residual support 0.425; these are force-based proxies,
not measured harness forces. Maximum torso tilt is 30.2 degrees, so this
is not accepted as satisfactory balance. Ground motion has 6/0 detected
right/left cycles and score -246.01, worse than the matched baseline's
7/5 cycles and -174.23. Neither establishes successful forward walking.

Both new neural construction checks have zero reference/accelerated voltage
difference over 200 fixed-mechanics steps. The preceding audit's production-
file statement applies to that audit: this later experiment adds an optional
three-line model-binding hook in runner.py. Existing families have no hook;
the original s3k score remains exactly -161.56754173676563. No production
builder, winner, database, dissertation chapter, or deck is modified by this
functional PF experiment. Full dynamic-controller backend parity is not claimed.

Robot launchers are now in `Code\MuJoCo_SNS\robot_torso\`:
`Launch_AH001_Replay.cmd` (double-click on easteregg2) and
`Launch_AH001_Replay.ps1` (`-RunNew` to simulate before replay).
They launch the saved supported hybrid MuJoCo demonstration, not Simscape.

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
