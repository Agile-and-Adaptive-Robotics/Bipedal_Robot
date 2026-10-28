# AH001 robot with neural control

**Current diagnostic fails leg-clearance checks.** In the verified 65% support
replay, ankle left/right order is reversed in 91.8% of post-ramp samples.
Median foot headings are approximately 28° right and -32° left relative to
the pelvis anterior axis. Leg mesh self-collision is disabled. These runs
establish numerical replay and assisted survival only; they are not valid
stepping or hardware-capability demonstrations. See `foot_heading_audit.json`.
Current tensions are ideal reference actuators; optimized BPA hardware and
its pressure/contraction/compliance behavior have not been imported.

The current MuJoCo 2.3.7 model is `ah001_neural.xml`. It uses **x anterior,
y superior, z right**, both spherical hips, both four-bar knees, welded CAD
feet, and a **ball-and-socket lumbar joint** with a rigid 0.25 m spine bar.
The current CAD export includes Ben's corrected left-knee pose. Source CAD
was not saved by this pipeline.

## Geometry, mass, and routes

The model has 17 independent degrees of freedom: six at the floating pelvis,
six at the hips, two constrained knee motions, and three at the lumbar ball.
Two knee-loop constraints have combined rank four. The companion URDF omits
the cut knee edges; use the MJCF for dynamics.

Total mass is 75 kg: 5.356664 kg for the lower body including welded feet,
0.12 kg for the provisional solid-Onyx spine, and 69.523336 kg balancing
payload. Half the payload attaches to the torso and half to the pelvis.
The torso payload center is 100 mm anterior and 75 mm superior to the lumbar
joint. The adjustable joint location is pelvis-frame (-160, 100, 0) mm,
an upper-rear-center approximation. Material densities and hardware ownership
are engineering assumptions pending physical measurements.

Routes come directly from the right side of
`Solid_Models/OpenSim/Gait2392_Robotbody/gait2392_robotbody.osim`, mirrored
left. The existing Routing27 actuator map selects **27 per side, 54 total**.
The complete source-route audit model, `ah001_routes_audit.xml`, retains
92 route definitions for comparison; these are not 92 robot actuators.

Edit `robotbody_routes.json` for attachment changes. Recorded mirrored
corrections in `route_sign_corrections.json` eliminate the reversed neutral
actions. No new hip-muscle consolidation was performed.

`dynamic_routes.py` evaluates moving-point splines and conditional membership
at runtime and includes moving-guide terms in length velocity and force.
A generic XML-only viewer cannot reproduce those dynamic routes. Use the
Python runner or replay script. The verification covers 492 source function
values, 25 native conditional checks, and 75 constrained poses, with maximum
full moment-arm error 1.38e-9 m. Unselected semitendinosus conditional routes
retain flagged small length jumps after retargeting.

The 54 actuators apply ideal reference tensions with a 30 ms activation
filter and MIF caps from the existing actuator map, including 2500 N per
erector spinae. This is not a calibrated pneumatic actuator simulation.
Twenty ankle/toe actuators lose joint leverage because the feet are welded.

## Current controller and evidence

`run_neural_robot.py` couples a 262-neuron, 230-input SNS network to the robot.
SNS membrane outputs drive the actuators. A bounded postural torque allocator
supplies descending motor-neuron currents; neural RG activity supplies
stepping references. Ia/II/Ib and foot-contact inputs are present. Motor target
current compensation makes this a hybrid engineering controller.

The recorded 20-second supported crouched stance is stable over its final
five seconds: knee excursions 1.35/0.95 degrees, pelvis tilt below 3.92 degrees,
and lumbar rotation below 0.59 degrees. Actual vertical support averages
94.68% body weight, and angular pelvis stabilization remains active.
This is assisted posture control, not independent balance. Disabling all
muscle forces under the same rig settings causes lumbar failure at 0.202 s.

Reference and exact accelerated neural backends produce bit-identical full
20-second stance logs (`controller_assessment.json`). Production human and
W2L baselines remain separate and unchanged.

The drive-3 stepping attempt produces alternating bilateral neural bursts,
but foot clearance is below 1 mm and sustained swing gates fail.
Its forward drift must not be described as walking. The current support
campaign tests explicit stance-load compensation and hard caps on **actual**
vertical support; lowering a nominal support number alone is insufficient.

The newer reduced-support checkpoint completes 40 seconds with a **5% minimum
allocated activation** and a hard **65% actual vertical support cap** after
15 seconds. It includes RG-timed lateral pelvis guidance and angular support.
Fresh reference and exact-fast runs are bit-identical for the complete 40 s
(`supported65_backend_verification.json`). Post-ramp foot clearances reach
55.2/21.8 mm, but lifting remains intermittent and physical cycles skip neural
cycles. This is reduced-support control, not regular stepping or forward walking.
The activation floor is adjustable with `--activation-floor`; zero remains the
general runner default, and the support-65 preset uses 0.05. A minimum activation
is not a coactivation index or a measured PAM stiffness.

Use `run_support65.ps1` on easteregg2 to reproduce this checkpoint with a new
timestamped output. Its `-ActivationFloor` and `-RgTauH` options permit explicit
comparisons. `robot_support65.gif` replays the verified post-ramp interval.

## 28 October: uncrossed defense diagnostic and attachment sensitivity

The earlier crossed-leg support-65 record is not accepted as physical stepping.
The separate `ah001_abduction_diagnostic_v2.xml` trial in
`runs/abd_defense_exact.*` maintains normal ankle ordering through 40 s, with
80% actual vertical assistance after 15 s. Its matched abductor-control removal
reverses ankle order in 67% of the 15–40 s samples. The diagnostic includes
provisional passive knee-stop surrogates, tightened knee-loop closure, and
opposite-leg mesh collisions. These assumptions are not measured hardware.
The corrected knee tangent preserves all hip rotation authority.

Fresh reference/exact-fast runs are bit-identical (`abduction_defense_verification.json`).
The final recorded GIF, still, placement recommendation, and slide verification
are in `Documentation/Reports and Papers/Dissertation/Defense/slideshow/AH001_20261028/`
relative to the repository root. They demonstrate supported hybrid control,
not calibrated BPA behavior, regular stepping, or independent balance.

Route sites are body-fixed virtual attachment/guide points, with the source
moving and conditional functions retained. Sites are not automatically snapped
to CAD surfaces. Offset sites assume ideal rigid brackets, with no bracket
mass or compliance. The foot weld is also provisional: it aligns the ankle
assembly top with the distal shank, not an audited mating datum.
`probe_foot_mount_w2lvar.py` preserves the baseline and tests a hypothetical
60 mm proximal foot translation. It changes neither tibial mesh nor tibial
inertia. Its w2lvar option tests timing inside the same hybrid controller;
the MN clamp still cancels direct PF/reflex contributions. Do not interpret it
as a comparison of full neural muscle-pattern control architectures.

## Run on easteregg2

Use PowerShell and the machine's explicit environment:

```powershell
Set-Location -LiteralPath 'D:\GitHub\Bipedal_Robot'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\verify_ah001.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\verify_dynamic_routes.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\replay_neural_robot.py' --view
```

The viewer replays recorded physics; it does not run a new simulation.
`neural_stance.gif` provides the same recorded motion without a viewer.
The runner saves all options, effective neural gains, source hashes, actual
rig forces, contacts, route states, and neural traces in `runs/`.
Use `assess_supported_step.py` to assess rhythm, sustained contacts,
geometry-based foot clearance, and progression separately.

Earlier `opensim_routes.json`, `source_route_inventory.json`,
`ah001_with_opensim_routes.xml`, and z-up hinge previews are superseded.
Do not use them as the current robot model. Current campaign evidence and
machine-specific ZCode prompts are in
`../spinal/reports_20261027_robot_control/`.
