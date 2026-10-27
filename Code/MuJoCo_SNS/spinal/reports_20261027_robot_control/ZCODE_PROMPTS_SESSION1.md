# ZCode prompts from Session 1 on easteregg2

These prompts are for **easteregg2**, the heavy workstation currently running
this session. They are not instructions for EB475WS4 or the laptop. All paths
below are absolute Windows paths on easteregg2.

Common environment: repository `D:\GitHub\Bipedal_Robot`, simulation Python
`D:\Anaconda\envs\myo\python.exe`, MuJoCo 2.3.7, and neural/physics step 0.002 s
for the new CAD robot. The legacy W2L runner retains its original 1 ms physics
with a neural update every two substeps. Preserve the dirty working tree and
Session 2's report edits. Do not commit or push.

## Robot neural control continuation

Continue the AH001 robot neural-control work on **easteregg2**. First read:

1. `D:\GitHub\Bipedal_Robot\AGENTS.md`
2. `D:\GitHub\Bipedal_Robot\CHATGPT_REPORT.md`, specifically the block delimited
   `SESSION_1_SIMULATION_RECOVERY_20261027_START` and its matching end marker.
   Other sessions are editing this shared report. Change only your own bounded
   section or append a clearly labeled update; never rewrite it from an old copy.
3. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\README.md`
4. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\controller_assessment.json`
5. `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\dynamic_route_verification.json`

Use these current files:

- Model: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\ah001_neural.xml`
- Runner: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py`
- Editable routes: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\robotbody_routes.json`
- Dynamic geometry: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\dynamic_routes.py`
- Source inventory: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\robotbody_source_inventory.json`
- Sign corrections: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\route_sign_corrections.json`
- Original route source: `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\gait2392_robotbody.osim`
- Existing actuator map: `D:\GitHub\Bipedal_Robot\Code\Matlab\Mesh_Optimization\Routing27_v2\actuator_map_27.json`
- Current CAD/foot inventories: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_robot_control\`

The robot uses **x anterior, y superior, z right**, y-up gravity, both four-bar
knees, spherical hips, a lumbar **ball-and-socket joint**, and a rigid 250 mm
spine. The pelvis datum removes arbitrary global CAD placement. The refreshed
snapshot includes Ben's left-knee correction. Do not reuse the superseded z-up
files `opensim_routes.json` or `ah001_with_opensim_routes.xml` in that directory.

The selected 27 right-side actuator routes come directly from robotbody and
are mirrored left, giving **54 actuators**. The other source routes are retained
for separate audits, not additional actuators. Ben deferred further hip-muscle
consolidation. Do not reinterpret appliesForce flags or activate all 92 muscles.

The resolved foot source is
`D:\GitHub\Bipedal_Robot\Solid_Models\Prospective_Parts\Ty Ankle Foot\Foot and Ankle Design, Newest (7-17-20)\Foot Ankle.SLDASM`.
All 22 components resolved during this session. All internal foot joints and
ankle connections are currently welded. Therefore ankle-only actuators may
have zero joint leverage even though the feet now provide ground contact.

The model totals 75 kg. Payload is the remainder after CAD and spine masses,
split equally between torso and pelvis/ilium. Its torso center is 100 mm
anterior and 75 mm above the lumbar joint. Erector reference MIF is 2500 N per
side; the full 27-actuator design capacities are in the map above. Current
actuators are first-order reference tensions, not calibrated BPA force models.
Do not label their capacities as measured hardware performance.

Moving paths use exact source cubic functions. Conditional membership changes
at runtime. Preserve the extra guide-motion Jacobian in force and velocity
calculations; simply moving MuJoCo sites omits this derivative. Source comparisons
and constrained finite differences must pass after any attachment change.
Five previously flagged neutral actions are corrected; the refreshed audit
must continue to report zero reversals. Attachment shifts are explicitly logged.

The controller is a hybrid layered SNS with descending force allocation and
compensated MN currents. Its construction and numerical execution do not by
themselves establish autonomous spinal control. The verified standing comparison
uses substantial vertical and angular pelvis assistance. Report actual support
forces and support torques separately. Do not count the direct-force diagnostic,
muscle-off ablation, initial settling displacement, contact chatter, or a long
duration as successful neural walking.

Reproducible PowerShell commands on easteregg2:

```powershell
if ((hostname).Trim().ToLowerInvariant() -ne 'easteregg2') { throw 'Wrong machine' }
Set-Location -LiteralPath 'D:\GitHub\Bipedal_Robot'
$simulationPython = 'D:\Anaconda\envs\myo\python.exe'
& $simulationPython 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\verify_dynamic_routes.py'
& $simulationPython 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py' --mode stand --duration 20 --backend exact_fast --output 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\zcode_stand_20s'
& $simulationPython 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py' --mode step --duration 20 --stand-time 5 --backend exact_fast --output 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\zcode_step_20s'
```

Review the JSON and raw trajectories generated by each command. Establish a
stable supported posture, test muscle-off and perturbation controls, then obtain
repeatable bilateral knee and foot cycles with meaningful clearance and contact
timing. Separate supported marching from forward locomotion and harness-free
walking. Use staged neural-only rhythm tests before tuning full robot trials.
Do not optimize another survival-only scalar that rewards collapse or rigid
support. Save every tuned parameter and its loader, use independent processes,
and re-check backend trajectory equivalence after changes.

## Di Russo original baseline and gait2392 adaptation

Continue Di Russo replication on **easteregg2**, with repository
`D:\GitHub\Bipedal_Robot` and Python `D:\Anaconda\envs\myo\python.exe`.
This is the separate task that remains unreplicated; do not treat robot-control
results as evidence for the paper.

Read the exact recovery report:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_codex\REPORT.md`.
The original ZCode report reviewed by this session is outside the checkout at
`D:\GitHub\_zcode_backups\diss_campaign_20261024\SESSION_REPORT_20261026.md`.

Original author sources:

- `https://github.com/DiRussoAndrea/Spinal_controller`
- `https://gitlab.com/simgait/SCONE`

Staged exact locations:

- Scenario: `D:\GitHub\Bipedal_Robot\tmp\codex_dirusso_original\spinal_controller_main.scone`
- Parameters: `D:\GitHub\Bipedal_Robot\tmp\codex_dirusso_original\import\initializations\spinal.par`
- Model: `D:\GitHub\Bipedal_Robot\tmp\codex_dirusso_original\import\models\gait9dof18musc_Thelen_20170320.osim`
- Initial state: `D:\GitHub\Bipedal_Robot\tmp\codex_dirusso_original\import\initializations\f0918_init_walk_20170205.sto`
- Custom fork: `D:\GitHub\Bipedal_Robot\tmp\codex_simgait_scone`
- Custom controller: `D:\GitHub\Bipedal_Robot\tmp\codex_simgait_scone\src\sconelib\scone\controllers\SegmentalController.cpp`
- Fork build guidance: `D:\GitHub\Bipedal_Robot\tmp\codex_simgait_scone\README.md`
- Parameter inventory: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_codex\author_parameter_inventory.json`
- Source-model check: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_codex\inspect_author_model.py`

The custom fork revision used is
`58ea68762791e133b0b18f2f0db83e3eb2a18fe3`. Its documented supported build
route is Ubuntu 18.04; Windows builds require additional work. Do not invent a
working SCONE executable path: the original SegmentalController runtime has
not been installed or replayed on easteregg2. The sparse checkouts may need
their remaining tracked files fetched before a build.

First reproduce the optimized native scenario on its original 9-DoF, 18-muscle
model. Record the executable, dependency versions, initial state, parameters,
objective terms, support forces, and trajectories. OpenSim 4.4.1 loading the
model is already verified; that does not establish controller compatibility.

Audit all 256 parameters and the actual custom controller, including delayed
state histories and update ordering, muscle-specific delays, fiber-based
normalization, output saturation, antagonist-IaIN Renshaw wiring, balance gates,
and heel reset semantics. Use original neural/actuator traces as regression
targets before adapting to gait2392. The neumovelab repository is a comparative
reimplementation, not the original author source.

The adaptation at
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\dirusso_gait2392\`
does not yet reproduce the original controller. Its rejected search is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_codex\optimize_dirusso.py`.
It requires `--debug-surrogate` deliberately; its duration/speed objective chose
a somersault. Do not resume it as a valid gait optimization. Reproduce the native
baseline, port the full objective and equations, then validate pose, contacts,
progression, fresh-process replay, and a perturbation basin. If the native
runtime cannot be obtained, report that exact prerequisite rather than hiding
it behind another guessed-controller optimization.

## Existing human and robot families

Ben confirmed the later campaign means `s3k`, `w2lvar`, and `syn6`. Work remains
on easteregg2 in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal` with
`D:\Anaconda\envs\myo\python.exe`. Read
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\DESIGN.md` before changing
the runner, builders, or parameters. Follow Ben's current task sequencing.

Preserve the production default census and the documented bit-exact anchor.
New robot controller networks have their own census and must not replace those
contracts. Verify knob delivery, dynamic action signs, per-term metrics,
sensitivity, fresh-process basins, and support forces for every model/family.
Never overwrite an existing study with incompatible distributions or schemas.
Do not report suspended cycling as harness-free walking.
