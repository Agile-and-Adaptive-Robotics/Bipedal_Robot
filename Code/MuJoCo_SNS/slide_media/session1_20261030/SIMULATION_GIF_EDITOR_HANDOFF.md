# Simulation GIF handoff for dissertation and slide-deck editors

**Machine: easteregg2**  
**Prepared by: ChatGPT/Codex, Session 1**  
**Date: 30 October 2026**  
**Repository: `D:\GitHub\Bipedal_Robot`**  
**Python: `D:\Anaconda\envs\myo\python.exe`**

This document identifies the three supplied GIFs, their programs and parameter
records, the tests being performed, and the supported interpretation of each.
The GIFs were copied without re-encoding; original captions and timing remain.
No dissertation text or slide deck was edited while preparing this handoff.

## 1. Current human Gait2392: supported ground stepping

**Full GIF path:**

```text
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\slide_media\session1_20261030\Gait2392_supported_stepping_64pct.gif
```

[Open GIF](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/slide_media/session1_20261030/Gait2392_supported_stepping_64pct.gif)

### Test and interpretation

An isolated **corrected curriculum stage 4** test: harness-supported ground
stepping in place. This is the human Gait2392 model, not the CAD robot.
It uses the experimental nonspiking biw2lvar-Hunt-PT controller, with bilaterally
tied learned conductances, four-cell RF/hamstring PF circuitry, separate
sartorius/gracilis pathways, and an engineering Peterka/Ting-inspired task
balance pathway. Motor-neuron outputs drive the native muscle actuators.
This run does not use the CAD robot's runtime posture/force allocator.

Measured vertical assistance is **63.7263% BW**, despite a **50% BW prescribed
upward force**: the vertical spring/damper supplies additional net assistance.
The feet therefore bear approximately 36% BW on average over the same window.
Horizontal, rotational, and out-of-plane restraints remain active.

The run passes the movement, once-counted ground-contact, frontal-separation,
sustained foot-mesh-clearance, and no-hyperextension gates. Five right and four
left qualifying swings occur in the 5–13 s analysis window. Maximum knee
extension angles over the entire recorded run are **−1.201° right and −1.460°
left**, with 0° neutral and negative angles representing flexion. A fresh
process repeats the recorded positions, controls, and neural traces exactly.
The same saved neural configuration passes air stages 1 and 2. Standing
stage 3 still fails; the full five-stage curriculum is incomplete. This is
not evidence of independent walking or successful further support weaning.

The starting ankle angles and body height are aligned using repaired model
geometry before constructing harness references. The input ankle/height values
below are pre-alignment values, not the final aligned pose. Existing ankle
spring rests are also rebased under the runner's current reference convention.
The earlier 43.61%-assistance result used deeply penetrating initial feet and
is superseded; do not substitute that number for this GIF's 63.73%.

Suggested caption: **“Human Gait2392 model performing supported stepping in
place under the experimental SNS controller, with approximately 64% measured
vertical assistance and protective knee-extension stops.”**

### Programs and evidence, all on easteregg2

- Simulation wrapper: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py`
- Underlying MuJoCo runner: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\runner.py`
- SNS controller: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\controller.py`
- Accelerated neural calculations: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\sparse_pairwise_eval.py`
- Recorded-motion GIF renderer: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\render_ground_preview.py`
- Initial foot alignment: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\ground_alignment.py`
- Saved run configuration: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\contact_feasible_start_v1\ground_guide.json`
- Saved results: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\contact_feasible_start_v1\ground_guide.result.json`
- Recorded trajectory: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\contact_feasible_start_v1\ground_guide.npz`
- Rendered-frame provenance: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\contact_feasible_start_v1\preview\ground_guide\provenance.json`
- Fresh-repeat checkpoint: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\aligned_support_ladder_v1\CHECKPOINT.json`
- Frozen weights/force-profile configuration: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\aligned_support_ladder_v1\packaged_configuration.json`

### Parameters

MuJoCo 2.3.7 with implicitfast; 2 ms physics and neural steps; accelerated
sparse/pairwise neural backend. SNS-Toolbox 1.5.2 is the installed network
construction library. The existing converted human model uses native world
x anterior, y left, z superior, with an explicit anatomical frame adapter.
Do not label those raw coordinates as the CAD robot's y-up frame.

This is the rigid-tendon muscle-surrogate run. The new compliant-tendon
prototype is separate and was not used in this GIF. Shared static force-law
v3 parameters cover 18 bilateral families/36 muscles, not all 92 muscles.

Key parameters: 1.2 shared tonic RG current; sensory timing scale .01;
crossed hindfoot-to-contralateral-RG-F gain 6; landing-event current amplitude
2; biarticular-pair inhibition .3; added E/E and F/F commissural gains zero.
PT gains kp .5, kd .2, ka 0, delay .1 s, bilateral recruitment enabled.
The saved circuit supplies direct joint/ankle PF gains; these are included
in the circuit-options record below rather than inferred from defaults.

Protective knee limits are −90° to 0°, with a 2° pre-neutral margin,
solref [.004, 1], and solimp [.999, .9999, .001, .5, 2]. No position or
velocity clipping is used. Any positive knee angle rejects a run.
Vertical guide stiffness is 5,000 N/m, with damping 1,260 N·s/m under the
baseline scaling. Its reference lowers 40 mm during 2–4 s. Restriction=1
adds 800 N·m/rad and 40 N·m·s/rad to the declared out-of-plane coordinates;
baseline rig restraints also remain. These are training assists.

**Full saved configuration:**

```json
{
  "biarticular_pair_gain": 0.3,
  "commissural_gain": 0.0,
  "learn_reciprocal": true,
  "weights_directory": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\hunt_training\\runs\\symmetric_bothcadence_20261030_s1_v2\\tonic_0p02_exact\\calibrated_pf_drive_v1\\knee_0p16\\weights",
  "backend": "sparse",
  "diagnose_swing": true,
  "stage": 4,
  "restriction": 1.0,
  "height_offset": -0.02,
  "start_hip_deg": -7.0,
  "start_knee_deg": -10.0,
  "start_ankle_deg": 0.0,
  "balance": {
    "kp": 0.5,
    "kd": 0.2,
    "ka": 0.0,
    "delay_s": 0.1,
    "replace_inherited": false,
    "bilateral_recruitment": true
  },
  "force_controlled_support": false,
  "support_fraction": 0.5,
  "regional_contacts": true,
  "cross_contact_gain": 6.0,
  "cross_contact_event_amplitude": 2.0,
  "diagnostic_scope": "Symmetric swing-expression test; early stages require revalidation.",
  "force_law_parameters_file": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\physics_fidelity\\shared_force_law_functional_v3.json",
  "timing_afferents": true,
  "timing_scale": 0.01,
  "knee_training_limits_deg": [
    -90.0,
    0.0
  ],
  "commissural_e_gain": 0.0,
  "standing_solver": "coupled_symmetric",
  "vertical_damping_multiplier": 4.2,
  "harness_lowering_m": -0.04,
  "harness_lowering_start_s": 2.0,
  "harness_lowering_duration_s": 2.0,
  "vertical_stiffness_scale": 0.01,
  "rg_tonic_current": 1.2,
  "align_ground_pose": true
}
```

**Saved circuit options:**

```json
{
  "require_bilateral_symmetry": true,
  "direct_ankle_pf_gain": 0.5238095238095238,
  "direct_joint_pf_gain": 0.5238095238095238
}
```

**Complete learned-parameter files:**

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\calibrated_pf_drive_v1\knee_0p16\weights\timing_weights.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\calibrated_pf_drive_v1\knee_0p16\weights\pf_weights.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\calibrated_pf_drive_v1\knee_0p16\weights\bi_weights.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\calibrated_pf_drive_v1\knee_0p16\weights\mn_weights.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\calibrated_pf_drive_v1\knee_0p16\weights\reflex_weights.json`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\symmetric_bothcadence_20261030_s1_v2\tonic_0p02_exact\calibrated_pf_drive_v1\knee_0p16\weights\circuit_options.json`

**Recorded alignment details:**

```json
{
  "scope": "Explicit initial pose alignment; not a hold/balance/walking result.",
  "flatten_ankles": true,
  "target_gap_m": -0.001,
  "initial_region_gaps_m": [
    -0.03492652120707534,
    -0.05624129280601519,
    -0.03636609198660396,
    -0.0576336320895777
  ],
  "final_region_gaps_m": [
    0.00044187236275077663,
    0.0004370247659194695,
    -0.0009999999999999937,
    -0.0009951524031686762
  ],
  "original_native_ankle_deg": [
    0.0,
    0.0
  ],
  "common_native_ankle_deg": 5.621369783258007,
  "vertical_joint": "pelvis_ty",
  "vertical_qpos_address": 1,
  "vertical_coordinate_shift_m": 0.03890394692036618,
  "world_normal_displacement_m": 0.03890394611365443,
  "rig_reference_contract": "Construct/capture the support rig AFTER this pose. If aligning an already-built rig, shift its vertical qpos_spring by exactly vertical_coordinate_shift_m too; do not silently retain the old reference.",
  "pipeline": "Native q values -> repaired probe -> aligned native q values copied back -> capture -> final repaired rig/ref construction.",
  "probe_raw_hip_axis": [
    [
      0.0,
      0.0,
      1.0
    ],
    [
      0.0,
      0.0,
      -1.0
    ]
  ]
}
```

## 2. CAD robot: 65% support-cap hopping-like motion

**Full GIF path:**

```text
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\slide_media\session1_20261030\Robot_supported_hopping_65pct.gif
```

[Open GIF](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/slide_media/session1_20261030/Robot_supported_hopping_65pct.gif)

### Test and interpretation

Historical **reduced-support neural-control diagnostic**, not a completed
five-stage curriculum. The 75 kg CAD robot has 27 reference-tension actuators
per side, 54 total; the feet are welded. Its controller is a **hybrid layered
SNS plus engineering posture/force allocator**: 262 neurons and 230 neural
inputs, RG-driven movement references, sensory inputs, and compensated
motor-neuron target currents. The ideal tension actuators have a 30 ms
activation filter and force caps from the existing muscle map. They do not
constitute a calibrated pneumatic/BPA hardware model.

The simulation lasts 40 s. After a 5 s standing period, the **actual vertical
support cap** reduces to 65% BW over 5–15 s. This GIF shows the post-ramp
15–39.75 s interval. Nominal support .8 is not the post-ramp measured value;
whole-run mean assistance is about 73.37% because it includes startup/weaning.
The 5% floor is a minimum allocated activation, not an antagonist coactivation
index or measured BPA stiffness. Lateral and angular assistance remain active.

The recording completes and the reference/exact-fast runs match. Foot lifting
is intermittent, and physical cycles skip neural cycles. Hopping-like motion
is a descriptive visual label, not a validated hopping gait. These historical
robot clips predate the current strict zero-hyperextension gate; do not claim
that they passed that later constraint.

Suggested caption: **“CAD robot model showing assisted SNS-driven motion with
a 65% vertical-support cap, lateral/angular guidance, and ideal reference-tension
actuators.”** Avoid calling this regular walking or demonstrated BPA performance.

### Programs and parameter evidence

- Simulation: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py`
- Machine-guarded reproduction preset: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_support65.ps1`
- SNS network construction: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\build_network.py`
- Exact accelerated neural evaluator: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\exact_fast_eval.py`
- Moving/conditional route update: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\dynamic_routes.py`
- MuJoCo model: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\ah001_neural.xml`
- Route data: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\robotbody_routes.json`
- Replay/GIF utility: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\replay_neural_robot.py`
- Stepping assessment: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\assess_supported_step.py`
- Recorded trajectory: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\verified65_exact_fast.npz`
- Recorded parameters/results: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\verified65_exact_fast.json`
- Assessment: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\verified65_exact_fast_assessment.json`
- Backend repeat verification: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\supported65_backend_verification.json`
- Archived runner used: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\source_snapshots\70ecee4442ea02b5e7e8aa2a448685bb624abbfc5cf98af4ae0b0214d89f34a1.py`

### Parameters

MuJoCo implicitfast, 2 ms steps, native robot anatomical frame x anterior,
y superior, z right, gravity [0, −9.81, 0]. Total model mass 75 kg;
ball-and-socket lumbar joint; rigid spine .25 m. Payload is the mass remainder,
not an additional 75 kg. Half the payload is carried through the torso, and
half through the pelvis/ilium support. The reference payload offset is
[.100, .075, 0] m from the lumbar mount in anatomical coordinates.

Important differences from the older 85% clip: event-based pattern,
explicit swing unloading, lateral pelvis guidance, horizontal stiffness
200 N/m, and right knee-reference amplitude .85 rad versus the base .6 rad.
These historical robot reference amplitudes are not the human model's signed
knee angles. Full recorded options follow; any null field means no explicit
override in that archived option record.

```json
{
  "mode": "step",
  "duration": 40.0,
  "stand_time": 5.0,
  "support": 0.8,
  "horizontal_stiffness": 200.0,
  "lateral_shift": 0.03,
  "lateral_stiffness": 1000.0,
  "lateral_damping": 120.0,
  "knee_amplitude": 0.6,
  "knee_amplitude_right": 0.85,
  "knee_amplitude_left": null,
  "hip_flexion": 0.35,
  "hip_extension": -0.12,
  "pattern": "event",
  "unload_swing": true,
  "activation_floor": 0.05,
  "walk_drive": 3.0,
  "hip_kp": 25.0,
  "hip_kd": 5.0,
  "knee_kp": 10.0,
  "knee_kd": 1.5,
  "load_compensation": 1.0,
  "backend": "exact_fast",
  "motor_tau": 0.005,
  "activation_tracking_gain": 0.0,
  "stable_pd_horizon": 0.0,
  "rig_angular_stiffness": 300.0,
  "rig_angular_damping": 60.0,
  "frontal_assist": 0.0,
  "disable_muscles": false,
  "vertical_stiffness": 800.0,
  "support_cap": 1.5,
  "support_cap_final": 0.65,
  "support_wean_start": 5.0,
  "support_wean_duration": 10.0,
  "diagnostic_direct": false,
  "output": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\robot_torso\\runs\\verified65_exact_fast"
}
```

**Effective neural gains saved with this run:**

```json
{
  "rg_mutual_inh": 4.0,
  "rg_adapt_inh": 2.5,
  "rg_to_pf": 2.4,
  "pf_recip_inh": 4.0,
  "pf_to_mn": 2.0,
  "posture_to_mn": 1.0,
  "ia_to_mn": 0.6,
  "ia_to_antagonist": 0.4,
  "ii_to_mn": 0.4,
  "ib_to_mn_inh": 0.35,
  "ib_group_exc": 0.5,
  "ib_exc_to_mn": 0.6,
  "bal_lat_to_abd": 0.8,
  "bal_trunk": 1.0,
  "descend_to_rg_e": 1.7,
  "descend_to_rg_f": 1.4,
  "posture_to_rg_e": 0.8,
  "drive_to_pf": 0.05,
  "posture_gain": 1.0,
  "phase_reset_e": 0.0,
  "phase_reset_f": 0.0,
  "f1_kneext_inh": 0.0,
  "f1_anklepf_inh": 0.0,
  "renshaw": 0.0,
  "ankle_post_walk_trim": 1.0,
  "heel_rge": 0.15,
  "toe_rge": 0.05,
  "ib_rge": 0.0,
  "contact_onset": 0.0,
  "contra_swing": 0.0,
  "contra_kinh": 0.0,
  "full_rules": 0.0,
  "pm_gain": 0.0,
  "pm_T": 1.2,
  "pm_ws": 0.0,
  "pm_add": 0.0,
  "pm_aff": 0.0,
  "ia_in": 0.0,
  "c1_gain": 0.15,
  "v3_gain": 0.0,
  "rg_weak_exc": 0.0,
  "joint_pf": 0.0,
  "heel_pf_layer": 0.0,
  "toe_df_inh": 0.0,
  "heel_in_f_exc": 0.0,
  "ia_pf_f": 0.0,
  "ii_pf_f": 0.0,
  "ib_e_central": 0.0,
  "ia_f_central": 0.0,
  "ii_f_central": 0.0,
  "ii_e_central": 0.0,
  "ia_f_contra_f": 0.0,
  "v3_to_ibexc": 0.0,
  "aff_e_rg": 0.0,
  "aff_f_rg": 0.0,
  "aff_e_pf": 0.0,
  "aff_f_pf": 0.0,
  "vest_ext": 0.0,
  "vest_flex_inh": 0.0,
  "vest_prop": 0.0,
  "coact": 0.0,
  "syn6": 0.0,
  "syn6_brainstem": 0.0,
  "gamma_l": 0.0,
  "gamma_r": 0.0,
  "alpha_l": 0.0,
  "alpha_r": 0.0,
  "alpha_to_mn": 0.5,
  "contra_flex": 0.0,
  "contact_ank_mn": 0.0,
  "iain_pf_gate": 0.0,
  "pf_mixed_recip": 0.0,
  "ii_set_pose": 0.0,
  "gamma_coact": 0.0
}
```

The archived JSON does not serialize every intrinsic neuron time constant.
It records hashes of the imported network/parameter sources. Preserve the
archived runner and parameter record; current preset defaults alone are not
a complete reconstruction of historical neural parameters.

## 3. CAD robot: earlier 85% support-cap motion

**Full GIF path:**

```text
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\slide_media\session1_20261030\Robot_earlier_motion_85pct.gif
```

[Open GIF](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/slide_media/session1_20261030/Robot_earlier_motion_85pct.gif)

### Test and interpretation

Historical **increased-lift/reference-amplitude diagnostic** during the support
and minimum-activation campaign. Same CAD robot and hybrid controller family
as clip 2. The matched record is lift_floor05_cap85: .6 rad knee amplitude,
.35/−.12 rad hip flexion/extension references, .05 activation floor,
load compensation 1, and a final actual-support cap .85. This is the earlier
direct-reference pattern, not the later event-pattern/swing-unloading preset.
The entire 0–24.75 s simulation is shown, including startup and cap reduction
during 5–15 s. Consequently, actual support can exceed 85% early in the GIF.

Provenance qualification: the original GIF creation invocation was not found.
The associated record was identified from its 25 s duration and matching
displayed support at 2.25 s (90.8%), 9.75 s (88.9%), and 17.25 s (85.0%),
distinguishing the nearby saved candidate records. Treat this as a matched
historical recording, not a freshly reproduced/render-verified trial.

The saved assessment reports bilateral neural rhythm and antiphase, but the
bilateral foot-lifting and forward-progression gates fail. This clip is useful
as an example of assisted exploratory motion or a limitation, not evidence
of successful walking. Its knee motion was not certified under the later
zero-hyperextension requirement. No full-curriculum claim applies.

Suggested caption: **“Earlier assisted-motion trial of the CAD robot model
using increased movement references and an 85% final support cap; regular
bilateral stepping was not achieved.”**

### Programs and records

- Simulation program: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py`
- SNS construction: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\build_network.py`
- Accelerated evaluator: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\exact_fast_eval.py`
- Replay utility currently available: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\replay_neural_robot.py`
- Recorded trajectory: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\lift_floor05_cap85.npz`
- Saved options/results: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\lift_floor05_cap85.json`
- Assessment: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\runs\lift_floor05_cap85_assessment.json`
- Original GIF: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\supported_step_85_floor05.gif`

### Full saved options

```json
{
  "mode": "step",
  "duration": 25.0,
  "stand_time": 5.0,
  "support": 0.8,
  "horizontal_stiffness": 20.0,
  "knee_amplitude": 0.6,
  "hip_flexion": 0.35,
  "hip_extension": -0.12,
  "activation_floor": 0.05,
  "walk_drive": 3.0,
  "hip_kp": 25.0,
  "hip_kd": 5.0,
  "knee_kp": 10.0,
  "knee_kd": 1.5,
  "load_compensation": 1.0,
  "backend": "exact_fast",
  "motor_tau": 0.005,
  "rig_angular_stiffness": 300.0,
  "rig_angular_damping": 60.0,
  "frontal_assist": 0.0,
  "disable_muscles": false,
  "vertical_stiffness": 800.0,
  "support_cap": 1.5,
  "support_cap_final": 0.85,
  "support_wean_start": 5.0,
  "support_wean_duration": 10.0,
  "diagnostic_direct": false,
  "output": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\robot_torso\\runs\\lift_floor05_cap85"
}
```

**Effective neural gains saved with this run:**

```json
{
  "rg_mutual_inh": 4.0,
  "rg_adapt_inh": 2.5,
  "rg_to_pf": 2.4,
  "pf_recip_inh": 4.0,
  "pf_to_mn": 2.0,
  "posture_to_mn": 1.0,
  "ia_to_mn": 0.6,
  "ia_to_antagonist": 0.4,
  "ii_to_mn": 0.4,
  "ib_to_mn_inh": 0.35,
  "ib_group_exc": 0.5,
  "ib_exc_to_mn": 0.6,
  "bal_lat_to_abd": 0.8,
  "bal_trunk": 1.0,
  "descend_to_rg_e": 1.7,
  "descend_to_rg_f": 1.4,
  "posture_to_rg_e": 0.8,
  "drive_to_pf": 0.05,
  "posture_gain": 1.0,
  "phase_reset_e": 0.0,
  "phase_reset_f": 0.0,
  "f1_kneext_inh": 0.0,
  "f1_anklepf_inh": 0.0,
  "renshaw": 0.0,
  "ankle_post_walk_trim": 1.0,
  "heel_rge": 0.15,
  "toe_rge": 0.05,
  "ib_rge": 0.0,
  "contact_onset": 0.0,
  "contra_swing": 0.0,
  "contra_kinh": 0.0,
  "full_rules": 0.0,
  "pm_gain": 0.0,
  "pm_T": 1.2,
  "pm_ws": 0.0,
  "pm_add": 0.0,
  "pm_aff": 0.0,
  "ia_in": 0.0,
  "c1_gain": 0.15,
  "v3_gain": 0.0,
  "rg_weak_exc": 0.0,
  "joint_pf": 0.0,
  "heel_pf_layer": 0.0,
  "toe_df_inh": 0.0,
  "heel_in_f_exc": 0.0,
  "ia_pf_f": 0.0,
  "ii_pf_f": 0.0,
  "ib_e_central": 0.0,
  "ia_f_central": 0.0,
  "ii_f_central": 0.0,
  "ii_e_central": 0.0,
  "ia_f_contra_f": 0.0,
  "v3_to_ibexc": 0.0,
  "aff_e_rg": 0.0,
  "aff_f_rg": 0.0,
  "aff_e_pf": 0.0,
  "aff_f_pf": 0.0,
  "vest_ext": 0.0,
  "vest_flex_inh": 0.0,
  "vest_prop": 0.0,
  "coact": 0.0,
  "syn6": 0.0,
  "syn6_brainstem": 0.0,
  "gamma_l": 0.0,
  "gamma_r": 0.0,
  "alpha_l": 0.0,
  "alpha_r": 0.0,
  "alpha_to_mn": 0.5,
  "contra_flex": 0.0,
  "contact_ank_mn": 0.0,
  "iain_pf_gate": 0.0,
  "pf_mixed_recip": 0.0,
  "ii_set_pose": 0.0,
  "gamma_coact": 0.0
}
```

Both robot clips use the same basic 75 kg torso/payload specification above.
Source hashes in the archived records identify their original runner versions;
the current runner/replay utilities have evolved. Do not regenerate old results
with current defaults and assume they are identical.

## Media properties and identity checks

The copied slide assets have the following verified properties. All loop
durations match their GIF frame delays. The Gait2392 clip is 1× playback of
5.0–12.9 s; the robot clips have .25 s frame delays and 1× recorded playback.

```json
[
  {
    "path": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\slide_media\\session1_20261030\\Gait2392_supported_stepping_64pct.gif",
    "width_px": 880,
    "height_px": 680,
    "frames": 80,
    "loop_duration_s": 8.0,
    "sha256": "ce83843cbf052c4080312aa45c1edecd7cbeae445c11997eba9600f625e20cd9"
  },
  {
    "path": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\slide_media\\session1_20261030\\Robot_supported_hopping_65pct.gif",
    "width_px": 800,
    "height_px": 660,
    "frames": 100,
    "loop_duration_s": 25.0,
    "sha256": "df323d39578ebe0d23d42c68dd241e5d5dd8a3bde66cd4e76a41da459215f68b"
  },
  {
    "path": "D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\slide_media\\session1_20261030\\Robot_earlier_motion_85pct.gif",
    "width_px": 800,
    "height_px": 660,
    "frames": 100,
    "loop_duration_s": 25.0,
    "sha256": "71243c129bd495dd9c0927e67ce8c861a5fac0dfed31e4a82e7c005cb498968c"
  }
]
```

## Editorial boundaries

- Use **models**, not physical robot experiments. All three are simulated.
- Distinguish the human model's approximately 64% **measured assistance** from
  the robot clips' **final actual-support caps** of 65% and 85%.
- Keep hybrid allocator terminology for the CAD robot. Do not describe those
  clips as the same controller or actuator physics as the current Gait2392 run.
- All runs retain substantial support/restraints. No harness-free walking,
  completed five-stage curriculum, validated hopping gait, or measured BPA
  locomotor performance is established by these GIFs.
- Clip 1 is the current protected-knee ground illustration; clips 2 and 3 are
  historical assisted-motion diagnostics. Preserve those different test statuses.

Further current evidence:

- [TRAINING_RECORD.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/TRAINING_RECORD.md)
- [SESSION1_PHYSICS_GROUND_CHECKPOINT_20261030.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/SESSION1_PHYSICS_GROUND_CHECKPOINT_20261030.md)
- [CONTACT_INITIALIZER_AND_ALIGNMENT.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/physics_fidelity/CONTACT_INITIALIZER_AND_ALIGNMENT.md)
- [README.md](D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/README.md)
