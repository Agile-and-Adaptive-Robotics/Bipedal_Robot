# Contact handover repair and harness weaning

Machine: **easteregg2**. Session 1, 3 November 2026. Selected condition: **clearance77p5_damped**. No manuscript or Overleaf edits.

## Result and acceptance

The selected CAD model completes 28 s. Actual vertical harness support is **77.5% BW during 16–28 s**, after an 80%-to 77.5% cap ramp during 8–16 s. This is a 2.5 percentage-point reduction from the previous 80% selection. The feet bear approximately 22.5% BW: mean right/left normal forces 78.58/86.97 N. The whole-trial mean includes the initial high-support period and must not be substituted for the final-rung value.

There are **7 right/6 left** qualifying unloaded whole-mesh-clear intervals and **12 alternations** in the final 12 s. Acceptance requires whole-foot mesh clearance>5 mm, ipsilateral normal load<7.36 N, contralateral load>15 N, and an interval>=80 ms, with at least 3 intervals per side and 5 alternations. Merely surviving or oscillating does not pass. Knee checks cover every physics step, not just the 10 ms trace.

Minimum knee flexion is **5.10 degrees**. Peak whole-mesh clearances are 35.40/53.54 mm right/left. Minimum sampled lateral world-z mesh gap is 47.23 mm. Period is 1.86 s. The 28 s world forward displacement is -74.7 mm: guided stepping approximately in place.

Harness attitude/lateral assistance remains active. Mean absolute x/y/z rig moments in the final window are 10.80/0.93/110.08 Nm. Horizontal x/z RMS forces are 1.35/2.70 N. Reducing vertical support did not remove this angular guidance. Retained soft contacts produce transient mesh penetration of 16.13/12.91 mm right/left. No contact coefficients or physical joint/muscle/tendon model were changed.

## Contact handover and muscle-mediated clearance

The old binary contact rule assigned about 73.6 N planned support as left contact first reached 1.25 N. Its next phase switch assigned the full 147 N despite only 30.7 N measured load. The old left knee crossed extension at 6.71 s. That exact saved condition now completes 18 s under the new handover, with 6 R/7 L clear unloaded intervals and at least 2.53 degrees per-physics-step knee flexion.

For each leg, define u=clip((F-F0)/(F1-F0),0,1), c=u^2*(3-2u). F0 is the saved contact threshold, nominal 1 N; F1=5 N. Phase preferences p are normalized **before** multiplying by c. Requested support is L=max(W-F_harness,0)*p*c. At each 10 ms allocation update, L_filtered+=alpha*(L-L_filtered), alpha=1-exp(-dt/0.010 s). A current-contact ceiling L_filtered<=max(W-F_harness,0)*c removes unsupported feedforward as contact vanishes. Confidence is never renormalized into full load. No extra support force is applied by this computation.

A smooth swing-clearance task restores right-foot unloading on longer/lower-support trials. The target is 0.04*sin(pi*u_swing)^2 metres, using the existing RG-event-derived swing progress. The task uses the actual minimum height of every foot mesh. Its virtual force is window*clip(400*(target-height)-24*vertical_speed,0,50) N, with window=sin(pi*u_swing)^2. The vertical Jacobian maps that task into hip/knee muscle torque targets; root and distal foot DOFs are excluded. The existing bounded allocator and SNS MN outputs deliver the command. This is a virtual task force, not an externally applied foot force or a kinematic foot-position correction. The label_damped refers to velocity feedback in this task, **not a new physical damping coefficient**.

The final controller retains 262graded cells, 230inputs, and 54reference-tension actuators. SNS RG activity sets timing; the engineering allocator supplies posture, stance-load, and clearance MN target currents. Force output retains the 30 ms activation filter and mirrored 27-per-side robotbody roster. KompareXXVII selected pressure-actuator sizing/maps are not used in this experiment. This is not a new Hunt or Rustuna fitting campaign. All scalar policies are applied bilaterally; geometric and measured state differences remain.

The predictive guard stops before one degree from full extension, using current angle plus one physics-step positive velocity. It also records maximum knee angle after every integration substep. Failed guard runs are rejected, without resetting or clamping qpos. Every retained completed body trial remained flexed, including rejected cases; process-start logging errors are separate from scientific trials.

## Verification and scope

Fresh repeat exact: **True**. Reference-backend exact: **True**. Differences cover all saved arrays, including contact, actuator, neural, and clearance traces. The 1 ms test keeps neural dt 2 ms and allocator/log dt 10 ms; it tests integration sensitivity rather than exact trajectory equality.

Local ±1% checks: **5/8 completed pass; expected 8, status complete**. These perturb 22shared legacy controller scalars, six handover/clearance settings, and three passive families, with bilateral homologues tied. Fixed neural matrices, harness protocol, and the safety bound are outside this perturbation scope. This is finite local verification, not a global stability proof. Full per-case outcomes are in FINAL_VERIFICATION.json. The exact formerly failed 80% case is also retained separately; that result must not be conflated with a local perturbation at the new 77.5% condition.

## Alternatives retained

A 40 N/50 ms contact transition delayed support and caused flexion collapse. At5 N,10 ms transfer preserved6R/7L nominal 18 s stepping;25 ms transfer also prevented the old failure but weakened sustained right unloading. Lower75/70/60% rungs reached the predictive guard. Changing stance flexion, lowering harness height, and trial flexor braking did not produce an accepted lower-support condition. A measured-force-plus-allowance ceiling caused crouch/height failures. A 77.5% run without the clearance task survived but had0R/6L qualifying intervals, so it was rejected as dragging. The600 N/m clearance setting passed nominal77.5% but had only1.17degree knee margin;400 N/m plus24 N s/m velocity feedback was selected for5.10degree margin and bilateral clearance.75% with either clearance setting still fails. Predictive flexor braking and the measured-force ceiling are **not enabled** in the selected controller.

Four initial brake probes stopped before integration because a logging-array construction error was caught. Their console/status files remain; corrected_v2 trials were rerun. None was promoted as successful evidence.

## Mechanics, execution, and media

The anatomical toe-base model and routes remain exactly those of the prior anatomical revision:75 kg,23independent mechanical DoF, fixed midfoot, collective toe hinge, and retained joint coefficients. Prior geometry/virtual-work checks and definitions are linked from the earlier EXPERIMENT_RECORD. No new compliant tendon or physical damping option is adopted.

Execution uses D:\Anaconda\envs\myo\python.exe, MuJoCo2.3.7, SNS-Toolbox1.5.2 exact_fast, Numba dynamic routes, bounded allocation, and recorded-state MuJoCo/Pillow rendering. Selected 28 s took52.78 s with concurrent jobs. Fresh/reference elapsed times use concurrent instrumentation and are not presented as a controlled speed ratio. The earlier paired Gait benchmark remains 6.82x full-process/18.31x neural-call acceleration in its own declared scope.

The GIF uses 20–28 s at 12 fps, 96frames; it shows the final77.5% rung at recorded cadence. The print still is from 22.9 s. Camera/material/route-highlighting changes affect display only. No pose, activation, playback speed, or physics correction is made for the GIF. Editable controller and frozen-source hashes are retained. All trials and sources are in this owned external folder, rather than disposable repo tmp.


## Completed body trials

| Condition | Duration s | Final-window gate | Clear R/L | Outcome |
|---|---:|---|---|---|
| brake60_g8_v2 | 12.810 | fail | [0, 0] | predictive knee flexion safety margin |
| brake70_g4_v2 | 16.380 | fail | [0, 0] | predictive knee flexion safety margin |
| brake70_g8_v2 | 15.502 | fail | [0, 0] | predictive knee flexion safety margin |
| brake80_g4_v2 | 28.000 | fail | [0, 6] | completed |
| clearance70 | 12.230 | fail | [0, 0] | predictive knee flexion safety margin |
| clearance75 | 17.238 | fail | [1, 1] | predictive knee flexion safety margin |
| clearance75_damped | 16.552 | fail | [1, 0] | predictive knee flexion safety margin |
| clearance77p5 | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_basin00 | 18.880 | fail | [2, 1] | predictive knee flexion safety margin |
| clearance77p5_damped_basin01 | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_basin02 | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_basin03 | 24.466 | fail | [5, 4] | predictive knee flexion safety margin |
| clearance77p5_damped_basin04 | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_basin05 | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_basin06 | 19.046 | fail | [2, 2] | predictive knee flexion safety margin |
| clearance77p5_damped_basin07 | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_fresh_repeat | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_physics1ms | 28.000 | pass | [7, 6] | completed |
| clearance77p5_damped_reference | 28.000 | pass | [7, 6] | completed |
| clearance80 | 28.000 | pass | [7, 6] | completed |
| failed01_full10_tau10 | 18.000 | pass | [4, 7] | completed |
| failed01_full5_tau10 | 18.000 | pass | [6, 7] | completed |
| failed01_full5_tau25 | 18.000 | pass | [6, 7] | completed |
| forcecap20_70 | 10.532 | fail | [0, 0] | pelvis height failure |
| forcecap20_80 | 21.426 | fail | [0, 0] | pelvis height failure |
| forcecap40_70 | 14.126 | fail | [0, 0] | pelvis height failure |
| nominal_full10_tau10 | 18.000 | pass | [4, 6] | completed |
| nominal_full5_tau10 | 18.000 | pass | [6, 7] | completed |
| nominal_full5_tau25 | 18.000 | pass | [3, 7] | completed |
| old_failure01_smooth80 | 5.182 | fail | [0, 0] | pelvis height failure |
| smooth75 | 5.086 | fail | [0, 0] | pelvis height failure |
| smooth77p5 | 28.000 | fail | [0, 6] | completed |
| smooth80 | 5.176 | fail | [0, 0] | pelvis height failure |
| wean60 | 12.214 | fail | [0, 0] | predictive knee flexion safety margin |
| wean60_flex35 | 12.514 | fail | [0, 0] | predictive knee flexion safety margin |
| wean70 | 14.056 | fail | [0, 0] | predictive knee flexion safety margin |
| wean70_flex30 | 14.076 | fail | [0, 0] | predictive knee flexion safety margin |
| wean70_flex35 | 14.436 | fail | [0, 0] | predictive knee flexion safety margin |
| wean70_lower40 | 14.056 | fail | [0, 0] | predictive knee flexion safety margin |
| wean75 | 19.924 | fail | [0, 2] | predictive knee flexion safety margin |
| wean75_brake4 | 22.152 | fail | [0, 2] | predictive knee flexion safety margin |
| wean75_flex30 | 20.338 | fail | [1, 2] | predictive knee flexion safety margin |
| wean75_lower20 | 19.892 | fail | [0, 2] | predictive knee flexion safety margin |
