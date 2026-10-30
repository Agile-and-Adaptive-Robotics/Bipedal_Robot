# Bounded simulation jobs for ZCode and other models

Session 1. Machine: **easteregg2**. Repository: **D:\GitHub\Bipedal_Robot**.
Python: **D:\Anaconda\envs\myo\python.exe**.

Copy one complete prompt per worker. These jobs produce reviewable results
without giving every worker ownership of the controller. Run independent jobs
in separate output folders. Leave controller design, changes to acceptance
criteria, and interpreting conflicting literature to a coordinating session.

| Job | Suitable bounded work | Dependency |
| --- | --- | --- |
| A | Small configuration sweep | Existing fitted weights |
| B | Force/current analysis and plots | Saved recordings; no simulation |
| C | Reproducibility and basin checks | One configuration passing stages 1 and 2 |
| D | Recorded-state animation | Existing stage-2 recording |
| E | Backend parity and timing | Existing family audit files |
| F | Simulink inventory and isolated runs | Local MATLAB and model dependencies |
| G | Prototype a corrected curriculum adapter | Coordinator review before use |
| H | Prepare one isolated Hunt training job | Source/data audit before execution |
| I | Ground-stage protocol and bounded validation | Same candidate passes stages 1 and 2 |

Use A, B, D, and E as the first jobs for Flash or another inexpensive model:
their inputs and checks are explicit. C needs careful bookkeeping. F, G,
and H need more code review because they cross program or training interfaces.

## A. Frontal-clearance sweep with frozen neural weights

```text
Machine: easteregg2. First run hostname and stop if it differs.
Repository: D:\GitHub\Bipedal_Robot
Python: D:\Anaconda\envs\myo\python.exe
Read D:\GitHub\Bipedal_Robot\AGENTS.md, following only easteregg2's machine
bullet, and the Python, MuJoCo, and SNS skills in
C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills.

Perform a bounded configuration sweep. Do not redesign the circuit or edit
production files, trained weights, manuscript files, or CHATGPT_REPORT.md.
Do not interact with Ben's browser or connectome-editor tabs.

Inputs:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\TRAINING_RECORD.md
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\native_joint_pf_mono_20261030_v1\release_candidate\stage2.json

The matched abductor-scale-8 configuration passes deafferented air but fails
the afferented frontal-clearance check. The scale-1 configuration passes
afferented air. Keep the weights directory and every other supplied setting
fixed. Make separate stage-2 configurations with axis_feedback 0.25, 1.5,
and 4.0, and hip_abduction_goal 0.08. Check the training record first: reuse
an already completed identical test rather than duplicating it.

Write all new configurations, logs, and results below
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\frontal_sweep_20261030\
Use a unique suffix if that directory already contains another worker's work.
Run each configuration in a fresh process with the existing evaluator:
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py' '<absolute configuration path>'
The evaluator's default sparse backend is already accelerated. Do not increase
the 2 ms physics/neural timestep. Use at most 3 simultaneous trials.

Report every result, including failures: knee cycles, phase error, period and
CV, minimum ankle/knee widths, crossing fraction, effort proxy, and assistance.
Do not call air stepping walking. Do not alter acceptance thresholds.
If one configuration passes stage 2, rerun that same complete configuration
at stage 1 in a fresh process. The evaluator disables sensory pathways there.
Do not advance to stage 3 or 4. Return the exact configurations, source and
weight hashes, and a short result table. Write only a job-local report.
```

## B. Diagnose the command-to-force chain from existing data

```text
Machine: easteregg2. Repository: D:\GitHub\Bipedal_Robot.
Interpreter: D:\Anaconda\envs\myo\python.exe. Verify hostname first.
Read AGENTS.md and the machine's Python/MuJoCo/SNS skills. Do not install
packages, edit neural weights, rerun training, or touch the live editor.

Compare these matched air-stage recordings:
OLD: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\swing_release_20261030_v1\air.mechanics.npz
NEW: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\native_joint_pf_mono_20261030_v1\air_retry\stage2.mechanics.npz
Each has a neighboring .result.json describing the actual configuration.
Read swing_diagnostics.py and analyze_swing_release.py from
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training.
The equality-aware mechanics implementation is:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\functional_gait2392\mechanics.py

Use t >= 5 and t < 13 s. Resolve every channel by saved names. Define RG
swing from RG_F > RG_E, and distinguish that neural phase from actual foot
unloading. Compare PF voltages, MN voltages, excitation commands, native
activation states, PF/reflex currents, and active/passive knee torques for
rectus femoris, the long hamstrings, and biceps femoris short head.
Use the saved equality-aware rho matrix. Do not substitute raw actuator
moment arms. Do not call all positive swing torque erroneous: it can brake
the leg. Relate the torque timing to knee velocity and actual motion.

Write plots and a report under
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\force_chain_20261030\
Use Arial, readable labels, and distinguish left/right by line style as well
as color. Retain numeric sources and units. State which recordings are
downsampled diagnostics. No new walking or neural-pass claim is authorized
from current traces alone. Return a short causal comparison and limitations.
```

## C. Fresh-process reproduction and ±1% basin checks

```text
Machine: easteregg2. Repository: D:\GitHub\Bipedal_Robot.
Python: D:\Anaconda\envs\myo\python.exe. Check hostname first.
Read AGENTS.md and
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\TRAINING_RECORD.md.

Dependency gate: identify ONE explicitly recorded complete configuration
that passes both early physical stages. Separate passing configurations for
the two stages do not meet this dependency. If none exists, report that fact
and stop; do not invent a winner or change thresholds.

Freeze the winning source/configuration and copy its weight files into a
new job folder below D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\basin_20261030.
Keep the originals unchanged. First reproduce both nominal stages in fresh
Python processes through the existing hunt_training\physical_eval.py.
Then run 8 deterministic perturbation seeds, independently multiplying every
fitted synaptic edge and sensory gain by 0.99 or 1.01. Perturb the explicitly
tuned scalar configuration gains too. Keep booleans, stage numbers, timestep,
measured metadata, neuron labels, and topology unchanged. Document exactly
which fields were perturbed, their count, and seeds. Do not perturb metrics
inside a weight file as though they were parameters.

Use at most 3 simultaneous fresh-process evaluations. Record per-stage
movement, frontal-clearance, cadence, and assistance metrics for every seed.
Preserve failed trials. Never aggregate away NaNs or frozen runs. Do not
advance the curriculum or retune the winner. Return a job-local report and
the exact reproduction/perturbation artifacts, with source and weight hashes.
Do not edit CHATGPT_REPORT.md or interact with Ben's editor.
```

## D. Plot and animate an existing successful air run

```text
Machine: easteregg2. Repository: D:\GitHub\Bipedal_Robot.
Python: D:\Anaconda\envs\myo\python.exe. Verify hostname first.
Read AGENTS.md plus the Python and MuJoCo skills. Do not retrain anything,
edit the recorded trajectory, or interact with the user's editor/browser.

Source:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\native_joint_pf_mono_20261030_v1\air_retry\stage2.npz
Read the adjacent stage2.result.json before describing it. It is neural AIR
stepping with full vertical support, not supported ground walking.

Run the recorded-state renderer:
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\render_air_preview.py' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\native_joint_pf_mono_20261030_v1\air_retry\stage2.npz'

It writes preview\stage2_air_neural.gif and .png beside the recording.
Confirm the MuJoCo model has the same qpos size as qfull. Render saved full
states without reseeding, correcting the poses, or interpolating a different
motion. Inspect the resulting GIF and still. Keep anterior/oblique views,
leg identification, support status, and time labels readable. Preserve the
provenance JSON. If rendering fails, fix only a job-local renderer copy and
explain the change; do not alter the simulation or its assessment.
Add a plot of hip/knee angles and RG outputs using saved channel names. The
main runner's q channels are already degrees. Its ctrl values are excitation
commands; do not label them native muscle activation states.
```

## E. Compare speedup kernels on s3k, w2lvar, and syn6

```text
Machine: easteregg2. Repository: D:\GitHub\Bipedal_Robot.
Python: D:\Anaconda\envs\myo\python.exe. Verify hostname first.
Read AGENTS.md and the Python/SNS/MuJoCo skills. No package installation,
production edits, larger timestep, fastmath, or changed channel equations.

Read the existing reference procedure:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\verify_preallocated_families_20261028.py
and backend files under
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\
The family audits are in spinal\reports_20261028_family_audit.

Make a job-local adaptation of the reference procedure under
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\backend_parity_20261030\
Use absolute paths to the original family audits. Redirect outputs into your
job folder rather than overwriting existing audit results. Test every existing
s3k/w2lvar/syn6 human/robot audit combination. Compare native and preallocated
voltages, outputs, and NaP h gates on identical seeded input sequences before
timing. Include a warmed sparse-pairwise comparison where supported.

These kernels reject spiking and electrical-synapse builds. Do not apply a
graded kernel to a spiking model or silently replace it with graded neurons.
Warm JIT separately. Time serially with the same BLAS-thread settings and
report compile/setup time separately. Neural-only benchmarks must remain
labeled neural-only. Do not infer whole-model real-time speed from them.
Return exact errors, timings, build counts, dependencies, and hashes. If any
parity check fails, preserve it and do not declare that backend usable.
```

## F. Simulink inventory and isolated existing-model runs

```text
Machine: easteregg2. Repository: D:\GitHub\Bipedal_Robot.
MATLAB: D:\Program Files\MATLAB\R2025a\bin\matlab.exe.
First check hostname. Read AGENTS.md and the MATLAB skill at
C:\Users\Ben Bolen\Documents\GitHub\ZCode_Skills\matlab\SKILL.md.

Work from D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape.
Read README_SNS_Simscape.md and the corresponding demo builder/runner before
running it. Inventory existing graded and spiking demo SLX files, library
dependencies, solver settings, logged signals, and MATLAB release metadata.
Run headlessly with matlab -batch. Do not launch or manipulate Ben's GUI.

Use SimulationInput overrides and a new output folder below
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\simulink_20261030\
for copies and logs. Do not rebuild or overwrite SNS_Library.slx, any existing
demo, manuscript figure, or existing winner. Test existing library/unit models
and a knee-reflex graded/spiking pair where their dependencies load. Preserve
each model's documented solver and step size. Report a version/license failure
as such, rather than changing mechanics or claiming the model is absent.

The new helpers run_sns_model_extended.m and prepare_spiking_library_extended.m
in the SNS_Simscape root provide a tested private-library route. The shared
library currently lacks the saved spiking demo's required blocks; build only
the private copy in your new job directory, using a fresh MATLAB process.
For optional Accelerator/Fast Restart tests, measure actual supported runs;
keep cold compilation separate and compare the logged results. Do not assume
that an editable neural Simulink model is a closed-loop Gait2392 simulation.
The MuJoCo bridge and true Simscape body are distinct routes. Return a table
of what loaded, ran, failed, and why, with the output paths and source hashes.
Do not edit CHATGPT_REPORT.md or interact with the live browser/editor.
```

## G. Prototype a corrected five-stage adapter for another graded family

```text
Machine: easteregg2. Repository: D:\GitHub\Bipedal_Robot.
Python: D:\Anaconda\envs\myo\python.exe. Check hostname and read AGENTS.md,
spinal\DESIGN.md, and the Python/SNS/MuJoCo skills before working.

Create an opt-in prototype below
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\corrected_curriculum_adapter_20261030.
Use s3k first. Do not edit runner.py, params.py, build_network.py, default
topologies, trained winners, or current Hunt files. Do not retrain networks.

Study the CURRENT corrected evaluator and force gates:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\ground_contact_check.py
The historical _curriculum*.py scripts call their stage 4 --no-ground. That
is not Ben's corrected harness-supported GROUND-stepping stage.

Build only a wrapper/logging adapter selecting the existing graded s3k net
and the existing exact/preallocated backend. Keep MuJoCo 2.3.7,
implicitfast, and 2 ms steps. Record each ground contact once, with signed
world-frame three-axis foot forces, per-foot normal load, and the harness
force/torque estimate. Preserve full qpos and native neural names.

Prove a logging-only run leaves qfull, ctrl, and neural trajectories identical
to the uninstrumented reference. Explain frame adapters explicitly. Do not
use legacy doubled heel+toe totals to certify ground load. Do not count
frozen poses, sliding contacts, teacher-forced motor outputs, or an external
stride clock as neural stepping. Return source, parity evidence, and a
proposed run table for coordinator review. Do not run an optimization campaign
or claim any curriculum pass before that review.
```

## H. Prepare and run one isolated Hunt fit, without changing the accepted circuit

```text
Machine: easteregg2. First run hostname and stop if it differs.
Repository: D:\GitHub\Bipedal_Robot
Python: D:\Anaconda\envs\myo\python.exe
Read AGENTS.md, the python, mujoco, neuroscience, and sns-toolbox skills, and:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\TRAINING_RECORD.md
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\train_feasible_hunt.py
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\fit_native_joint_pf.py

First prepare a source/data audit in:
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\hunt_fit_plan_20261030
Resolve every actual input path, output path, environment switch, resume guard,
and fitting budget. Confirm the accepted step-A winner, the fast/slow feasible
teacher datasets, timing weights, and current fitted weights exist. Teacher data:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\swing_release_20261030_v1\feasible_targets
Current candidate:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\native_joint_pf_mono_20261030_v1

Identify a SINGLE justified fitting change from the failure diagnostics. Do not
redesign the network, combine sartorius/gracilis with RF, change the four-cell
RF/HAM topology, change acceptance gates, or refit A/B merely to improve C.
Record Hunt A/B/C/D explicitly: accepted intrinsic rhythm, accepted expected-
afferent timing, PF/MN feedforward, and residual sensory correction.

If the existing driver cannot direct ALL outputs to a new campaign without
overwriting input weights, stop after the audit and explain the required adapter.
Otherwise create a new campaign named zcode_hunt_fit_20261030_v1 under the
existing ignored hunt_training\runs directory. Never reuse an existing campaign.
Keep accepted inputs read-only, use the existing compiled fitting kernels with
one BLAS thread, and cap this job at one fit/30 minutes. Preserve incomplete
output as incomplete if it times out. Do not mark it TRAINING_COMPLETE by hand.

Return actual targets, native dynamic PF waveform errors, active-channel motor
R²/range, coupled rhythm, and physical stage-1/stage-2 results for the SAME
configuration. A low aggregate command RMSE is not a successful fit. Do not run
stage 3 or later. Do not modify CHATGPT_REPORT.md, source winners, the browser
editor, production source, or dissertation files. Write a separate RESULT.md
and machine-readable manifest in your private job folder.
```

## I. Prepare corrected ground validation, gated on the same early-stage winner

```text
Machine: easteregg2. First run hostname and stop if it differs.
Repository: D:\GitHub\Bipedal_Robot
Python: D:\Anaconda\envs\myo\python.exe
Read AGENTS.md and these current sources:
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\TRAINING_RECORD.md
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\ground_contact_check.py
D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\balance.py

Output only under:
D:\GitHub\Bipedal_Robot\tmp\zcode_jobs\ground_protocol_20261030
First identify a single candidate that passes BOTH corrected stages 1 and 2 in
fresh processes with identical weights/gains/support policy. Separate stage
winners do not qualify. The current gain-8 candidate has not met this condition.
If none qualifies, do NOT run later stages. Produce explicit consumed-key JSON
recipes and a dependency table for review, with the unmet gate stated clearly.

If one qualifies and its evidence is supplied, first run only stage-3 assisted
standing using the existing evaluator and accelerated backend. Log actual foot
forces in three world coordinates, actual harness force/torque estimates,
posture, commands, native activations, and failed terms. Use a small bounded
assistance comparison; do not add pose tracking, inverse control, or a stride
clock to obtain a pass. Keep stepping/load/crossing thresholds unchanged.

Stage 4 must include a ground plane and repeated alternating stance/swing.
Prepare a conservative harness ladder: initial sagittal restraint and high
support, then separate lowering/restraint/support reductions. Define actual
force and unloading checks at each rung. Do not execute stage 4 until the same
candidate passes standing and coordinator review. Stage 5 follows a passed
ground solution; ablate only pathways that are actually active.

Do not run historical _curriculum*.py stage 4, which removes ground. Do not
call a squared-command cost metabolic power. Do not change shared training
records, CHATGPT_REPORT.md, model files, weights, the editor, or dissertation.
Return RESULT.md, configurations, exact commands, and a compact results table.
```
