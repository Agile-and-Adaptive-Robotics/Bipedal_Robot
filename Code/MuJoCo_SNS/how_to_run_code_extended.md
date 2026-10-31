# Extended walker simulation guide

Session 1, easteregg2, 2026-10-30. Companion to `HowToRunCode.md`; the older
guide is preserved. Commands below use PowerShell unless a block says MATLAB.
They describe implemented entry points, including their limitations. A successful
command, an oscillating neural trace, and a passing walking assessment are three
different results.

Read this guide by **model, program, neuron type, and then operation**. The
standalone worker prompts are in `ZCODE_SIMULATION_PROMPTS.md`, beside this file.
Run artifacts belong in an ignored directory, not beside manuscript figures.
No command here authorizes changing dissertation results, committing, or pushing.

## 0. Select the machine and establish one output directory

Start with `hostname`. The executable paths in this section are for
**easteregg2**, whose repository is `D:\GitHub\Bipedal_Robot`.

```powershell
$machine = hostname
if ($machine -ne 'easteregg2') { throw 'Read this machine''s AGENTS.md entry before selecting executables.' }
$repo = 'D:\GitHub\Bipedal_Robot'
$py = 'D:\Anaconda\envs\myo\python.exe'
$matlab = 'D:\Program Files\MATLAB\R2025a\bin\matlab.exe'
$mj = Join-Path $repo 'Code\MuJoCo_SNS'
$spinal = Join-Path $mj 'spinal'
$hunt = Join-Path $spinal 'hunt_training'
$sns = Join-Path $repo 'Code\Matlab\SNS_Simscape'
$run = Join-Path $repo ('tmp\walker_runs\' + (Get-Date -Format 'yyyyMMdd_HHmmss'))
New-Item -ItemType Directory -Path $run -Force | Out-Null
$env:CONDA_PREFIX = 'D:\Anaconda\envs\myo'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:NUMEXPR_NUM_THREADS = '1'
Set-Location $repo
```

The confirmed environment is Python 3.10, MuJoCo 2.3.7, SNS-Toolbox 1.5.2,
and MATLAB R2025a Update 1. Do not upgrade numerical packages during a campaign.
Use the explicit Python executable; bare `python` can invoke the Store stub.
For Python work spanning several statements, save a `.py` file and execute it.
Multiline `python -c` quoting is unreliable on these machines.

On another machine, read only its matching AGENTS.md machine entry and substitute
its repository and environment paths. The Simulink–MuJoCo bridge was demonstrated
on **DESKTOP-5Q16KE9**, the laptop; that does not establish an installed bridge on
easteregg2. No compiled `.mexw64` bridge was found in this checkout.

Use one folder per worker and candidate. Keep logs, JSON configurations, NPZ
recordings, MAT results, and generated cache files under `tmp\walker_runs` or the
existing ignored `hunt_training\runs` tree. Preserve selected evidence with its
source hash and parameters before cleaning scratch files. Do not delete existing
runs or change ignore rules broadly to hide source files.

Ignored weights, converted meshes, and datasets do not travel with a Git pull.
Before moving a job to another machine, verify every input path there and retain
the selected data/weights with their manifest. A manifest is not a data backup.

### Model and capability map

| Model or controller | Body/program | Graded neurons | Spiking neurons | Current acceleration route |
| --- | --- | --- | --- | --- |
| s3k | Converted human Gait2392, Python MuJoCo | Implemented | Separate stock hybrid mirror | New launcher: native, preallocated, or sparse for graded |
| w2lvar | Same converted Gait2392, Python MuJoCo | Implemented | No verified w2lvar spiking selector | Same graded launcher |
| syn6 | Same converted Gait2392, Python MuJoCo | Implemented | No verified syn6 spiking selector | Same graded launcher |
| biw2lvar-Hunt-PT, four-cell RF/HAM | Same human model, Hunt evaluator | Implemented, experimental training | No faithful port implemented | Existing fixed-tau sparse/Numba evaluator and native reference |
| W2L original/split/contact variants | Virtual-walker geometry, `w2l_mujoco` / `w2l_cpg` | Several distinct transcriptions | Split-RG hybrid version exists | Supported graded compiled kernel, verified per variant |
| W2L split-RG with BPAs | Virtual walker, twelve BPA actuators | Implemented | Implemented hybrid split-RG | New single-trial BPA launcher; graded preallocation or native |
| Li-derived walker | Virtual walker, Python MuJoCo | Graded adaptation | Do not infer a faithful spiking port from the paper | Existing diagnostic scripts; no new verified fast adapter |
| Di Russo adaptation | Gait2392, Python MuJoCo | Experimental adaptation | No verified port | Existing `dirusso_fast.py` API and parity script |
| AH001 CAD robot | Robot torso runner, 54 routes | Hybrid SNS plus engineering allocator | No implemented spiking replacement | Existing `exact_fast` backend |
| SNS knee/cup/reflex, Deng CPG, planar air-step demos | MATLAB/Simulink | Implemented | Selected demos only | Normal mode; test Accelerator/Fast Restart separately |

Names are not interchangeable. The W2L contact-driven network, the split-RG BPA
network, the human w2lvar network, and the AH001 controller are different builds.
For example, the human production s3k/w2lvar/syn6 censuses are 410/888/794 neurons;
the split-RG BPA transcription is 87 neurons, and the contact-driven W2L build is
95. Count the actual compiled candidate when options add conditional topology.

## 1. Human Gait2392: s3k, w2lvar, and syn6 in Python/MuJoCo

### 1.1 Model, coordinates, and neuron selection

The default body file is:

`D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\mjc\gait2392_simbody\gait2392_simbody_cvt3.xml`

This legacy converted human model uses a **z-up world**. Its anatomical adapter
maps OpenSim `(x anterior, y superior, z right)` to world `(x, -z, y)`.
Knee flexion is negative. Use saved coordinate definitions when interpreting hip
angles and moment arms; naming alone does not establish a torque sign. AH001's
y-up model is described separately in section 4. Do not rotate old geometry or
change gravity to make different baselines look uniform.

The new launcher is `Code\MuJoCo_SNS\run_spinal_extended.py`. It selects the
existing network without editing `params.py`, the builders, or `runner.py`.

```powershell
& $py "$mj\run_spinal_extended.py" --family s3k --backend preallocated --output "$run\s3k.npz" -- --eval
& $py "$mj\run_spinal_extended.py" --family w2lvar --backend preallocated --output "$run\w2lvar.npz" -- --eval
& $py "$mj\run_spinal_extended.py" --family syn6 --backend preallocated --output "$run\syn6.npz" -- --eval
```

These are baseline evaluations, not claims that old winners walk under current
physics. Each output must have a new name. A neighboring `.launch.json` records
the family, backend, model, flags, source hashes, build census, actuator names,
and elapsed time. Preserve both files.

`--family spiking --backend native` selects the **stock hybrid spiking mirror**:
spiking afferent/RG/PF/interneuron layers, analog motor neurons, and filtered
readouts. It is not a spiking version of the current Hunt controller.

```powershell
& $py "$mj\run_spinal_extended.py" --family spiking --backend native --output "$run\spiking_stand.npz" -- --stand-eval 0.04
```

The mirror advances its neurons in 0.5 ms inner steps, four per 2 ms physics
step. Raw LIF membrane voltages, spike events, filtered readouts, and graded MN
commands are different signals. Check `neuro_names` and the builder's readout
mapping before labeling a plot. A filtered RG readout is not a spike raster.
The graded fast kernels reject spiking/electrical networks; selecting spiking
with `preallocated` is intentionally an error.

### 1.2 Acceleration and a reproducible comparison

Use `--backend native` for the SNS reference, `preallocated` for the fixed-tau
graded evaluator, or `sparse` for the pairwise-reduction evaluator. The launcher
copies the currently applied conductances into its fast representation before
each forward call. New topology created during a run needs its own adapter.

```powershell
& $py "$mj\run_spinal_extended.py" --family s3k --backend native --output "$run\native.npz" -- --eval
& $py "$mj\run_spinal_extended.py" --family s3k --backend preallocated --output "$run\fast.npz" -- --eval
```

Compare `t`, `qfull`, `ctrl`, and `neuro` from these fresh processes, together with
the manifests. The launcher's preallocated route was checked for all three graded
families on complete 15.3 s trajectories: all three arrays matched exactly.
This establishes those configurations' numerical parity, not every future edit.
The three-family verification is `tmp\test_extended_launchers.py`; its raw results
are disposable, not a new production acceptance contract.

Benchmark serially on an otherwise idle machine. Time these separately:

1. Python startup/model loading and first Numba compilation.
2. Warm neural updates on identical saved inputs.
3. Full closed-loop physics with rendering disabled.
4. Optional plotting and animation export.

Report simulated seconds, wall seconds, `simulated/wall` real-time factor, and
`native wall/fast wall` speedup. A neural-only ratio must not be reported as the
whole simulation speedup. Concurrent verification runs are not a controlled
performance benchmark. Save model/source hashes, thread settings, and parity error.
Leave MuJoCo at 2.3.7, `implicitfast`, and the established 2 ms human-model step.

`w2l_cpg\preallocated_exact_eval.py` and
`w2l_cpg\sparse_pairwise_eval.py` contain the supported kernels. They preserve
the fixed sodium-channel inactivation time constant and prescribed reduction
order. The older `w2l_cpg\fast_eval.py` sparse evaluator has failed full-trajectory
parity; it is not an approved interchangeable speedup. A larger physics step,
`fastmath`, altered summation order, adaptive tolerances, contact damping, or a
different sodium-channel time constant changes the experiment.

To extend acceleration to another graded model, create a local wrapper, retain
the native path, assert the neuron/synapse types and fixed-tau semantics, then
compare saved random-input neural traces and full closed-loop trajectories.
Reapply dynamic gains and reset all voltage/gating state. Do not patch installed
SNS-Toolbox globally. The Di Russo-specific API is not a generic SNS evaluator.

The core operation inside an isolated wrapper, after building a supported
graded `net`, is:

```python
from preallocated_exact_eval import PreallocatedExactForward
raw = net.compiled
fast = PreallocatedExactForward(raw)
# Before each update, if gains can change:
fast.gmax[:] = raw.g_max_non[fast.post, fast.pre]
fast.forward(input_vector)
# For a saved trace, assign into preallocated storage or copy:
voltage_trace[k] = fast.V
```

This fragment is not a complete walker loop: use the constructed launcher rather
than omitting sensor updates, state resets, muscle updates, or output scaling.
Never append the same reused voltage buffer to a list without copying it.
The Hunt controller has its own gain synchronization; preserve that implementation.
Spiking acceleration would additionally need event, threshold/refractory,
synaptic decay, and readout-state parity. No such generic fast replacement is
constructed by this fragment.

### 1.3 Time, fitted weights, and older command-line flags

`--eval` runs the legacy stand–walk–stand schedule. **`--time N` changes the final
standing interval, not the total simulation duration**: `--eval --time 0.8`
still records 15.3 s. For a genuinely short standing smoke test, use
`--stand-eval 0.8`. Do not shorten an evaluation below its gate's scoring windows.

Pass existing runner flags after the separator `--`. Examples:

```powershell
& $py "$mj\run_spinal_extended.py" --family s3k --backend preallocated --output "$run\push.npz" -- --eval --push 20 --push-time 7 --push-axis y
& $py "$mj\run_spinal_extended.py" --family s3k --backend preallocated --output "$run\air_noaff.npz" -- --eval --no-ground --no-afferents
```

Here world `y` is lateral, not anatomical superior. Inspect the pulse duration in
`runner.py` before reporting an impulse. `--harnessK`, `--rig-scale`,
`--leg-damping`, `--stance-width`, `--no-harness`, and `--trunk-load` alter support
or mechanics; record them. A script that removes ground is an air experiment.

`--fitted` reads the legacy fitted file; `--best`, `--best5` through selected
versioned flags read older saved winners. There is **no literal `--bestN` flag**.
Do not use the old v10 recipe as a current walking result. Also inspect
`connectome_gains.json`, which can be applied automatically by the stock builder.
Record every consumed parameter file, not only the command line.

### 1.4 Edit and display this network

Before touching production `runner.py`, `build_network.py`, or `params.py`, read
`spinal\DESIGN.md`. Keep the production census and documented anchor evaluation
unchanged. Family builders include `build_network_w2lvar.py`,
`build_network_syn6.py`, and `build_network_spiking.py`; verify actual filenames
with `rg --files` before making a worker edit. Put experimental wrappers and
parameter sets in an isolated campaign. New gain knobs must be saved and loaded.

`runner.py --scope` uses the live neural scope; `--view` opens the MuJoCo viewer,
and `--realtime` slows display to wall time. Use them for inspection, not speed
measurements. For instance, append `-- --eval --scope --view --realtime` to a new
launcher invocation. Closing a viewer can interrupt a run; retain a headless
evaluation as the quantitative record.

`spinal\draw_circuit.py` supports `--which core|full|deng|weights|all`,
`--source params|fitted|best`, `--fmt pdf,svg,png`, and `--vclasses`. Its legacy
source selectors do not automatically draw the current Hunt candidate. It writes
to the shared `spinal\figures` directory; use a copied/staged version if those
outputs are being edited. For a faithful candidate diagram, export the actual
compiled graph, including enabled edges and the loaded weights.

### 1.5 Choose plotted signals

The new helper `plot_recording_extended.py` lists signals before plotting them:

```powershell
& $py "$mj\plot_recording_extended.py" "$run\s3k.npz" --kind joints --list
& $py "$mj\plot_recording_extended.py" "$run\s3k.npz" --kind neurons --list
& $py "$mj\plot_recording_extended.py" "$run\s3k.npz" --kind joints --signals 'knee*' 'hip_flexion*' --start 2 --end 13 --out "$run\joints.png"
& $py "$mj\plot_recording_extended.py" "$run\s3k.npz" --kind neurons --signals 'RG_*' --start 5 --end 9 --out "$run\rg.png"
```

`--signals` uses exact names or shell-style wildcards; PowerShell single quotes
keep the wildcard for Python. `--indices 0 1 2` selects columns by index.
`--kind controls` plots named command channels. In a mechanics diagnostic,
`--kind activation` selects native activation state. The stock recording's
`act` array is a named subset of **excitation commands**, despite its old name;
it must not be relabeled physiological activation. Full `ctrl` contains all
actuator commands, but their names live in the launch manifest rather than the
legacy NPZ. Use the manifest or an explicit export for all-actuator analysis.

Main recordings contain `t`, joint `q` in degrees, full `qfull` in native units,
neural `neuro` and `neuro_names`, `ctrl`, selected command `act`, and other logs.
Use `qfull` for replay. The legacy heel-plus-toe contact total can double-count
a contact; it cannot certify the corrected ground-load gate.

## 2. Human Gait2392: biw2lvar-Hunt-PT and four-cell RF/HAM

### 2.1 What this controller contains

This is an experimental graded SNS, not a renamed stock w2lvar winner. It has
bilateral RGs, joint PFs, RF/long-hamstring PFs, motor pools, sensory pathways,
causal phase filters, and the opt-in task-balance pathway. Each RF/HAM module is
**two excitatory PF cells and two inhibitory relays**. The RF PF receives hip
flexor/knee extensor drive; the HAM PF receives hip extensor/knee flexor drive.
Each excites its own relay, which inhibits its opponent. RF drives rectus femoris;
HAM drives semimembranosus, semitendinosus, and biceps femoris long head.
Biceps femoris short head remains a knee-flexor channel. Sartorius and gracilis
remain separate. `HUNT_BANK_*` labels mean engineering causal phase filters,
not an anatomical neuron class described by Hunt.

The present controller has no tested spiking mirror and no faithful Simulink
export. A picture of its topology is not an executable parameter-equivalent port.
Use its native graded reference or its supported graded accelerated evaluator.

### 2.2 Current training checkpoint and the immediate task

Read `spinal\hunt_training\TRAINING_RECORD.md` for the research record. It is
separate from the ZCode handoff `CHATGPT_REPORT.md`. Current campaign:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\runs\tonic_abduction_20261030_s1_v1\0p02`

A 0.02 conductance increment on existing tonic-PF outputs to hip-abductor MNs
passes **both air stages with the same saved parameters**, including fresh-process
repeats. Stage 1 has 3/3 knee cycles, 2.467 s period, and minimum ankle-origin
width 0.0991 m; stage 2 has 7/7 cycles, 1.214 s period, and width 0.1086 m.
This increment is not a 2% muscle-activation floor. Six scalar sensitivity
checks pass, but only one of two simultaneous synaptic-jitter checks passes.

The same weights pass assisted stage-3 standing with about 55% of body weight
on the feet and 45% carried by the vertical spring. The harness, joint
restraints, and inherited posture/reflex pathways remain active. Removing
the added Peterka/Ting-inspired cells also passes; their benefit is unproven.
**No candidate passes the corrected alternating ground-contact gate.**
See the training record for failed ground-transfer and symmetric-PF comparisons.
ZCode A's scale-8 frontal sweep fails all three tested cases; its job report
and ZCode B's force-chain diagnostics are under `tmp\zcode_jobs\`.

### 2.3 Hunt's four training operations

These operations train the circuit; they are not four alternative names for the
five physical validation stages:

| Operation | What is fitted | What remains fixed | Evidence required |
| --- | --- | --- | --- |
| A | Intrinsic RG timing, free period near twice intended afferented period | Body-independent training definition | Native free oscillation, period, bilateral phase |
| B | Expected-afferent timing gains | Accepted intrinsic RG parameters | Full coupled network recovers intended cadence with defined inputs |
| C | PF/MN feedforward command patterns | Accepted A/B | Dynamic waveform fit and physical motion, not only aggregate RMSE |
| D | Residual sensory correction | Accepted A/B/C | Closed-loop improvement, with causal sensory ablations |

The approximately 80/20 feedforward/feedback target is an operational training
choice. Nonlinear muscle force is not an additive 80/20 force decomposition.
Torque-teacher trajectories are external-control diagnostics, never neural
curriculum passes. Distinguish desired excitation from MuJoCo activation state
when fitting forces and replaying targets.

`train_feasible_hunt.py` is an experimental campaign driver. Its inputs currently
include simulation-derived teacher data, a shared accepted step-A file, and
campaign-specific timing/weight files; it is not a general one-command retrain
from arbitrary IK. Supported environment switches include:

| Variable | Meaning |
| --- | --- |
| `HUNT_FEASIBLE_CAMPAIGN` | New campaign directory name under `runs` |
| `HUNT_TIMING_SOURCE` | Existing accepted timing source to copy |
| `HUNT_FEASIBLE_TARGET=mono_command` | Monoarticular command targets rather than net-torque targets |
| `HUNT_FEASIBLE_MODES=afferented` | Fit the specified available cadence set |
| `HUNT_PF_INACTIVE` | Inactive PF target for the explicit experimental target variant |
| `HUNT_JOINT_PF_INHIBITION` | Opt-in construction inhibition scale, subject to loaded learned edges |
| `HUNT_FEASIBLE_PF_ONLY=1` | Diagnostic PF fitting path; inspect resume behavior before using |

Existing weight files can skip fitting. The completion guard and resume behavior
must be checked before calling a campaign complete. A PF-only resume with existing
weights is not equivalent to a clean PF-only fit. Copy the relevant driver into
a private campaign if modifying its assumptions. Preserve A/B and C/D provenance.

`fit_native_joint_pf.py` fits actual free joint-PF dynamics using compiled kernels.
Its `HUNT_NATIVE_PF_CAMPAIGN`, `HUNT_NATIVE_SLOW_WEIGHT`, and
`HUNT_NATIVE_TARGET=mono_command` switches select explicit experiments. It still
uses campaign-specific input datasets. Run the worker's source audit first.
Low whole-vector RMSE can hide inactive muscle channels, wrong cadence, and a
latched PF; inspect active-channel R², amplitude, antagonist inhibition, and
physical unloading. The separate record retains failed fits.

### 2.4 Corrected five-stage physical curriculum

| Stage | Physical experiment | Progress gate |
| --- | --- | --- |
| 1 | Deafferented, fully supported air stepping | Free rhythm, alternating movement, no crossing |
| 2 | Afferented, fully supported air stepping | Intended cadence, bilateral movement, phase, no crossing |
| 3 | Assisted ground standing with balance pathway | Posture, meaningful foot load, measured assistance, perturbation recovery |
| 4 | Harness-supported **ground stepping** | Repeated alternating stance/swing, real contacts, unloading, measured support |
| 5 | Accepted ground solution under perturbations/reduced dependence | Same ground gates plus the stated robustness or assistance reduction |

The older `_curriculum*.py` files explicitly remove ground in their stage 4.
Their numbering is historical and does not implement this table. Never describe
air stepping as harness-supported ground walking.

The corrected evaluator is `hunt_training\physical_eval.py`. It accepts one
JSON configuration path. Make a **new JSON file**, copied from the candidate
being evaluated, with an explicit absolute `weights_directory`, `backend`, and
`stage`. The JSON file's directory and stem determine the output prefix; a JSON
key named `output` does not redirect this evaluator. Then run:

```powershell
& $py "$hunt\physical_eval.py" "$run\stage2_candidate.json" *> "$run\stage2_candidate.log"
```

Do not construct an incomplete configuration from memory. Use the current
candidate's saved stage JSON, inspect the evaluator's key reads, and preserve all
timing, drive, gains, options, support, and axis settings. Unknown keys can be
silently unused; prove a perturbation or gain is consumed. Inspect a result's
`assessment.accepted`, failed criteria, and trace provenance, not only its score.

For the existing baseline early-stage retry only:

```powershell
& $py "$hunt\retry_air_stages.py" --campaign "$hunt\runs\native_joint_pf_mono_20261030_v1"
```

This writes fixed retry names inside that campaign. It is not an isolated sweep
and can replace previous retry evidence. For new work, copy its JSON recipes and
use unique output prefixes. `verify_release_candidate.py` is specifically the
gain-8 comparison, not a general promotion command. Prompt A uses new outputs.

The current ground gate includes at least three force onsets per foot, meaningful
load, stance/swing duration, duty range, bilateral alternation/coverage, tilt, and
2 cm frontal clearance. It counts each contact once and uses force hysteresis
and filtering. See `ground_contact_check.py` and `physical_eval.py` for thresholds.
A visible knee cycle with the foot continuously touching or skating is not a
completed ground step. Keep acceptance thresholds fixed while tuning.

### 2.5 Harness, balance, perturbations, and costs

The evaluator exposes restriction, height offset, force-controlled support,
support fraction, vertical/lateral stiffness, and lowering/travel settings.
Stage 4 may start with sagittal restraint and substantial support, then lower
toward contact, loosen lateral restraint, and reduce support in separately logged
rungs. Change one factor per rung. A different support setting is a new test,
not a retroactive label on an earlier air recording.

Use the signed harness force/torque estimate and foot-force vectors to identify
excess lateral push-off. The current harness wrench is calculated from applied
restraint terms, not a simulated single six-axis load cell with a mounting frame.
Log commanded assistance and actual assistance separately. For new sensor work,
sum all signed foot contact forces in a common frame, plus center of pressure
when meaningful. Vertical load alone cannot measure shear or twisting effort.

The Peterka/Ting-inspired branch is an engineering task-balance pathway alongside
local muscle feedback. It is not merely a renamed myostatic reflex. Record the
signals, delays, gains, and target tasks used; do not call it a literal biological
replication. Perturb only a passed stage-3 configuration before stage-4 promotion.
For stage-5 sensory ablation, verify the targeted contact/afferent edges are
enabled. Reducing an already-disabled heel PF edge is a no-op.

A metabolic-like penalty may regularize successful solutions, but activation²,
command², pressure use, and estimated muscle metabolic power are different costs.
Until a documented energy model is implemented, call a squared-command term an
effort surrogate. Preserve hard movement/load gates so a silent, frozen solution
cannot win by expending little effort.

The saved `.assistance.npz` file exposes the recorded forces directly:

```powershell
$assist = "$hunt\runs\native_joint_pf_mono_20261030_v1\air_retry\stage2.assistance.npz"
& $py "$mj\plot_recording_extended.py" $assist --kind harness --list
& $py "$mj\plot_recording_extended.py" $assist --kind harness --start 5 --end 13 --out "$run\harness_force.png"
& $py "$mj\plot_recording_extended.py" $assist --kind harness_moments --out "$run\harness_moments.png"
& $py "$mj\plot_recording_extended.py" $assist --kind foot_forces --out "$run\foot_vectors.png"
```

This stage-2 example is in air, so a lack of ground load is expected. Apply the
same selections to a future accepted ground recording rather than treating the
air plot as ground evidence. Harness moments use N·m; forces use N.

### 2.6 Select neural diagnostics and draw the actual candidate

In an evaluation JSON, `trace_cells` is an explicit list of existing compiled
cell names. `diagnose_swing` enables the mechanics diagnostic; `export_graph` enables
the compiled graph snapshot. Consult the saved graph's labels rather than
inventing RF/HAM names. Main and diagnostic NPZ files have different schemas.
The plot helper supports both named neural formats and mechanics activation.

Graph export records the **loaded simulator**. A connectome editor drawing is a
separate display. Ben's exports are under
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Mujoco_SNS_models`.
The `Hunt_fourcell_session1_*` exports describe the older failed four-cell
snapshot, not the newest monoarticular fit. Do not silently replace a diagram
with a newer model or present it as its exact weights.

For SNS-Toolbox/Graphviz drawing, `draw_network.py INPUT_DIRECTORY OUTPUT_DIRECTORY`
reads `INPUT_DIRECTORY\standing.graph.json`. It writes an HTML/SVG subnetwork
view. Ensure the output parent exists and use a new output directory. The graph
must actually be named `standing.graph.json`; a stage-2 snapshot requires an
explicit copied input under that name, with its provenance retained.

The connectome block editor and its usage are covered in section 7. **Do not
reload, replace, or import over the user's live tabs while the user is editing.**
Editor edge changes do not automatically change controller weights.

### 2.7 Animation, stills, and GIFs

`render_air_preview.py RECORDING.npz` replays full saved coordinates of the Hunt
air experiment, writes `preview\stage2_air_neural.gif`, a still PNG, and
provenance beside the recording. It does not rerun optimization or modify poses.

```powershell
& $py "$hunt\render_air_preview.py" "$hunt\runs\native_joint_pf_mono_20261030_v1\air_retry\stage2.npz"
```

Its filenames are fixed: copy the recording into a private output folder first
if an existing preview must be preserved. Label the result **fully supported air
stepping**, with its actual controller/candidate, rather than ground locomotion.
Use `--start 5 --end 13 --fps 12 --playback-speed 1` to choose the time window,
rendering density, and playback rate. GIF delays are derived from saved timestamps;
the older fixed 50 ms delay could play a subsampled recording too fast.
For another model, reuse a replay pattern only after verifying its XML, frame,
full-state dimensions, camera, and dynamic route update requirements. Never
generate missing coordinates from decorative interpolation for quantitative
evidence. Use fewer rendered frames rather than changing physics integration.

## 3. Virtual-walker W2L and Li models in Python/MuJoCo

### 3.1 Select the actual transcription

Root: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco`.
The virtual support/WalkingPath is part of the original geometry's operation.
Retaining a harness is not equivalent to demonstrating free standing.

| Entry | Neural/body scope | How to run or inspect |
| --- | --- | --- |
| `test_w2l_air.py` | Original or split-RG transcription, physical air test | `--net=orig` or `--net=split`, `--comm=1.0`, `--dur=20` |
| `test_w2l_air_afferented.py` | Scripted heel/toe feedback and causal pulse tests in air | `--causal=none`, `--causal=heel`, or `--causal=toe` |
| `test_w2l_ground.py` | Contact-driven supported ground diagnostic | Read its `KV` options and per-phase output logic before running |
| `test_li_stepping.py` | Graded Li adaptation, air hold then lowering | `--dur=20`; `--hold=3`, `--lower=1` are explicit timing options |
| `test_li_stepping_zdrop.py` | Separate Li lowering experiment | Read its header; not a neuron-mode switch |
| `test_bpa_stepping.py` | W2L with Festo BPA actuation | Legacy `--key=value` options; new launcher below avoids an Optuna study |

Examples of the **legacy script syntax**, for an owner-controlled output session:

```powershell
& $py "$spinal\w2l_mujoco\test_w2l_air.py" --net=split --comm=1.0 --dur=20
& $py "$spinal\w2l_mujoco\test_w2l_air_afferented.py" --dur=20 --causal=heel
& $py "$spinal\w2l_mujoco\test_li_stepping.py" --dur=20 --hold=3 --lower=1
```

These older scripts can write fixed artifacts or print only metrics. Read their
save calls before running, preserve existing outputs, and capture console logs.
They do not all implement argparse: `--help` may run the experiment instead of
displaying help, and an unsupported option can be ignored. Do not give them
flags copied from the human runner. A source/header check is the safe first step.

`build_w2l_orig_net.py`, `build_w2l_split_net.py`, and the contact-driven
`w2l_cpg\build_w2l_net.py` transcribe different template selections. The split-RG
build removes crossed RG-to-contralateral-PF drive, makes a right RG, and couples
the RGs through its documented commissural circuit. Match the source and census
to the claimed model. The Li Python version is an adaptation of the literature
circuit, not the original spiking calcium model.

### 3.2 BPA model: one reproducible graded or spiking run

The new `run_w2l_bpa_extended.py` calls the existing `run_bpa` protocol directly,
without starting Optuna or changing a database. Its default seed has two 20 mm
BPAs per route, 0.7 nominal retained support, drive E/F of 3/4, and 0.25 s RG
timing parameter. The protocol holds the model in air, lowers it, and then scores
retained-support behavior. Use at least 12 s for a smoke test; a 0.1 s run has
no valid scoring windows. Default is 20 s.

```powershell
& $py "$mj\run_w2l_bpa_extended.py" --neurons graded --backend native --duration 20 --output "$run\bpa_native.npz"
& $py "$mj\run_w2l_bpa_extended.py" --neurons graded --backend preallocated --duration 20 --output "$run\bpa_fast.npz"
& $py "$mj\run_w2l_bpa_extended.py" --neurons spiking --backend native --duration 20 --output "$run\bpa_spiking.npz"
```

Optional flags include `--support 0.7`, `--drive-ext 3`, `--drive-flex 4`,
`--tau 0.25`, `--comm 1`, `--diameter 20`, and `--count 2`. Each run writes its
manifest and the existing recording schema: time `t`, six joint angles `q`,
foot loads `f`, height/progression `h`, pressure `p`, and left RG-extensor `rgE`.
This is not the human full-`qpos` schema, so it needs its own replay adapter.
Pressure is in kPa. The seed uses a 620 kPa supply scale and 0.6 cap, giving a
372 kPa ceiling; it is not the ideal-tension AH001 actuator model.

Here the body step is **1 ms**, with neural calls at **2 ms**; the spiking network
performs its own smaller internal steps. Do not impose the human 2 ms physics
rule on this established virtual-walker protocol. Preserve each model's clock.
Native-versus-fast graded comparisons must include pressure and joint trajectories,
not only the RG period. These old protocol metrics are not the corrected Hunt
ground-stepping gate. A finite run with only one loaded foot is not walking.
The legacy `hip_rng_L` expression also has a sign error; use `max-min` of the
saved angle over the specified window when calculating a positive hip RoM.
Retain the raw legacy metric so the correction is explicit.

Plot the compact BPA schema with the same helper:

```powershell
& $py "$mj\plot_recording_extended.py" "$run\bpa_fast.npz" --kind joints --signals 'knee_*' --out "$run\bpa_knees.png"
& $py "$mj\plot_recording_extended.py" "$run\bpa_fast.npz" --kind pressure --list
& $py "$mj\plot_recording_extended.py" "$run\bpa_fast.npz" --kind pressure --indices 0 1 --out "$run\bpa_pressure.png"
& $py "$mj\plot_recording_extended.py" "$run\bpa_fast.npz" --kind foot_forces --out "$run\bpa_contacts.png"
```

New BPA manifests supply sorted pressure-channel names. Older recordings fall
back to explicit column numbers. The compact neural log contains only the left
RG-extensor trace; `--kind neurons` does not recover unrecorded neurons.

### 3.3 Edit, display, perturb, and train

Edit an experimental copy of the relevant template/gain mapping and preserve
per-edge provenance. `comm=0` in the split build removes commissural topology;
it is a genuine ablation. Fixed-tau sodium-channel semantics must survive a port.
Spiking adaptation parameters are not numerically interchangeable with graded
sodium-channel parameters just because both use a flag called `tau`.

Use editor templates for inspection, then compare the compiled census and enabled
edges. An editor JSON and a mined AnimatLab template are not necessarily the same
schema. Use the actual build's `voltages`, output mapping, and input-port names
when selecting signals. The BPA adapter currently saves a compact diagnostic,
not all neural voltages; add logging in an isolated copy if those are required.

The afferented-air script's heel/toe pulses test a neural sensory pathway. A body
push, platform tilt, lowered support, and a sensor pulse are different experiments.
The ground script has its own perturbation and contact-mode options; inspect
their consumed keys before use. Never call these legacy milestones the corrected
five-stage Hunt curriculum without implementing and validating the adapter.

`_curriculum_bpa.py` and `_curriculum_bpa_spiking.py` can start shared Optuna
studies. They are not a single-trial launcher. Give a worker a private database,
fixed seeds, explicit bounds, and a limited trial count before authorizing a sweep.
No supported menu converts one of these networks into current Hunt training.

## 4. AH001 CAD robot with torso: Python/MuJoCo

### 4.1 Model, neuron mode, and actuator scope

Root: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso`.
Start with its `README.md` and `FOOT_MOUNT_W2LVAR_RESULTS.md` where present.
The working runner is `run_neural_robot.py`, and `ah001_neural.xml` is the base
model. Routes are the right-side **27 per side** from `gait2392_robotbody`, mirrored
to produce **54 actuators**. The 92-route human audit file is not this robot.

AH001 uses Ben's anatomical frame **x anterior, y superior, z right**, gravity
`[0,-9.81,0]`, a ball-and-socket lumbar joint, and a 0.25 m rigid spine. The
provisional 75 kg total includes split pelvis/torso payload. Record the exact
model variant: knee limits, loop assumptions, foot mount, and self-collision
changes materially affect the result. Do not silently shorten the tibia.

Its actuator forces are presently ideal bounded reference tensions with an
activation time constant, not calibrated pneumatic hardware. Its 262-neuron,
230-input SNS supplies rhythm in combination with an engineering posture/force
allocator that provides MN target currents. This controller is distinct from the
current four-cell human controller. No spiking AH001 replacement is implemented.
Welded feet leave distal actuator channels without independent ankle/toe leverage.

### 4.2 Run and select acceleration

```powershell
$robot = Join-Path $mj 'robot_torso'
& $py "$robot\run_neural_robot.py" --model "$robot\ah001_neural.xml" --mode step --duration 20 --support 0.8 --activation-floor 0.05 --backend exact_fast --output "$run\ah001"
```

This explicit configuration is a run recipe, not a claim that the base model
passes an uncrossed ground-stepping gate. Reproduce a saved successful diagnostic
with its **exact saved flags and XML**, rather than borrowing only its support
percentage. The old 65% support recording crossed the legs; the later uncrossed
defense diagnostic used additional engineering assumptions and approximately
80% actual support. Nominal support and measured support differ.

Use `--backend reference` for comparison and `--backend exact_fast` for the
existing fast implementation. Preserve the muscle/route update loop. A generic
MuJoCo viewer does not reconstruct moving and conditional path points.
Do not install the human SNS fast wrapper around this runner.

### 4.3 Neural edits, floors, harnesses, and perturbations

Important runner flags include `--pattern direct|event`, `--unload-swing`,
`--walk-drive`, `--rg-tau-h`, `--motor-tau`, `--activation-tracking-gain`,
`--abductor-floor`, `--abduction-gain`, `--abduction-damping`,
`--abduction-force-limit`, `--stance-width`, and hip/knee control gains.
Preserve whether a change acts through MN currents, a target allocation, or an
external rig force. `--diagnostic-direct` and `--disable-muscles` are diagnostic
controls; they cannot be reported as successful muscle-driven neural walking.

The default activation floor is zero; a 0.05 floor was used in one support
preset. A minimum excitation command is not a measured 5% antagonist
coactivation index or guaranteed joint stiffness. Compare paired muscle force
and moment arms if evaluating coactivation. More abductor activation can change
hip rotation and sagittal torque as well as width.

Harness options include lateral/horizontal/vertical stiffness, angular support,
support caps, and `--support-cap-final`, `--support-wean-start`, and
`--support-wean-duration`. They expose an assistance ladder, not an already-passed
Hunt curriculum. Preserve the measured rig-force trace and foot loads at every
rung. For a new physical push protocol, add a timed wrench in an isolated wrapper
and log it; do not invent an unsupported `--push` flag.

### 4.4 Display, plots, stills, and GIFs

```powershell
& $py "$robot\replay_neural_robot.py" --run "$run\ah001.npz" --view
& $py "$robot\replay_neural_robot.py" --run "$run\ah001.npz" --gif "$run\ah001.gif" --snapshot "$run\ah001.png" --start-time 5 --time 12
& $py "$robot\plot_neural_control.py" "$run\ah001.npz"
```

The plotting script writes evidence PNG/PDF/alt text in the robot directory
using the recording stem; use a unique stem to preserve earlier evidence.
`render_defense_robot.py` accepts `--run`, `--output`, `--start`, `--end`,
`--still-time`, `--fps`, `--title`, and `--scope-label` for a labeled two-view
export. Use a caption naming actual support, controller, route/model variant,
and whether it is air stepping or ground contact. Render from the saved run;
do not remove crossing or alter poses in the renderer.

### 4.5 Simscape availability

The AH001 imported Simscape model and CAD assembly are separate assets, not an
automatic MATLAB version of the above neural/route controller. Locate
`mdl_humanoid_lower_ah001_imported.slx`, its DataFile, and geometry dependencies
before opening it. On easteregg2, probe a short runtime rather than assuming a
Simscape license or complete import. A passive CAD import may animate joints
without implementing any of the neural muscles. The bridge/laptop distinction
in section 6 applies. Do not claim functional parity from matching appearances.

## 5. Di Russo Gait2392 adaptation: Python/MuJoCo

Root: `spinal\dirusso_gait2392`. This is an experimental adaptation with failed
movement results in the project record, not a verified reproduction of the paper.

`dirusso_fast.py` is an importable module, not a configurable executable launcher.
The implemented API is:

```python
run_trial_fast(params, wsup=0.7, dur=10, W_ff=None, gains=None,
               balance_gains=None, record=False, disable_strikes=False)
evaluate_batch(jobs, n_workers=4)
```

Batch jobs use `(params, wsup, dur, W_ff, gains)` tuples. Supply the correct saved
parameter object and feedforward matrix; do not call a default matrix a new fit.
The module uses a compiled controller, vectorized sensor work, and model reuse.
`verify_fast_parity.py` compares the reference path and runs a small parallel
benchmark, but writes a fixed report in that directory. Stage a copy or archive
the existing report first. Its neural speed ratio is not a blanket physics ratio.

No implemented flag selects a spiking Di Russo controller or current Hunt PF
modules. New ports require their own neuron, connection, and force validation.
Use the same signal-selection, source-provenance, and ground-gate principles as
above, with a schema-specific plotting/replay adapter.

## 6. MATLAB/Simulink: SNS demos, walker circuits, and model coupling

### 6.1 Programs and model map

Root: `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape`.
Read `README_SNS_Simscape.md`. Add the root and the specific `demos`, `results`,
or `gui` folder needed. Do not use `genpath` over archived/development copies,
which can shadow current builders. Builders rebuild files; runners execute them.

| Model | Scope and neuron choice | Builder/runner family |
| --- | --- | --- |
| `KneeReflexDemo` | Graded muscle-reflex example | `sns_build_demo`, `sns_run_demo` |
| `KneeReflexDemo_Spiking` | Hybrid spiking reflex example | `sns_build_knee_spiking`, `sns_run_knee_spiking` |
| `BeerCupReflexDemo` / `_Spiking` | Graded or hybrid cup task | `sns_build_beer_demo` / `sns_build_beer_spiking` and matching runners |
| `BPACPGLegDemo` | Graded BPA CPG leg example | `sns_build_cpg_demo`, `sns_run_cpg_demo` |
| `SNS_SpikingRG_MN` | Spiking RG and graded MN demonstration | `sns_build_rgmn_spiking`, `sns_run_rgmn_spiking` |
| `SNS_Deng_RG`, `SNS_Deng_PF`, `SNS_Deng_CPGDemo` | Deng-style graded persistent-Na modules | `sns_build_deng_cpg`, `sns_run_deng_demo` |
| `SNS_AirStep` | Planar, fixed-pelvis, two-link BPA air stepping | `sns_build_air_step`; `dev\run_air_step_20261025.m` |
| `SNS_Li_CPG` | Graded neural adaptation, scripted sensory smoke input | `sns_build_li`; `dev\run_li_smoke_20261024.m` |
| `SNS_W2L_CPG` | Saved W2L neural model in `results` | Audit its actual blocks and source before running |
| `SNS_SpinalNetwork` | Saved legacy 410-neuron neural export | Separate from the current Hunt candidate |
| `Gait2392_simbody_simscape` | Imported body model under `osim_import` | Body import, not an automatically coupled walker |

The planar air-step demo, standalone neural circuit, imported Simscape body,
and closed-loop Python Gait2392 walker are different demonstrations. Preserve
those distinctions in filenames, figures, and claims.

### 6.2 Run an existing graded model without rebuilding it

The new `run_sns_model_extended.m` loads an existing model, applies temporary
stop-time/mode settings, and saves `out` plus metadata in a new directory. It
does not rebuild or save the source model and refuses an existing result file.
Use a fresh MATLAB batch process for reproducibility:

```powershell
& $matlab -batch "addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape'); run_sns_model_extended('KneeReflexDemo','D:/GitHub/Bipedal_Robot/tmp/walker_runs/knee_example','5','normal');" *> "$run\matlab_knee.log"
```

Choose a new output folder when repeating. The helper preserves the saved solver
and step size. Its optional fourth argument is `normal`, `accelerator`, or
`rapid-accelerator`; compilation/compatibility must be tested for each model.
Two 0.02 s smoke tests succeeded here for the graded knee and Deng CPG demos.
The knee example also completed a five-second run. These are runtime checks,
not full behavior validation.

For interactive use, in MATLAB's Command Window:

```matlab
sns = 'D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape';
addpath(sns,fullfile(sns,'demos'),fullfile(sns,'results'));
open_system('KneeReflexDemo');
```

The helper refuses to run an already-loaded model with unsaved changes. Save
your intended edits to a separate model first; it will not close a user's
pre-existing model. A MATLAB batch process has its own loaded-model state.

### 6.3 Run spiking demos using a private library copy

The current shared `SNS_Library.slx` lacks `SpikingNeuron` and
`HybridSpikingSynapse`, so the saved spiking knee demo initially failed with
`Simulink:Libraries:MissingBlockInLib`. A model opening successfully does not
prove its referenced library is complete.

The new `prepare_spiking_library_extended.m` copies the library and its spiking
builder into a **new private directory**, builds the missing blocks there, and
leaves the shared library untouched. Run it in a fresh batch process:

```powershell
& $matlab -batch "addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape'); lib=prepare_spiking_library_extended('D:/GitHub/Bipedal_Robot/tmp/walker_runs/spike_library_example'); run_sns_model_extended('KneeReflexDemo_Spiking','D:/GitHub/Bipedal_Robot/tmp/walker_runs/spike_knee_example','5','normal',lib);" *> "$run\matlab_spiking.log"
```

Both directories must be new. The fifth runner argument explicitly selects
that library folder. This repaired route passed both a 0.02 s smoke test and a
five-second spiking knee run on R2025a. Small fixed-point gain-quantization warnings remain; preserve them in
the log and assess their effects before claiming numerical equivalence.
The smoke check does not validate every spiking demo or establish equivalence
to a graded network. Run each model's own tests after a circuit edit.

Do not run `sns_build_library`, `sns_build_spiking`, or other builders against
the shared library during someone else's editing session. These are writing
operations. A full rebuild may also need actuator/Deng modules; stage the
relevant builder family and inspect its dependencies rather than assuming one
builder constructs all block types.

### 6.4 Edit the neural network and select outputs

Double-click a neuron or synapse to edit its mask parameters, or make the change
in an isolated builder copy so it is reproducible. A graded synapse takes the
presynaptic voltage and emits conductance terms `[g; g*Esyn]`. The postsynaptic
neuron consumes applied current and its synaptic ports. Use `SynSum` for more
incoming synapses than the neuron provides. Do not add an external Vpost wire
to a one-input synapse or replace conductance summation with a current sum.
Preserve nF, uS, mV, nA, and the model's time units. AnimatLab and SNS-library
saturation conventions differ.

For the MATLAB GUI, add `gui` and run `app = SNSGui;`. Its JSON loader uses
`neurons`, `synapses`, `inputs`, and `outputs`; select nodes for properties and
use its Connect mode for source/target wiring. The GUI's muscle/sensory drawing
elements are not automatically exported as a complete neuromechanical model.

`sns_build_from_json(jsonPath,mdlName,slxDst)` builds that supported neural JSON
schema. The stock legacy exporter does not include all persistent-Na gating
parameters needed for the latest Hunt circuit. The connectome editor's
`nodes/edges` display JSON is another format. Do not rename fields and claim a
faithful port. Verify capacitance, leak, reversal potentials, gating equations,
delays, spike readouts, and motor scaling in addition to graph connectivity.

For results, load the helper's MAT file in MATLAB and list available outputs:

```matlab
r = load('D:/GitHub/Bipedal_Robot/tmp/walker_runs/knee_example/KneeReflexDemo_run.mat');
r.out.who
% Select an actual listed timeseries, or an element of r.out.logsout.
% Example after confirming that this model exported log_th:
if isprop(r.out,'log_th') || any(strcmp(r.out.who,'log_th'))
    th = r.out.log_th;
    figure; plot(th.Time,rad2deg(th.Data)); xlabel('Time (s)'); ylabel('Angle (deg)');
end
```

Different demos export different names. For a Dataset, use
`r.out.logsout.getElement('actual_signal_name').Values`. For interactive signal
selection, mark the desired wire for signal logging, give it a useful name,
simulate a **copied model**, and open Simulation Data Inspector. Choose signals
and axes there. For spiking neurons, plot spike output separately from membrane
voltage and filtered MN drive. Save a figure and the underlying MAT data.

### 6.5 Perturbation, acceleration, and training

An existing From Workspace/Step/input signal can apply a defined pulse; use
`Simulink.SimulationInput.setVariable` for a temporary variable override when the
model actually consumes it. Mechanical torque, neural current, tendon load,
and altered sensory reliability are different perturbations. Log pulse onset,
duration, direction/frame, magnitude, and recovery. Do not relabel a prescribed
sensory input as measured foot contact.

Normal mode is the reference. Test Accelerator with unchanged solver and step,
separating cold compilation from warm execution. Use Fast Restart only for
supported parameter changes; topology, sample-time, and some mask changes require
recompilation. Check trajectories after every speedup. A compiled mode that fails
is not evidence that a coarser solver is an equivalent substitute. Graphical
scopes and animation can dominate wall time, so disable them in benchmarks.

The Hunt four-step/five-stage pipeline currently runs in Python. The Simulink
demos do not acquire that curriculum through a menu or filename. A Simulink port
requires expected-afferent datasets, actual neuron dynamics, input/output parity,
and the same physical validation gates. Do not use `SNS_AirStep` as stage-4
ground evidence.

### 6.6 Animation and GIFs

`demos\sns_animate_demo.m` and `demos\sns_animate_demo_v2.m` render the supported
demo's saved MAT data as a two-dimensional mechanism. Read each function's input
signature and its expected signal fields before calling it. They are not
universal renderers for imported Gait2392 or AH001 CAD geometry.

For the graded knee MAT file made by section 6.2, this MATLAB example stages the
legacy renderer and its expected input name without replacing shared results:

```matlab
sns = 'D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape';
stage = 'D:/GitHub/Bipedal_Robot/tmp/walker_runs/knee_gif_example';
assert(~isfolder(stage),'Choose a new animation directory');
mkdir(fullfile(stage,'demos')); mkdir(fullfile(stage,'results'));
copyfile(fullfile(sns,'demos','sns_animate_demo_v2.m'),fullfile(stage,'demos'));
copyfile('D:/GitHub/Bipedal_Robot/tmp/walker_runs/knee_example/KneeReflexDemo_run.mat', ...
    fullfile(stage,'results','sns_demo_results.mat'));
addpath(fullfile(stage,'demos'),'-begin');
sns_animate_demo_v2('knee');
```

The GIF is under the staged `results/animations/v2` folder; the PNG still is
under `results/pictures/v2`. This legacy renderer uses a fixed 0.03 s frame
delay, so its playback is not necessarily real time. Its decorative BPA activity
is derived from joint angle, not recorded muscle activation. Use it as a joint-
motion illustration, and use the saved neural/force plots as the control evidence.

For a Simscape model, Mechanics Explorer animation requires a successfully
compiled Multibody model and its geometry. Playback or video export availability
depends on that environment. A bridge model can instead save full MuJoCo state
and use the corresponding Python replay. Record frame rate, time range, and
camera separately from the physics step. Keep a PNG still and source MAT/NPZ
with every defense GIF. Do not create GIF frames by rerunning a stochastic model
with different parameters.

### 6.7 MuJoCo bridge and Simscape imports

The bridge build driver is
`Code\Matlab\SNS_Simscape\mujoco_bridge\matlab\laptop_build_bridge.m`.
It was built and demonstrated on the laptop with R2025b. This easteregg2 checkout
has no discovered compiled MEX binary; do not present its bridge examples as
ready to run here. Rebuilding/installing a blockset is separate work from running
an existing simulation and can change MATLAB paths.

`e1_sns_reflex_loop.m` and `e2_cpg_open_loop.m` are development examples, not
read-only viewers: inspect their model creation/save calls first. In particular,
the latter saves a network and rebuilds an experiment model. Stage copies for
testing. Confirm actual controller/body sample times and actuator mapping.

For Simscape imports, verify the `.slx`, DataFile, geometry folder, and all
external references. A license probe is less useful than a small actual runtime
test. Do not interpret a licensing failure as a control failure. Successful
`smimport` or visible CAD geometry alone does not create muscles, sensors,
neural control, or valid contact.

## 7. Connectome block editor and browser display

Use `spinal\connectome_block_editor.html`, the v3.3 block editor, not the stale
small `connectome_editor.html`. The embedded template blob is primary;
`connectome_templates.json` is a fallback. `ben_connectome_canonical.json` is
the separate 110-neuron/216-edge canonical graph, not the full Hunt network.

The currently authorized local page is
`http://127.0.0.1:8765/connectome_block_editor.html`. If it is already open and
being edited, leave that page and its saved models alone. To serve a private
copy for a new session, copy the single self-contained HTML into a narrow
directory, then serve **that directory**, rather than the whole public repo:

```powershell
$editorDir = Join-Path $run 'editor'
New-Item -ItemType Directory -Path $editorDir | Out-Null
Copy-Item -LiteralPath "$spinal\connectome_block_editor.html" -Destination $editorDir
& $py -m http.server 8766 --bind 127.0.0.1 --directory $editorDir
```

That last command stays running; use another shell for simulation. Open
`http://127.0.0.1:8766/connectome_block_editor.html`. A different port is a
different localStorage origin, with a different saved-model collection. A file
URL and HTTP origin also have separate storage. Export/download important work
to disk; do not rely only on localStorage.

Select the intended tab/subnetwork, select a node or edge to inspect its actual
properties, and use the editor's focus view for readable blocks. The four
`Hunt Session1` template entries are display projections of their documented
snapshot. Import a new candidate into a new model/tab, not over a user's live
edits. Preserve edge signs, full numeric gains, and source/target identity;
compact canvas labels are rounded for readability.

Ben's downloaded exports normally land in
`Neuromechanical_Models\Mujoco_SNS_models`; look there first. Diagram exports
do not automatically update Python/Simulink source. A simulation edit needs an
explicit, validated importer or changes in the relevant builder/weights file.
After authorized testing, remove only the agent's clearly identified `TEST_*`
localStorage models from the tested origin. Never clear all browser storage.

## 8. Delegate bounded jobs to ZCode or another model

Use `ZCODE_SIMULATION_PROMPTS.md` as copyable job packets. Each packet states
**easteregg2**, exact repository/interpreter paths, inputs, permitted scope,
output directory, and stopping criteria. Give one worker one job and one output
folder. Start with analysis, logging, parity, and rendering jobs; they have
concrete checks and do not require redesigning the controller.

The useful immediate split is Prompt A (bounded frontal-feedback trials) and
Prompt B (force/current analysis), using separate folders. Prompt C's basin test
waits until **the same configuration** passes stages 1 and 2. Do not let a worker
chain failed stages or silently weaken thresholds. Read its evidence before
promoting a candidate. A coordinator can then authorize a restricted stage-3
standing/balance sweep and later actual ground stepping.

Keep optimization databases local, name studies explicitly, and cap workers
initially at four on easteregg2. Every worker process gets its own MuJoCo data,
network state, random seed, and output prefix. Avoid nested BLAS/process
oversubscription. Warm JIT compilation before timing. Verify finalists in fresh
processes; a repeated in-process attractor is not a robust solution. Use the
documented ±1% basin probe only after movement gates pass.

## 9. Minimum record and troubleshooting

For each scientific result, retain: machine/environment; source/model/config
hashes; all loaded weights; seed; frame adapter; neuron type and backend;
physics/neural step; initial conditions; assistance and external drive;
actual foot/harness forces; movement criteria and failed terms; wall-time scope;
raw recordings; plots; and replay provenance. Name any external teacher or
allocator. Keep one separate training record, with failures as well as successes.

| Symptom | First check |
| --- | --- |
| `--help` starts a long run | Legacy script uses handwritten argument parsing; read its source |
| Short run crashes in a max/min reduction | No samples in the protocol's scored time window; use its full schedule |
| Low RMSE but no movement | Inactive channels dominate, PF latched, or command/activation target confused |
| Neural alternation but crossed legs | Hip frontal/rotation torque, routing signs, load-dependent moment arms, width gate |
| Claimed support differs from GIF | Plot actual signed rig force divided by model weight, with transients shown |
| Too much foot load | Check once-per-contact accounting, frames, duplicated heel/toe sums |
| Fast run differs from reference | Gain refresh, reset state, sodium-channel tau, summation order, cached parameters |
| Spiking model opens but will not simulate | Missing library blocks or sample-time/readout mismatch; use the private-library recipe |
| MATLAB function resolves to an old copy | `which functionName -all`; remove broad archived paths |
| Browser edits disappear | Different origin/port, reload/import, or localStorage-only persistence; recover disk export |
| GIF looks good but force gate fails | It remains a failed physical trial; caption the limitation |

The acceptance target remains actual alternating ground stepping with measured
assistance. No harness-free walking has been established in this program's
records. Do not replace that statement with a successful animation or a numerical
speedup claim.

## W2L virtual walkers: direct AnimatLab → MuJoCo, added 31 October

Machine: **easteregg2**. The source-specific port, consumed configuration keys,
observer definitions, and complete experiment record are in
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\README.md`.
The original circuit has 80 neurons/150 chemical connections; the newer editable
project has 96/204. Select sources explicitly: these are separate circuits.
Their IntegrateFireSim equations are ported directly, with a Numba implementation
and a scalar reference; this runner does not offer the Gait2392 spiking selector.

### Run the original circuit with physical guides

This 20 s run uses a doubled neural time scale, rails lowered 20 mm, 30 mm
guide clearance, and a phase-triggered yaw perturbation. Use a new tag for each
experiment so saved trajectories remain available.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUN_GUIDE_30MM_CLOCK2.json' --duration 20 --tag my_guided_clock2_run
```

The resulting folder is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\my_guided_clock2_run`.
To compare slower timing, add `--neural-time-scale 4` and use a new tag. Explicit
CLI flags override configuration values. To compare identical calculations
without Numba, add `--reference`; exclude first-use compilation when timing.
The current matched 4 s air benchmark is 17.208 s reference versus 0.970 s
Numba, including physics and observers, excluding import, warm-up, and export.

### Render the saved motion

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\render_trial.py' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\my_guided_clock2_run' --start 2 --stop 5
```

This saves `SUPPORTED_STEPPING_REPLAY.gif` and `SUPPORTED_STEPPING_STILL.png`
inside that folder. Rendering replays saved states; it does not rerun control.
The runner logs root orientation/angular velocity, pelvis rotational energy,
whole-model COM, rail/guide/harness reactions, regional foot forces, clearance,
double support, belt-relative strides, and neural/muscle states. The record
defines filtering and qualification thresholds. Retain raw force channels.

### Edit and run the native AnimatLab companion

Open
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_GuidedT2_20261031\W2L_VirtualWalker_GuidedT2.aproj`
in AnimatLab for the 10 mm guide, doubled-clock version. The matching `.asim`
and `_Standalone.asim` are beside it. Hierarchical RG, hip, and knee/ankle pages
remain editable. Its separate record is
`D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Guided_20261031\EXPERIMENT_RECORD.md`.
Clock changes scale membrane and calcium-gate time constants, preserving physics
and muscle timing. After a GUI edit, export and verify the standalone graph:
the supplied newer project's GUI and old standalone export differ.

Current original-circuit 20 s supported tests pass the declared clearance and
alternation checks at clock scales 2 and 4; scale 6 loses consistent alternation.
The newer source circuit has a latched RG. A separately saved diagnostic with
symmetric zero tonic currents and only the two mirrored V3-to-contralateral-RG-IN
outputs disabled passes the bounded air and guided-ground checks; see the
experiment record for exact configuration and repeat status. This is a distinct
controller variant. These tests are kept separate from the Hunt training and
Gait2392 curriculum records.
