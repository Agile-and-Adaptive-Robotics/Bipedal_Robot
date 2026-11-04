# Fresh walker replays for the defense

Machine: **easteregg2**, 2 November2026. All three requested preserved configurations were run sequentially, with at most one heavy simulation in this agent's lane. No fitting, neural edits, new damping/contact/tendon mechanics, manuscript edits, or commits. Root owns the new CAD-foot revision separately.

| Model/condition | Fresh result | Measured net vertical assistance | Ground-clear unloaded swings |
|---|---|---:|---|
| Human Gait2392, graded biw2lvar-Hunt-PT + autogenic RC, stage4 |16s, exact saved capture/model identity |48.848356%BW during5–13s; horizontal restraint115.560892N RMS |5right/4left |
| Canonical NEWER108 graded virtual W2L, Rustuna winner |10s, all numeric trace arrays identical to retained fresh winner |39.237920%BW whole trial |9left/8right |
| NEWER W2L with192central+512sensory spiking cells and80graded cells, pooled-home comparator |10s, full trajectory plus central/sensory spike arrays identical |42.906486%BW whole trial |10left/10right |

All retain their recorded guides/rails/harness and flexed-knee guards. These are supported ground/treadmill steps. The Gait case steps in place. Spiking W2L moves+0.503636m relative to the belt while root world displacement is−0.442774m. Belt-relative progression is not independent overground walking.

These W2L conditions are **not a matched graded-versus-spiking experiment**. They retain different accepted initial poses and the spiking comparator's explicit pooled hip-physics option; the support difference cannot be assigned to neuron type alone. Both W2L examples have12uniarticular muscle routes; they are not the separate Li–Deng biarticular rat. No new walking success is assigned to s3k, Syn6, rat122/134, or central-spiking Gait.

## Open the GIFs and stills

### Human Gait2392

GIF: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\gait_hunt_rc\capture\plots\PRESENTATION_REPLAY.gif

Still: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\gait_hunt_rc\capture\plots\PRESENTATION_STILL.png

Model: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\gait_hunt_rc\capture\model.mjb

Controller/configuration: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\gait_hunt_rc\repeat.json

Program: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py

Fresh full-file identity: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\gait_hunt_rc\FRESH_COMPARISON.json

Assistance decomposition: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\gait_hunt_rc\capture\plots\SUPPORT_DECOMPOSITION_5_TO_13.json

The measured model has75.1646kg anatomical mass plus3.2kg existing numerical follower masses. The body-weight denominator is768.756726N. No mass was changed here. The 1,033 actual cells are graded; four additional logging names are aliases. The four-cell RF/HAM modules are retained; sartorius/gracilis are separate. All gains, weights, actuator force law, mechanics, and source guards are inherited unchanged from the accepted configuration.

### Graded W2L108

GIF: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\w2l108_graded\W2L_VirtualWalker_Newer\fresh_graded_ground\SUPPORTED_STEPPING_REPLAY.gif

Still: same folder, SUPPORTED_STEPPING_STILL.png

Model: same folder, model.xml

Settings and full fresh comparison: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\w2l108_graded\FRESH_CHECK.json

Source configuration: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA108_GROUND_WINNER.json

Program: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py

Executed external wrapper: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\source\fresh_w2l108.py

The wrapper consumes the expected canonical source hash and original full initial qpos. It changes only output_root and tag. Four mirrored ankle phase assignments and the saved recruitment/contact gains are retained. The actual controller is the direct translated source-equation/Numba runtime, not a generic replacement SNS network. The source clock is already included in its parameters; runtime scale remains1.

### Combined spiking W2L

GIF: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\w2l704_spiking\W2L_VirtualWalker_Newer\fresh_spiking_ground\SUPPORTED_STEPPING_REPLAY.gif

Still: same folder, SUPPORTED_STEPPING_STILL.png

Model: same folder, model.xml

Parameters, exact invocation inputs, implementation hashes: same folder, SPIKING_COMPARISON.json

Full motion/spike parity: same folder, FRESH_PARITY.json

Actual spike files: same folder, central_spikes.npz and sensory_spikes.npz. These store12central populations/192cells and16sensory populations/512cells, respectively. The remaining80cells are graded.

Neural checkpoint: D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_RESET_REST_CANDIDATE.json

Pose: D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\stance_home_v7\standing_home.npz

Program: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_spiking_central.py

Executed flags: --mode supported --duration 10 --configuration [checkpoint above] --sensory-spiking --initial-qpos [pose above] --pool-inherited-hip-physics --out D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\w2l704_spiking --tag fresh_spiking_ground. Use a new tag for a new trial.

This explicitly retains the saved hip-flexor/extensor sigmoid field averages and legacy hip-extensor Ib conversion. It does not silently alter source AnimatLab files. It is the existing pooled-home comparator, not an unqualified all-spiking circuit or a new training result.

## Runtime and media

Scientific Python: D:\Anaconda\envs\myo\python.exe, MuJoCo2.3.7. Gait uses the verified sparse/accelerated SNS path and compiled body repairs; W2L uses the verified preallocated Numba source neural/Hill equations. Gait16s took100.67s including its instrumented evaluation construction; the spiking W2L warmed10s loop took5.01s. Construction/JIT/rendering scopes differ, so these are individual scoped timings, not a cross-engine speed benchmark.

Fresh Gait GIF:107frames at20fps (~5.35s). W2L GIFs:51frames at80ms (~4.08s). Each replays actual saved states at the recorded cadence. All three stills were visually inspected. Colors/cameras/materials and explanatory banners are display changes only.

Machine-readable registry and file hashes: D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\MODEL_REPLAY_MANIFEST.json. This is an evidence/GUI-input registry, not a completed universal launcher or a license to overwrite existing outputs.

All owned jobs are complete. Replay programs reject existing result tags or captures; use a new owned destination for future runs. Do not launch the old time-bounded campaign scripts.
Exact spiking command on easteregg2, using a new tag to preserve results:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_spiking_central.py' --mode supported --duration 10 --configuration 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_RESET_REST_CANDIDATE.json' --sensory-spiking --initial-qpos 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\stance_home_v7\standing_home.npz' --pool-inherited-hip-physics --out 'D:\Users\Ben Bolen\ToolchainDefense_20261102\other_models\w2l704_spiking' --tag user_repeat_01
```

The two private fresh_* wrappers preserve their already-used fixed destinations and intentionally refuse overwrite. They are execution records, not generic GUI launch endpoints. The underlying W2L run() function accepts output_root; the GUI should supply a new run directory/tag. Gait's existing configuration/capture recipe in the earlier campaign's RESULTS_AND_FILES.md likewise supplies a new destination without changing the saved controller.
