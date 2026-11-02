# Virtual W2L walkers: native circuit and MuJoCo port

## Opt-in spiking experiments on easteregg2

The separate current record is `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\EXPERIMENT_RECORD.md`; the checkpoints beside it record exact hashes, frozen air/ground parameters, measured forces, and fresh-process parity. Raw outputs stay outside this public repository. The graded results below remain unchanged.

The combined variant has 704 actual spiking cells (192 RG/PF plus 512 II/contact) and 80 graded cells. At its earlier supported home it passes 10/10 swings, fresh saved-array/event parity, and 8/8 bounded 16-scalar perturbations; see `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\COMBINED_SPIKING_CHECKPOINT.json`. A subsequently verified CoM-over-feet home reduces net assistance to 42.25% BW with unchanged learned gains (9/9 swings, fresh exact); see `COM_HOME_CHECKPOINT.json` in the same external folder. This corrects an initial CoM about 69 mm behind the feet. Static feasibility forces are never applied to the SNS.

An explicit source-based pooled hip-force/legacy-sensor comparator passes 10/10 at 42.91% BW and +0.504 m belt-relative progression; fresh exact, plus 8/8 new-home closed-loop ±1% checks on 16 learned scalars. All eight retain forward belt-relative progression (0.367–0.519 m/10 s), flexed knees, and 41.04–43.63% BW support. Its actual source field changes are in `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\POOLED_PHYSICS_CHECKPOINT.json`. Both new-home variants still fall after complete horizontal-rail release. Neither is independent walking or a completed physical curriculum, and the original robust checkpoints are preserved.

The sensory-spiking W2L hybrid uses 512 LIF cells and 92 graded cells and passes the supported-ground gate plus 8/8 bounded encoding/gain perturbations. A distinct adaptive central-spiking variant replaces 12 RG/PF cells with 192 actual Ca/LIF cells; 96 other source cells remain graded. The original central winner passes nominal air/ground and 7/8 perturbations. A separately saved RG reset-margin candidate passes both air checks, 10/10 qualified ground swings, a byte-identical fresh repeat, and 16/16 bounded perturbations, including eight held-out seeds. Its exact configuration and reproduction command are in `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_RESET_REST_CHECKPOINT.json`. Nominal net support is 46.02% BW. Lowering the rails another 5–15 mm retains stepping but increases measured assistance, so it is not successful support weaning. These are explicit neural variants of the same canonical NEWER108 source, not additional active AnimatLab model families or a completed physical curriculum.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_spiking_w2l.py' --case hybrid --tag YOUR_FRESH_SENSORY_TAG
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_spiking_central.py' --mode supported --configuration 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_RESET_REST_CANDIDATE.json' --tag YOUR_FRESH_CENTRAL_TAG
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_spiking_central.py' --mode supported --configuration 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_RESET_REST_CANDIDATE.json' --sensory-spiking --initial-qpos 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\stance_home_v7\standing_home.npz' --pool-inherited-hip-physics --tag YOUR_FRESH_POOLED_TAG
```

`run_spiking_gait.py` tests sensory encoding on the corrected Hunt Gait2392 ground checkpoint. `run_spiking_standing.py` tests a separately calibrated syn6 standing controller; it must not be presented as the same Hunt controller. MATLAB wrappers `verify_existing_spiking_reflex.m` and `verify_spiking_reflex_pulse.m` load an isolated external library snapshot containing the existing spiking extension. None of these options changes the canonical GUI into a spiking model.

**Latest supported-ground result:** the actual108-neuron model now passes after correcting its physical starting posture with bounded near-flat IK. Source-phase alone yields7/8 clear swings,92.9% alternation; four explicitly mirrored PF→ankle-MN target swaps yield9/8,100%,42.24%BW net assistance. A24-trial Rustuna calibration retains9/8,100% while reducing whole-trial net assistance to39.24%BW; the selected result repeats byte-identically in a fresh process. Knees remain flexed. See the top of EXPERIMENT_RECORD.md for forces, controls, source hashes, and bounded claims. These runtime variants have not yet been embedded back into the editable AnimatLab geometry/network.

Current repeat command on **easteregg2**:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --newer --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\RUSTUNA108_GROUND_WINNER.json' --duration 10 --tag YOUR_UNIQUE_TAG
```

The consumed JSON records actual initial qpos, four ankle phase swaps, shared PF→MN recruitment, and absolute bilateral toe-contact→ankle-MN gain. It preserves the source30mm passive guides/rails−20mm and uses runtime clock1 because source time constants already include clock2. No inverse-force targets or static tension preload are used in these SNS trials. Rat→Gait transfer is paused.

**Current model is canonical NEWER108/216**, with12 new independent homonymous II relays plus the four existing phase-II channels. Source clock2 is already in its actual parameters: use runtime neural scale1. Actual108 baseline and matched RG/PF/MN −10 nA/100 ms fast-MuJoCo tests are completed, with joint motion and measured current waveforms. See the latest section at the top of EXPERIMENT_RECORD.md and `runs\CANONICAL108_STIMULUS_TRIPLET.json`. The separate rat controller has120 neurons; W2L plots must not be labelled rat.

**Use the canonical NEWER model for the paper and current work:** `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj` and its synchronized standalone. NEWER has bilateral RG, separate bilateral heel/toe contacts, four legacy phase-II channels, and twelve independent homonymous II channels. The original80 model is only a historical comparison. Six generated sibling families were removed after immutable input/evidence archival. The successful MuJoCo runtime configuration above includes a physical home/phase/tuning overlay that is not yet embedded in the reference GUI.

The prior NEWER96/204 diagnostic has actual passing air and supported-ground results. One positive V3-output scale 0.01 candidate also passed (ground8/8 clear swings,100% alternation,73.53% BW assistance). Its snapshot precedes the current reference-numeric canonical update, so it is a benchmark, not a passing claim for an untested changed export. See the top of EXPERIMENT_RECORD.md for exact settings and source scope. Six historical source snapshots live under ignored `runs\source_inputs\`, with SHA-checked lookup; they are not alternate active models.

Machine: **easteregg2**. Python: `D:\Anaconda\envs\myo\python.exe`.
This is a separate experimental package; it does not replace the Gait2392 controller or modify the supplied AnimatLab mechanics.

The two requested source walkers are imported directly from their standalone AnimatLab XML, including bodies, passive support wings, horizontal rails, treadmill, joint limits, twelve uniarticular muscle routes, and twelve stretch-receptor routes. Neural outputs drive the original sigmoid/Hill actuator dynamics in a closed feedback loop. There is no posture allocator, force-target optimizer, or replacement generic W2L network. `EXPERIMENT_RECORD.md` gives results and limits.

## Historical original80 checkpoint, not the paper selection

PowerShell, from any directory:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\ORIGINAL_TIMING2_SAGITTAL.json' --duration 10 --tag my_supported_repeat
```

This selects the original 80-neuron/150-edge circuit, shared neural time scale 2, horizontal rails lowered 20 mm, and an explicit continuous sagittal harness. It is rail-supported treadmill stepping, not independent standing or overground walking. The harness has no vertical component; the rails still carry about 86% of body weight.

## Physical contact guides with free play

The optional guide carriage has eight collision boxes, 2 kg mass, and a passive fore–aft slider with finite ±1.9 m travel. It moves through contact, not by following or teleporting to the robot. Side walls, wing posts, and overhead caps allow actual clearance before contact. The robot root stays free. `contact_guides.py` defines the apparatus; `verify_guides.py` checks contact-free motion inside its slack and restoring force/torque signs outside it. These are engineering assumptions, not calibrated hardware properties.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --neural-time-scale 2 --rail-height-offset -.02 --sagittal-harness 0 --guide-clearance .01 --yaw-perturb-nm 5 --duration 10 --tag my_physical_guides
```

Here `.01` means 10 mm neutral side/fore–aft slack. The overhead gap is initially also 10 mm and grows as the model settles onto the horizontal rails; `GUIDE_ASSESSMENT.json` reports the settled vertical separation separately. No continuous harness is enabled. The optional disturbance is a 5 N·m world-Y torque for 100 ms, triggered at the first left RG-flexor crossing of −59 mV after 2 s. Its actual onset and torque are logged. Omit `--yaw-perturb-nm` for an undisturbed run.

`ORIGINAL_TIMING2_PHYSICAL_GUIDES.json` is the versioned undisturbed 10 mm guide checkpoint. Explicit CLI flags override values in a configuration file. For an extended-duration test, `--duration 20 --belt-motor-end-s 20` keeps the belt moving until 20 s; the source default remains a motor stop at 10 s. `long_guide_comparison.py` records matched 20 s cases and checks finite belt, rail, and carriage margins. `newer_guide_pairs.py` separately tests the existing newer 96/204 candidates; their reduced body rotation is not a neural stepping pass.

`bounded_guides.py` runs the finite comparison; `extended_guide_probes.py` adds deliberately faster and slower clocks. All use the same original neural weights and physical muscle parameters. Larger neural time scale means slower neurons, not a larger numerical timestep.

## Inspect outputs

Every run saves its actual `model.xml`, `RESULT.json`, and `trajectory.npz` under `runs\W2L_VirtualWalker\<tag>\`. Raw data are gitignored. Results include source and implementation hashes, all neural voltages and tensions, pose, center of mass, angular velocity, pelvis rotational energy, regional foot forces, center of pressure, foot clearance, support reactions, and guide reactions. NaN center of pressure during no contact is intentional; it is not a simulation instability.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\assess_trial.py' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\my_physical_guides'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\draw_trial.py' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\my_physical_guides'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\render_trial.py' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker\my_physical_guides'
```

The renderer replays the measured trajectory; it does not rerun or synthesize neural control. Its title is selected from the actual support condition and reports measured vertical assistance. `draw_trial.py` shows short RG/PF/MN and joint-angle windows. Force plots retain raw contact chatter.

## Historical source imports and checks

The initially supplied newer standalone had96 neurons/208 edges and differed from its editable96/204 project. Its contact polynomial was experimentally measured injecting20 A, not20 nA. Those historical snapshots and their1 nA/20 nA diagnostics are archived under ignored `runs\source_inputs\`; their original generated sibling folders were removed. The active `--newer` source is now the synchronized canonical108/216 model described above.

`validate_kinematics.py`, `route_audit.py`, and `functional_sign_audit.py` check geometry, actual moment-arm signs, and neuron-to-muscle connections. `verify_checkpoints.py` checks air stages, a fresh-process repeat, and a 1 ms physics comparison. `benchmark_port.py` compares identical scalar Python and Numba equations, excluding loading, warm-up, saving, and rendering. It is not an AnimatLab-versus-MuJoCo or SNS-Toolbox dense-backend speed comparison.

The source adapter supports the supplied nonspiking IntegrateFireSim circuit. Unsupported spiking, delayed, electrical, noisy, and extra ion-channel configurations are rejected rather than replaced silently. The canonical anatomical frame is x anterior, y superior, z right; native source Z is reflected explicitly. Root-Y below 0.6 m, any anatomical knee hyperextension, or nonfinite state stops a run.

## Historical96 RG diagnostic

The newer source's RG latches even without physical feedback. Tied 0/0 tonic currents alone do not fix it. A bounded diagnostic additionally sets only the two mirrored V3→contralateral RG-extensor-IN output conductances to zero, retaining c1, all neurons, source intrinsic parameters, and both existing 10 nA/10 ms startup kicks. Source tonic 2/0 plus V3-zero does not restore all four RG cells either. The successful branch is explicitly the combined **tonic0/0 + V3-output0** condition; it is not an optimized or unchanged source model.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --source 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\source_inputs\W2L_VirtualWalker_Newer_GUI1nA_Standalone.asim' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\NEWER_GUI1NA_TIED0_V3ZERO_DIAGNOSTIC.json' --mode supported --guide-clearance .03 --duration 10 --tag my_newer_diagnostic_ground
```

For its air checks, use the same source/config, `--mode air --duration 8`, and omit the guide flag. Add `--feedback-off` for the deafferented check. The recorded cases passed both air gates and a bounded supported-ground gate (6/8 clear swings, 84.6% alternation, 35.2% BW mean assistance). Ground onset is transiently asymmetric; this is not complete Hunt training or independent walking. Source/default V3 gain remains 1 and was verified to reproduce the old trajectory bit-identically. All altered settings are saved, consumed by the loader, and labeled on rendered media.

Completed diagnostic outputs are at `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer_GUI1nA\newer_v3_diag_v1_ground_guides\`. The source native companion and its independent scope are documented in EXPERIMENT_RECORD.md. The default source models remain available separately.

## Body-balance and timing experiments

Current separate record: `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\BALANCE_PATHWAY_RECORD.md`; machine-readable results: `BALANCE_CHECKPOINT.json` in the same folder. The two new spiking task-feedback projections preserve air/supported stepping but do not pass rail release. Whole-neural-clock tests prove unchanged deafferented MN waveforms at matched biological times; faster timing passes air/ground but does not solve release, while slower timing fails the physical air gate. These behavioral clocks are distinct from computational speed-ups.

The frozen704-spiking/80-graded controller at the explicit0.2m/s belt condition passes10/9 swings at36.48%BW net assistance, a byte-exact fresh repeat, and8/8 held-out ±1%16-scalar checks. Geometry, pooled source Hill parameters, source ROM, learned gains, and physical guides are unchanged. This is an environment-condition finding, not neural training or harness weaning. Reproduce on easteregg2:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_spiking_central.py' --mode supported --configuration 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\CENTRAL_RESET_REST_CANDIDATE.json' --sensory-spiking --initial-qpos 'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\stance_home_v7\standing_home.npz' --pool-inherited-hip-physics --belt-speed 0.2 --duration 10 --tag my_pace02_check
```

## Source motor-output calibration

`D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\motor_output_fit_v1\MOTOR_OUTPUT_RECORD.md` and its checkpoint preserve the actual source-current/force/torque fits and bounded air→ground tests. Stronger calf outputs cause hopping; selective dorsiflexor output loses support transfer. A fitted knee-flexor gain2.3 gives a nominal34.39%BW result with exact fresh state/event/physical replay, but fails one of eight17-scalar checks and is NOT promoted. The robust0.2m/s baseline above remains preferred. Runtime optional output knobs default to their original values; actual gains/thresholds are saved. Source sigmoid laws are unchanged.

## Phase/load termination timing audit

The explicit experimental event pathway is documented at `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\load_termination_v1\CHECKPOINT.json`. It adds32 actual spiking cells and4 chemical connections to inherited PF inhibitory cells, with force-qualified, self-expiring pulses and RG phase gating. Zero-path numeric/spike parity and complete event/filter state replay pass; source models remain unchanged. Recorded timing shows the cue arrives after the relevant bursts have ended, so this prototype is not promoted and its gain sweep was not resumed. The verified704-spiking/80-graded reference above remains preferred.
