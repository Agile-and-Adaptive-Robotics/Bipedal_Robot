# Virtual W2L walkers: native circuit and MuJoCo port

**Current model is canonical NEWER108/216**, with12 new independent homonymous II relays plus the four existing phase-II channels. Source clock2 is already in its actual parameters: use runtime neural scale1. Actual108 baseline and matched RG/PF/MN −10 nA/100 ms fast-MuJoCo tests are completed, with joint motion and measured current waveforms. See the latest section at the top of EXPERIMENT_RECORD.md and `runs\CANONICAL108_STIMULUS_TRIPLET.json`. The separate rat controller has120 neurons; W2L plots must not be labelled rat.

**Use the canonical NEWER model for the paper and current work:** `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj` and its synchronized standalone. NEWER has bilateral RG, separate bilateral heel/toe contacts, and four physical II channels (left/right hip and knee). The original80 model is only a historical comparison. Generated sibling model families are being removed after immutable input/evidence archival.

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

## Other source branches and checks

`--newer` imports the supplied newer native standalone graph (96 neurons/208 edges), which differs from its editable project. The GUI-authoritative companion uses `--source` with its full standalone path and has 96 neurons/204 edges. The supplied newer contact polynomial was experimentally measured injecting 20 A, not 20 nA. Separate 1 nA and 20 nA candidate companions preserve distinct graph and gain choices; neither has passed this package's ground gate. Do not silently substitute them for the original source checkpoint.

`validate_kinematics.py`, `route_audit.py`, and `functional_sign_audit.py` check geometry, actual moment-arm signs, and neuron-to-muscle connections. `verify_checkpoints.py` checks air stages, a fresh-process repeat, and a 1 ms physics comparison. `benchmark_port.py` compares identical scalar Python and Numba equations, excluding loading, warm-up, saving, and rendering. It is not an AnimatLab-versus-MuJoCo or SNS-Toolbox dense-backend speed comparison.

The source adapter supports the supplied nonspiking IntegrateFireSim circuit. Unsupported spiking, delayed, electrical, noisy, and extra ion-channel configurations are rejected rather than replaced silently. The canonical anatomical frame is x anterior, y superior, z right; native source Z is reflected explicitly. Root-Y below 0.6 m, any anatomical knee hyperextension, or nonfinite state stops a run.

## Explicit newer RG diagnostic

The newer source's RG latches even without physical feedback. Tied 0/0 tonic currents alone do not fix it. A bounded diagnostic additionally sets only the two mirrored V3→contralateral RG-extensor-IN output conductances to zero, retaining c1, all neurons, source intrinsic parameters, and both existing 10 nA/10 ms startup kicks. Source tonic 2/0 plus V3-zero does not restore all four RG cells either. The successful branch is explicitly the combined **tonic0/0 + V3-output0** condition; it is not an optimized or unchanged source model.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\run_virtual_walker.py' --source 'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_GUI1nA_20261031\W2L_VirtualWalker_Newer_GUI1nA_Standalone.asim' --config 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\NEWER_GUI1NA_TIED0_V3ZERO_DIAGNOSTIC.json' --mode supported --guide-clearance .03 --duration 10 --tag my_newer_diagnostic_ground
```

For its air checks, use the same source/config, `--mode air --duration 8`, and omit the guide flag. Add `--feedback-off` for the deafferented check. The recorded cases passed both air gates and a bounded supported-ground gate (6/8 clear swings, 84.6% alternation, 35.2% BW mean assistance). Ground onset is transiently asymmetric; this is not complete Hunt training or independent walking. Source/default V3 gain remains 1 and was verified to reproduce the old trajectory bit-identically. All altered settings are saved, consumed by the loader, and labeled on rendered media.

Completed diagnostic outputs are at `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\runs\W2L_VirtualWalker_Newer_GUI1nA\newer_v3_diag_v1_ground_guides\`. The source native companion and its independent scope are documented in EXPERIMENT_RECORD.md. The default source models remain available separately.
