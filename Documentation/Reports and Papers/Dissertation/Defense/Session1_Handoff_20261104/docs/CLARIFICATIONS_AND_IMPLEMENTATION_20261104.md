# Controller/toolchain clarifications — Session 1, easteregg2

4 November 2026. Defense/editor handoff; no manuscript edits. Latest acceleration implementation, exact numerical comparisons, timings, source hashes, and configs: D:\Users\Ben Bolen\NeuroAcceleration_20261104\EXPERIMENT_RECORD.md.

## Orin implementation and dissertation coverage

Actual code: D:\GitHub\Bipedal_Robot\Code\Hardware\walker_io\README.md, orin_daemon.py, protocol.py, sensors.py, afferents.py, valves.py, and lab_console.py. README explicitly labels skeletons. Sensor/valve hardware TODOs and the compiled SNS hook remain incomplete; sns_rt.py is a proposed integration hook, not an implemented exported circuit. The 1 kHz description is a target, not measured performance.

Local source D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\70-overall.tex:95 states the Teensy/Orin hardware objective. No explicit walker_io implementation description was found in local 60-sns.tex or 95-AppendixD.tex. This is a local-source finding, not a live Overleaf audit. Defense description: a software interface scaffold linking sensory encoding, a future compiled SNS step, valve abstractions, and watchdogs, awaiting hardware integration.

## Six muscle families and conversion roles

The static bilateral surrogate fit covers bifemsh (biceps femoris short head), vas_med, vas_int, vas_lat, soleus, and tib_ant: six named families, twelve bilateral actuators. These are not six synergies or all-muscle dynamic validation. D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\shared_force_law_functional_v1.json and MOMENT_FORCE_AUDIT.md contain the scope/results.

There are actual source-matching gaps, including active/passive force differences that peak-force matching hides. This does not establish poor work by a particular agent. repair_robotbody_muscles.py builds the special robot-route OpenSim copy; it does not calibrate a human conversion. convert_to_mujoco.py contains a laptop-specific bone-mesh folder, C:\Users\Ben\Documents\GitHub\myoconverter\models\osim\Gait2354Simbody\Geometry, plus quick-preview conversion settings. 'Laptop geometry path' meant this asset-folder setting, not an observed foot scale error. Separate source corrections, actuator-law fitting, numerical regularization, and harness settings so each scientific condition is reproducible.

## FSA, Hunt, and Ben's curriculum

Szczecinski 2020 Table 2 derives compatible GLIF spiking parameters for specified encoded quantities and dynamic functions. Use it in place of corresponding graded-neuron derivations when replacing/adding spiking components. It can initialize the RG/PF/other functional modules it actually represents, and reduce subsequent fitting. It does not automatically implement an ionic burster or replace body-level training. Source: https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2020.577804/full.

Ben intended his stages to incorporate calibration as well as validation. Stage 1 includes intrinsic rhythm calibration and deafferented/suspended verification. Stage 2 includes expected-afferent timing, feedforward recruitment, and sensory waveform fitting in ordered blocks, followed by closed-loop afferented air motion. Feedforward mapping is also needed to observe meaningful mechanical motion in stage 1. Hunt A-D can therefore have fitting substeps before and within these environments; they are not rigidly separate physical stages. The current physical_eval stage selector only evaluates saved fitted parameters. intrinsic.py, sensory_timing.py, and motor_fit.py supply separate fitting components; a single orchestration button combining fitting plus physical acceptance is not implemented.

## Stage 3 postural pathway: current implementation and required replacement

Current hunt_training\balance.py sends directional PT_BAL commands and inhibitory relays directly to muscle motor neurons, using a fitted muscle-level task basis. This is NOT Ben's requested postural premotor circuit. Turning off RG tonic drive in standing does not create that circuit. Some joint-axis circuitry elsewhere does not remedy this direct pathway.

Required functional architecture:

CoM/trunk position and velocity + orientation + foot load/contact -> descending directional posture commands -> per-joint antagonistic postural premotor populations and inhibitory interneurons -> shared muscle MN/reflex subnetworks -> actuators.

In stage 3, phasic locomotor RG/PF commands are gated off and the postural premotor branch supplies the command function. Existing local afferent/reflex processing remains declared, with individual pathway signs checked. In stage 4, posture and locomotor commands blend at the premotor level with bounded saturation and recorded mode weights; do not sum two uncontrolled drives directly at MNs.

The postural input is not simply copied to every antagonist. Each directional command excites the correcting functional pool, inhibits the opposing pool through its interneurons, and may recruit coordinated pools across legs. Contralateral adductor/abductor combinations depend on which foot bears load and actual moment-arm/body-response signs. Test the mapping from loaded Jacobians and perturbation response rather than assigning signs from muscle names alone. Keep x anterior, y superior, z right explicit; yaw and sagittal/frontal rotation axes must not be confused.

Design contracts before claiming this works:
- Named represented errors, reference posture, sensor units, delays, gains, bandwidth, saturation, and mode gates.
- Explicit per-joint two- or four-cell premotor groups with named excitatory/inhibitory edges, source justification, and mirrored parameter families.
- Positive/negative directional-command probes with RG/PF gated off; verify correct torque sign, opposite-pool response, no excitation of both antagonists by accident, and declared coactivation floor.
- Separate intra-leg and cross-leg tests under bilateral, left-only, and right-only support. Gate unloaded-leg recruitment appropriately.
- Body verification logs premotor cells, MNs, joint motion, three-axis contact/support forces, harness wrench, knees, mesh clearance, and recovery.

No replacement postural module or new standing result is claimed in this pass. This document resolves the intended architecture and records the existing implementation gap for the next controller revision.

## What the +/-1% tests actually changed

Gait's eight-case suite uses ten shared scalars, tied across sides: soleus_uS, tib_ant_uS, knee_flexor_multiplier, Ia_homonymous, II_homonymous, Ib_autogenic, Ia_reciprocal_pool, rg_tonic_current, PT_kp, PT_kd. Each episode gets independent uniform factors in [0.99,1.01], fixed for that episode, cloned from the saved condition. This is neither changing every connection nor learning/noise during a run. Source: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\finalist_basin.py (SELECTED). Records: D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\training\finalist_ground9_rc_basin_v1\PLAN.json and BASIN_RESULT.json. Its reported 8/8 applies only to those tests.

CAD's 5/8 suite was different: shared engineering gains, cadence, motor timing, contact/clearance parameters, and three passive-joint families. It was not a test of every neural synapse. Records: D:\Users\Ben Bolen\CAD3DOF_20261102\contact_handover_20261103\EXPERIMENT_RECORD.md and SELECTED_CONDITION.json. CAD remains paused; no new joint/tendon changes or CAD tests here.

Parameter tolerance, ongoing neural noise, sensor errors, and physical disturbances are different experiments. Preserve analytical design weights during environmental tests. Keep an optional, separately labeled parameter-tolerance check to measure distance from a bifurcation; a failure warrants investigation rather than changing the official winner or declaring all networks unstable.

## Robustness and ongoing noise

The broader suite remains incomplete. Recent selected-scalar checks are not full IK-bank generalization, push recovery, or stochastic robustness. The 43-cycle gait bank supplies reference material, but current Hunt target preparation uses selected IK/GRF data. Hold out whole participants, speeds, and trials; do not randomly split neighboring samples of one cycle between fit and test.

Neural noise is an established modeling option, especially in spiking/rate-code comparisons. Szczecinski 2020 explicitly tests reset noise and populations; that is a concrete justification for matching its noise model, not a claim that every deterministic circuit must contain noise. Current W2L sensory initial-phase staggering is deterministic, not ongoing random noise.

A 1% fixed conductance change and 1% input-current noise are different perturbations. Declare physical units, baseline range, correlation time, independent versus common noise, and injection site. A useful colored-current test is eta_next=exp(-dt/tau)*eta + sigma*sqrt(1-exp(-2*dt/tau))*N(0,1); sigma is stationary current standard deviation. Use recorded/seeded streams, clean baseline, multiple severities and independent seeds. Do not perturb every conductance or threshold by default. For spiking models evaluate decoded output, gait/recovery, and failure rates, rather than requiring identical spike times in noisy runs.

Small noise can improve encoding/desynchronization or harm coordination; a small parameter-test pass does not establish either. Conversely, a failed +/-1% case does not prove a particular smaller current noise will fail. CAD's limited local margin needs improvement; Gait's 8-case pass is promising within its narrow scope. No broad ongoing-noise body suite was added or run in this acceleration pass.

## Scaled Gait2392 reuse and deployment

Ben is correct: same-topology Gait2392 models with the same muscle names/side/unit interface should reuse the same circuit and bindings. Scaling alone does not require a new topology. Geometry, strength, inertia, encoding ranges, and delays may need recalibration and physical checks for comparable behavior; distinguish software compatibility from a guarantee that all old weights yield the same gait.

Current acceleration helps fixed-memory CPU execution and later code generation. It does not yet implement GPU kernels or complete hardware I/O. Benchmark an ARM CPU controller first, then GPU tasks where workload warrants it. Full sensor-to-valve latency and jitter, startup/JIT removal, and watchdog behavior still need hardware tests. MuJoCo physics need not run in every physical robot control tick.

## Dissertation coverage and the two model-construction stories

Verified in local ProofFinal sources on easteregg2. Chapter 7, 70-overall.tex:95, mentions the hardware objective, Orin, Teensy-class hardware, and sensory channels. It does not describe the implemented walker_io scaffold. No walker_io, orin_daemon, or sns_rt description was found in Chapter 6 or the appendices. Ben's request to include the scaffold was not fulfilled. Do not substitute the intended hardware architecture for the code that actually exists. The README and individual modules are at D:\GitHub\Bipedal_Robot\Code\Hardware\walker_io\; hardware operations and SNS integration remain incomplete. Its illustrative loop budget is not a benchmark.

Chapter 6 contains pieces of the intended story: headless/GUI/API/SDK at lines 42-50; model-construction/interface diagram at 77-83; SNS-Simscape blocks at 165-188; CAD-derived Simscape mechanics at 210; runtime MuJoCo-Simulink bridge at 225-227; fitting/curriculum at 229 onward. However, this is not yet the complete worked walkthrough in GAIT_REPAIR_TRAINING_AND_DEPLOYMENT_20261103.md. A tools-centered narrative should connect the operations to examples rather than simply list platforms and results.

Recommended linked examples:
1. Human-model route: choose scaled Gait2392 and reference gait -> convert with MyoConverter -> verify/repair coordinates, muscle paths, moment arms, and force behavior against that exact OpenSim source -> bind an existing named controller -> choose declared actuator/passive-loss profiles -> apply appropriate graded/spiking FSA and body scaling -> calibrate through Hunt and the curriculum -> run accelerated headless trials -> display aligned RG/PF/MN/joint/load traces and a supported-motion animation -> evaluate held-out conditions.
2. Robot-model route: construct CAD geometry, joints, inertias, and actuator attachments -> import mechanics into Simscape Multibody -> assemble the SNS controller with the custom graphical editor and SNS-Simscape blocks -> connect activation/pressure/force and afferent interfaces -> test graphically with mechanical/neuronal measurements -> use the MuJoCo-Simulink runtime bridge or a separately checked MuJoCo model/controller build for accelerated trials -> compare the same named signals and mechanical criteria.

The neural editor, SNS-Simscape block library, connectome JSON tooling, and runtime bridge should each be named accurately. The bridge exchanges commands/measurements; it is not a lossless arbitrary Simscape-mechanism-to-MuJoCo exporter. A separately constructed CAD MuJoCo model is currently a separate conversion/build route. An example should identify which route it actually used. New universal GUI orchestration remains distinct from the components already demonstrated.

Scaffold placement recommendation: one concise Chapter 7 hardware-transfer paragraph describing sensor/afferent/valve/protocol/watchdog interfaces plus the unfinished SNS/hardware hooks; an Appendix D implementation table for filenames, units, and status. This supplements the worked Chapter 6 simulation story. No manuscript text was changed during this audit.

## Rybak drives: paper-to-code comparison

Primary sources: Rybak et al. 2024, eLife 13:RP98841, https://elifesciences.org/articles/98841 (Figure 4 and Table 1, Equations 13-15); the attached hemisection schematic also appears as Figure 1 in the 2025 follow-up, https://elifesciences.org/articles/103504, which builds on the 2024 model.

Paper roles:
- alpha_L/alpha_R: descending excitation of ipsilateral and contralateral flexor centers and selected ipsilateral interneurons. They also presynaptically reduce ipsilateral sensory-feedback effects.
- gamma_L/gamma_R: ipsilateral extensor excitation. In the 2024 intact parameter set, gamma is held constant, while alpha is varied with the simulated locomotion condition.
- Sensory pathways and drive interact. Equations 14-15 scale feedback as k*E/(1+k_PS*alpha_ipsi). Removing drives in the transected condition also removes their feedback suppression. That experiment is not identical to an intact voluntary walking model or simply multiplying all excitation by zero.

Current Gait families (code inspection, not new simulations):
- s3k: shared DRIVE -> both RG_E/RG_F. Optional GAMMA_s -> both ipsilateral RG centers; ALPHA_s -> ipsilateral MN pools.
- w2lvar: same shared/optional per-side drive contract, with its own additional commissural circuitry.
- syn6: same contract. A separate default-off syn6_brainstem branch sends shared-DRIVE-fed BS_GAMMA to RG_E and BS_ALPHA to RG_F/selected commissurals; it still does not implement four independently scheduled side-specific commands. Its source records earlier latching when this extra excitation was stacked.
- In all three, runner.py supplies optional GAMMA_c_s=drive*gamma_s and ALPHA_c_s=drive*alpha_s (apart from the explicit alpha-standing override). Thus available side-specific amplitudes do not imply independently controlled drive time courses or the Rybak targets.
- Current Hunt-PT+RC stage-4 condition retains a shared DRIVE motor-pattern gate and common rg_tonic_current=1.2 nA on both E/F centers of both sides, plus a brief initial antiphase kick. It is not the Rybak alpha/gamma architecture. Stage 3 sets that added tonic input to zero.

Evidence: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\params.py:337-358; build_network.py:456-470 and 1076-1085; build_network_w2lvar.py:350-361; build_network_syn6.py:326-355 and 445-462; runner.py:1543-1568; hunt_training\physical_eval.py:157 onward; hunt_training\controller.py:200. Saved last verified Hunt config: D:\Users\Ben Bolen\NeuroAcceleration_20261104\body\compact_final.json. Historical human wean checkpoint reports_20261016\wean_human\round8_best3.json contains nonzero left/right alpha/gamma values, so not every dissertation-era trial was a single shared amplitude. Those legacy names still do not have the paper's alpha/gamma meanings.

Design recommendation: retain archived key meanings and historical results. Use clearly named new command channels such as RG_F_drive_L/R and RG_E_drive_L/R with an explicit mapping to the cited architecture. Fit shared/mirrored connection parameters; allow independently specified left/right commands, initialized equally for symmetric tasks. Additional Rybak interneuron targets and sensory-gain modulation require explicit pathway contracts, not just renaming four variables. Keep descending posture commands separate from locomotor drives and route them through the requested postural premotor branch.

Weaning recommendation: distinguish artificial startup/bootstrap excitation, ongoing physiological locomotor commands, descending posture, and mechanical harness support. Reduce artificial tonic/bootstrap assistance where it is no longer needed. A healthy locomotion model can retain descending commands; zero supraspinal drive is a specific feedback-dependence/spinal-like experiment, not a universal success criterion. If testing that condition, preserve a defined target cadence/loading task, inspect the operating regime and feedback mapping, and measure all remaining assistance. Harness weaning remains a separate measured-force experiment. No drive/topology changes or training were performed in this comparison.
