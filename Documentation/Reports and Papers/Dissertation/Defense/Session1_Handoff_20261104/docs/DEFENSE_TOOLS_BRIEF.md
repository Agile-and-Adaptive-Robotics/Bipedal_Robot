# Defense brief: tools for neuromechanical experiments

Session 1 • easteregg2 • 3 November 2026. Prepared for the defense and GUI sessions; no manuscript or slide edits.

## Thirty-second contribution statement

I developed a toolchain for building and testing neuromechanical models across graphical and scripted environments. It connects OpenSim- and CAD-derived mechanics to synthetic neural controllers and selectable actuator models. The work includes conversion repairs, editable neural blocks, command and sensory interfaces, accelerated calculations, and repeatable training and verification procedures. The models themselves are reusable experimental tools: a researcher can change a circuit, actuator, or mechanical condition and measure the resulting neural and mechanical response. The contribution is the implemented, tested integration and experimental workflow built on established simulation packages.

## What counts as a tool

A reusable model, library, adapter, or experiment procedure is a research tool when its inputs, assumptions, outputs, and verification are explicit. A GIF demonstrates one experiment produced by that tool. “Developed and verified” is the supported claim here; worldwide priority requires a separate literature comparison.

| Tool | Project-developed capability | Example evidence | Main implementation |
|---|---|---|---|
| OpenSim conversion and repair | Prepare bilateral routes, correct executable-model faults, and check forces/actions | Converted Gait2392 and source comparisons | [repair_robotbody_muscles.py](<D:/GitHub/Bipedal_Robot/Solid_Models/OpenSim/Gait2392_Robotbody/repair_robotbody_muscles.py>); [runner.py: patch_xml/apply_start_pose](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/runner.py>) |
| CAD model and route assembly | Retain frames, inertias, knee loops, and moving/conditional muscle paths | CAD closed-loop and virtual-work checks | [build_ah001.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/build_ah001.py>); [DynamicRoutes](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/dynamic_routes.py>) |
| Graphical SNS library | Edit neurons, synapses, sensory encoders, and actuator blocks | Matched-input neural and actuator tests | [sns_build_library.m](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/sns_build_library.m>); [sns_build_actuators.m](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/sns_build_actuators.m>) |
| Neural/mechanical interfaces | Exchange actuator commands and sensory measurements | Bridge activation, rate, and sensor-readback tests | [run_bridge_tests.m](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/mujoco_bridge/matlab/run_bridge_tests.m>); [sns_build_from_json.m](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/sns_build_from_json.m>) |
| Selectable actuation | Connect characterized BPA forces and opt-in biological tendon models | Force-law parity and separate mechanics checks | [BPAMuscle](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/bpa_muscle.py>); [DynamicThelenBank](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/physics_fidelity/thelen_dynamic.py>) |
| Accelerated execution | Compile repeated operations, preallocate storage, and preserve numerical behavior | Paired fast/reference traces | [ExactFastForward](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/w2l_cpg/exact_fast_eval.py>) |
| Design, fitting, and curriculum | Apply functional targets, mirrored bounds, and bounded experiment requests | Recorded Rustuna studies and stage-specific checks | [hunt_trial_adapter.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/hunt_trial_adapter.py>); [physical_eval.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py>) |
| Verification and communication | Align signals, measure assistance/clearance, repeat trials, and export figures | Fresh-process, perturbation, and saved-state replay records | [AlignedCapture](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/verification_capture.py>); [verification_figures.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/verification_figures.py>) |

OpenSim, MyoConverter, MuJoCo, SNS-Toolbox, Simulink, Simscape, and the MathWorks bridge are third-party foundations. Published neural architectures remain credited to their authors. The project contributions are the implementations, adaptations, blocks, actuator integrations, repairs, and evaluation procedures. A current GUI can make those tools easier to configure and reuse; it does not need to invent another simulator.

## Performance and interface scope

The paired Gait benchmark is **199.624 s native versus 29.284 s accelerated: 6.82× for the full process**. Time within neural calls is **179.079 versus 9.780 s: 18.31×**. Recorded numeric trajectories match exactly. A roughly 200× kernel result must retain its own benchmark scope. Neither number implies that every complete walker, rendering task, or compliant-tendon configuration runs that much faster. [Paired evidence](<D:/Users/Ben Bolen/NeuroCampaign_20261102_9h/figures/BACKEND_AND_VALIDATOR_COMPARISON.md>)

Graphical editing supports circuit inspection and block-level changes. Headless execution runs without opening the graphical interface, making automated trials and optimization practical. An API exposes those commands; the GUI can invoke the same backend and display its records.

Fresh MyoConverter reconversion is not yet a universal button: this machine's clone lacks several historically documented patches, and the conversion driver contains a laptop geometry path. Existing repaired models and their recorded workflows are demonstrated. Likewise, the MuJoCo–Simulink bridge exchanges runtime signals; it does not establish lossless translation of arbitrary Simscape models.

BPA pressure-response/loading-direction variants and compliant biological tendons are separate options. The assumed loading-direction perturbation is not an identified hysteresis loop. The selected CAD262 stepping example uses reference tensions; the historical pressure-pair example uses two pressure-modeled BPA routes. State the exact actuator implementation for each result.

## Committee questions and direct answers

| Member | Likely question | Answer to develop |
|---|---|---|
| David Turcic | What is controlled, how were costs chosen, and what establishes stability? | Show the actual controller, force/support accounting, and objective. Safety and real unloading/clearance precede phase, assistance, and effort in the selected fit. Bilateral homologues share fitted parameters. Local linear checks and finite-time perturbation tests have declared scopes; they are not a whole-body stability theorem. |
| Radhika Reddy | What is neural or biomimetic, and what is an engineering addition? | Identify RG/PF/MN organization, implemented cell dynamics, feedback routes, and actuator conversion. Separate those from posture supervision and calibration targets. The replayed spiking W2L contains 704 spiking cells plus 80 graded cells; do not call it all-spiking or attribute differences between unmatched conditions to neuron type. |
| Bin Jiang | Are the fits unique, and does acceleration change the mathematics? | Explain functional transformations, units, the discrete update, bounds, and optimization terms. Nonconvex fits and synergy factorizations need not be unique. Numerical equivalence, fresh repeats, and local parameter tests answer different questions. Present the paired fast/reference evidence rather than a speed claim alone. |
| Alexander J. Hunt | What did the training workflow add? | Distinguish reused intrinsic/timing calibration from newly fitted feedback and motor settings. Show the same saved controller across named stage checks, the measured apparatus forces, and the failed support rungs. The executable procedures and diagnostics are reusable tools. |

Detailed equations: [controller mathematics](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/DEFENSE_CONTROLLER_MATH_20261102.md>). Chapter 6 already introduces this tools framing in Table `tab:sns_capability`, the construction/interface figure, and the library, bridge, and execution subsections of [60-sns.tex](<D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Dissertation/ProofFinal/chapters/60-sns.tex>). This observation is not a certification of the latest submitted PDF. Full script details and readiness findings are in the [evidence inventory](<D:/Users/Ben Bolen/ToolchainDefense_20261102/inventory/TOOLCHAIN_EVIDENCE_INVENTORY.md>).

## Workflow graphics and current demonstration files

- [Editable vector workflow](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/RESEARCH_TOOLCHAIN_WORKFLOW.svg>), [slide PNG](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/RESEARCH_TOOLCHAIN_WORKFLOW.png>), and [PDF](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/RESEARCH_TOOLCHAIN_WORKFLOW.pdf>).
- [Paired Gait acceleration graphic](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/PAIRED_GAIT_SPEEDUP.png>).
- [Corrected CAD anatomical-foot model and80% support GIF](<D:/Users/Ben Bolen/CAD3DOF_20261102/anatomical_revision/RESULTS_AND_FILES.md>).
- [Other three models, measured support, full paths, and exact repeat records](<D:/Users/Ben Bolen/ToolchainDefense_20261102/other_models/RESULTS_AND_FILES.md>).

The corrected CAD import provides a concrete tool example: topology inspection recovered actual toe hinges lost to welding, retained the neutral assembly and mass, and checked route signs/virtual work before dynamic experiments. Its nominal80% stepping result passes fresh/reference repeats, but only7/8 local perturbations pass.


## Latest CAD contact-handover condition

[Current full paths and records](<D:/Users/Ben Bolen/CAD3DOF_20261102/contact_handover_20261103/RESULTS_AND_FILES.md>): 77.5% actual vertical support during 16–28 s, with 7R/6L mesh-clear unloaded intervals and at least 5.10 degrees knee flexion. Nominal fresh/reference repeats match exactly; 1 ms check passes. Five of eight local ±1% cases pass. The GUI must display nominal success and local robustness separately. The _damped label denotes swing-velocity feedback, not changed physical joint damping.


## Source fidelity, calibration, robustness, and deployment workflow

[Detailed repair/acceleration/training/robustness and Orin handoff](<D:/Users/Ben Bolen/ToolchainDefense_20261102/GAIT_REPAIR_TRAINING_AND_DEPLOYMENT_20261103.md>). [Editable workflow SVG](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/FSA_HUNT_CURRICULUM_WORKFLOW.svg>), [slide PNG](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/FSA_HUNT_CURRICULUM_WORKFLOW.png>), and [PDF](<D:/Users/Ben Bolen/ToolchainDefense_20261102/figures/FSA_HUNT_CURRICULUM_WORKFLOW.pdf>). This is a component-grounded plan and defense walkthrough, not a newly completed end-to-end GUI or a new simulation campaign.


4 November clarification: see CLARIFICATIONS_AND_IMPLEMENTATION_20261104.md for Orin coverage, exact parameter-test scope, requested postural premotor architecture, noise gaps, and same-family scaled-model compatibility. Measured new acceleration: D:\Users\Ben Bolen\NeuroAcceleration_20261104\EXPERIMENT_RECORD.md.
