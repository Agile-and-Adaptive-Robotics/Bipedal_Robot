# GUI workflow handoff: neuromechanical experimentation

Session 1 • easteregg2 • 3 November 2026. For the GUI and defense sessions. This is a proposed orchestration contract around existing tools, not a completed universal launcher. It authorizes no manuscript changes.

## Workflow and backend responsibilities

**Import → inspect/repair → choose controller → choose actuation → configure experiment → execute → validate → replay/export.** Keep source identity, units, and selected physics visible throughout. A saved experiment manifest should drive both the GUI and headless execution.

| Stage | Inputs and outputs | Actual entry points and readiness |
|---|---|---|
| 1. Import | OpenSim `.osim` plus geometry, CAD exports, or a prepared-model manifest → versioned model artifact | [convert_to_mujoco.py](<D:/GitHub/Bipedal_Robot/Solid_Models/OpenSim/Gait2392_Robotbody/convert_to_mujoco.py>); [build_ah001.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/build_ah001.py>). Model-specific branches exist; fresh conversion needs preflight/packaging. |
| 2. Inspect and repair | Model plus anatomical landmarks → frames, coordinate map, attachments, masses, and checked mechanics | [prepare_opensim_routes.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/prepare_opensim_routes.py>), [robot_kinematics.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/robot_kinematics.py>), and [runner.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/runner.py>). Preserve source and generate repaired copies. |
| 3. Controller | Family, readable graph, neuron laws, and checkpoint → executable neural model | [export_network_json.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/export_network_json.py>) and [sns_build_from_json.m](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/sns_build_from_json.m>) support a declared graded schema. Arbitrary gated/spiking round-trip is not established. |
| 4. Actuation | Route roster, dimensions, force law, and calibration → actuator backend and sensor map | [BPAMuscleSystem](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/bpa_mujoco.py>), [actuator_interface.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination/actuator_interface.py>), and [thelen_dynamic.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/physics_fidelity/thelen_dynamic.py>). Expose biological, reference-tension, pressure-BPA, and tendon options separately. |
| 5. Experiment | Saved controller plus stage/stimulus/objective/bounds → immutable run request | [hunt_trial_adapter.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/hunt_trial_adapter.py>): `preflight`, `evaluate`, `handle_request`; [design_gate.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/design_gate.py>) checks prerequisites. |
| 6. Execute | Run request → bounded results and status | [physical_eval.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/physical_eval.py>), [run_virtual_walker.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/virtual_walker_20261031/run_virtual_walker.py>), [rustuna_myo_worker.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/virtual_walker_20261031/rustuna_myo_worker.py>). Dispatch to the compatible runner; do not substitute a generic backend silently. |
| 7. Validate | Aligned neural/mechanical traces → model-specific verdict and diagnostics | [AlignedCapture](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/verification_capture.py>), [support_mass_accounting.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/fsa_design_20261101/support_mass_accounting.py>), and existing clearance/contact validators. Execution success and accepted stepping are separate fields. |
| 8. Replay/export | Saved states and verdict → plot/GIF/still, provenance, and replay command | [verification_figures.py](<D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/hunt_training/verification_figures.py>): `choose_window`, `plot`, `render`. Rendering changes the view, not the simulated movement. |

The Simulink branch uses SNS-Simscape blocks and the existing [MathWorks bridge integration](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/mujoco_bridge/BRIDGE_REPORT.md>) for runtime command/sensor exchange with MuJoCo. Simscape mechanics can also run with the block controller. Draw these as alternative execution branches. [mjcf2urdf.py](<D:/GitHub/Bipedal_Robot/Code/Matlab/SNS_Simscape/osim_import/mjcf2urdf.py>) ignores tendons, muscles, and equality constraints; it does not justify a lossless two-way mechanical-conversion arrow.

## Preflight and result rules

- Check the machine, scientific environment, geometry files, compiled bridge availability, licenses, and source hashes before enabling execution. The available MyoConverter clone needs patch reconciliation; its current conversion driver also has a laptop path and fast/unvalidated defaults.
- Show anatomical versus native axes, anatomical coordinates versus auxiliary follower coordinates, and independent DoF versus raw solver coordinates. Separate anatomical, numerical, and apparatus masses when defining body weight.
- Keep passive joint damping, contact parameters, biological tendon compliance, BPA series compliance, and pressure dynamics as separate selections. Record assumed versus identified parameters, including hysteresis/loading-direction behavior.
- Preserve recorded initialization and actuator ordering. Require explicit mirrored bounds for a bilateral fit and a new owned output folder; never overwrite a selected run.
- Limit heavy simulations to three concurrently, expose cancellation/resume status, and keep compact outputs outside the repository.
- Report actual measured vertical assistance and horizontal/rotational guidance, whole-mesh clearance, unloading, knee safety, repeatability, and parameter-test scope. A viewer opening successfully is not a validation result.

## Minimal manifest contract

The existing registry is [MODEL_REPLAY_MANIFEST.json](<D:/Users/Ben Bolen/ToolchainDefense_20261102/other_models/MODEL_REPLAY_MANIFEST.json>). It records programs, configurations, traces, media, hashes, and support scopes. A proposed GUI run request can reference that registry instead of copying file paths into every control:

```json
{
  "schema": "neuro_experiment_request_v1_proposed",
  "registry": "D:/Users/Ben Bolen/ToolchainDefense_20261102/other_models/MODEL_REPLAY_MANIFEST.json",
  "model_id": "gait2392_hunt_pt_rc_s4",
  "operation": "replay",
  "physics_policy": "preserve_recorded",
  "parameter_policy": "bilateral_ties",
  "output_root": "D:/Users/Ben Bolen/ToolchainDefense_20261102/gui_runs/example_001",
  "overwrite": false,
  "fresh_process": true,
  "acceptance_profile": "model_specific_recorded_gate"
}
```

This example is a proposed wrapper schema, not an existing executable configuration. Its adapter must resolve and validate the registry entry. The resolved manifest also needs environment/version, frame/coordinate map, initialization, exact controller and actuator law, timestep/backend, protocol, objective/bounds when fitting, acceptance thresholds, source hashes, time/force denominators, and replay/media outputs. Store actual measured outcomes separately from requested settings.

## Three fresh examples already available

| Registry ID | Demonstrated condition | GUI display that preserves its identity |
|---|---|---|
| `gait2392_hunt_pt_rc_s4` | Graded Gait2392; 16 s; 48.85% net vertical assistance during 5–13 s | Guided steps in place; 115.56 N horizontal restraint RMS; 92 actuator entries. |
| `w2l108_graded` | Canonical newer W2L; 10 s; 39.24% net vertical assistance | 108 graded cells, 12 uniarticular muscles, retained rails/guides, and treadmill-relative progression. |
| `w2l704_spiking` | Combined central/sensory-spiking W2L; 10 s; 42.91% net vertical assistance | 704 spiking cells plus 80 graded cells; retained pooled-home condition and guides. |

These are accepted separate saved conditions, not a matched neuron-type experiment. The two W2L examples are not the biarticular rat. Their exact program/configuration/media paths and invocation details are in [other-model results](<D:/Users/Ben Bolen/ToolchainDefense_20261102/other_models/RESULTS_AND_FILES.md>). Reuse those entries as initial GUI examples; provide a new destination and tag rather than rerunning the fixed-destination execution-record wrappers.

Further detail: [defense tools brief](<D:/Users/Ben Bolen/ToolchainDefense_20261102/DEFENSE_TOOLS_BRIEF.md>) and [full evidence inventory](<D:/Users/Ben Bolen/ToolchainDefense_20261102/inventory/TOOLCHAIN_EVIDENCE_INVENTORY.md>). The existing GUI session owns frontend implementation; this handoff changes no GUI, manuscript, or Overleaf project.

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
