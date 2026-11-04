# Evidence inventory: neuromechanical research tools

Machine: easteregg2. Session 1 inventory audit, 2026-11-02. Read-only source/report inspection; no model or manuscript changes and no simulations. This supports the defense framing and the GUI session's workflow design.

## A defensible contribution statement

The contribution is a configurable neuromechanical experimentation toolchain: it brings OpenSim- and CAD-derived mechanics into executable simulations, connects synthetic neural controllers and alternative actuator laws, and supplies accelerated, repeatable training and verification workflows. The reusable models are experimental instruments within that toolchain. A GIF is a result produced by a tool, rather than the tool itself.

Use “developed and verified” or “project-developed tools.” This inventory establishes implementation and bounded demonstrations. It does not establish worldwide novelty or priority over all published systems. OpenSim, MyoConverter, MuJoCo, SNS-Toolbox, Simulink, Simscape, and the MathWorks bridge remain credited third-party foundations. Deng, Li, Szczecinski, Rybak, and other source neural architectures remain credited scientific designs; implementing or adapting them is not inventing their neuron equations.

Current dissertation source already gives this distinction: D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal\chapters\60-sns.tex, lines48–73 and165–227. Table tab:sns_capability groups models/tools. Line73 identifies model implementations, custom blocks, actuator integrations, and evaluation procedures as contributions. This inventory does not certify correspondence to the latest submitted PDF or authorize manuscript edits.

## Implemented tools and script/function map

### 1. OpenSim model preparation and MyoConverter repair workflow

Input: human gait2392_simbody.osim or Ben's gait2392_robotbody.osim. Output: generated symmetric robot .osim, converted MJCF, repaired executable model, and force/action audits.

- D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\repair_robotbody_muscles.py: mirrors right GeometryPath to left and re-enables disabled muscles, preserving source. Handles static, conditional, and moving path points. This creates gait2392_robot.osim; it is not a repair of every human Gait conversion defect.
- D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\convert_to_mujoco.py: conversion driver for human/robot sources; strip_missing_meshes() handles decorative mesh omissions. MyoConverter performs the actual conversion.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\runner.py: patch_xml(), apply_harness(), apply_start_pose(), _project_followers(), and solve_standing_activations(). Bounded conversion repairs include force-scale fixes, prescribed ground contact set, moving-point follower treatment, pose initialization, solver/time-step settings, and explicit restraint mechanics. Some historical options are imperfect; audit the selected captured model, not a generic README default.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\MOMENT_FORCE_AUDIT.md, STANDING_INITIALIZER_AUDIT.md, and CONTACT_INITIALIZER_AND_ALIGNMENT.md: later evidence of force-fit asymmetry, source-curve comparisons, passive-force/sign corrections, equality projection, and ground alignment.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\README.md: records local MyoConverter patches in optimization/utils/UtilsOpensim.py, optimization/model_states/OsimMuscleStates.py, conversion_steps/O2MStep3.py, and optimization/utils/UtilsForceOpt.py. They catch failed equilibrium poses, retain usable masks, skip unusable force fits, and fix Windows worker serialization. These patches are in an external MyoConverter checkout, not all bundled here. Placeholder force fits must stay labelled; conversion completion is not dynamic equivalence to OpenSim.

Status: implemented and demonstrated for named converted models. The repair workflow and diagnostics are tools; MyoConverter itself is not a new project tool.


### Fresh reconversion readiness audit (important)

The historical README patch descriptions are not proof that the currently available external clone carries those patches. On easteregg2, D:\GitHub\myoconverter exists, but source inspection on2026-11-02 finds:

- D:\GitHub\myoconverter\myoconverter\optimization\utils\UtilsOpensim.py:963 directly calls equilibrateMuscles in getMuscleForceMaps without the described per-pose catch.
- D:\GitHub\myoconverter\myoconverter\conversion_steps\O2MStep3.py:325 sets res_opt=None in a fallback, while line336 still calls len(res_opt).
- D:\GitHub\myoconverter\myoconverter\optimization\utils\UtilsForceOpt.py:144–146 creates a shared MjModel and passes it to a multiprocessing initializer, rather than the described path-only Windows worker fix.
- D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\convert_to_mujoco.py:26 hardcodes the laptop geometry path C:\Users\Ben\Documents\GitHub\myoconverter\models\osim\Gait2354Simbody\Geometry. Current defaults include speedy=True and validation=False.

Therefore the demonstrated tool is the recorded conversion/repair workflow and existing repaired model implementations, not a currently universal one-click fresh reconversion service. The GUI should expose preflight/dependency status, make geometry/environment paths configurable, and package/check the required source patches before enabling this branch as turnkey. No source code was changed in this inventory audit.
### 2. CAD-derived neuromechanical model assembly and route transfer

- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\build_ah001.py: CAD mass/geometry and audited joint-frame assembly, including both closed-loop four-bar knees; cad_to_isb() and inertia_for_group(). MJCF is authoritative for loop dynamics; companion URDF cuts loop edges.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\add_torso.py: adjustable lumbar/spine/payload construction.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\inspect_source_routes.py: source route/function and actuator-roster extraction.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\prepare_opensim_routes.py: basis(), alignment(), apply(), distal_shank_center(), main(); maps body landmarks and source route points into CAD body frames.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\build_opensim_routes.py: build(specification,output,include_reference=False), creates routed actuators/sites.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\dynamic_routes.py: spline_value_derivative() and DynamicRoutes.update()/coordinates()/apply_force_correction(). Evaluates moving and conditional sites and includes extra virtual-work derivatives in force and sensory velocity. This is a reusable mechanics adapter, not simply a mesh renderer.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\robot_kinematics.py: knee_pose() and anatomical_knee_tangent(), preserving four-bar constraints and hip authority.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\verify_dynamic_routes.py and verify_knee_tangent.py: numerical derivative/function checks.
- Latest articulated-foot extension is session-scoped, not yet a polished generic importer: D:\Users\Ben Bolen\CAD3DOF_20261102\source\export_step_geometry.py and build_articulated.py; selected outcomes indexed at D:\Users\Ben Bolen\CAD3DOF_20261102\stage4\RESULTS_AND_FILES.md. It imports STEP solids/mass/inertia, normalizes reference frames, and constructs articulated feet; further anatomical joint correction is underway.

Status: verified model-specific CAD route, frame, inertia, loop, and neural integration workflows. Do not label this a general-purpose arbitrary-CAD-to-MuJoCo converter or a lossless Simscape compiler.

### 3. SNS-Simscape graphical neural and actuator library

Actual path is D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape, not a guessed Neuromechanical_Models path.

- sns_build_library.m: generates SNS_Library.slx with neuron, synapse, sensory, activation, and summation blocks. Synapses emit [g,g*E]; neurons sum internally.
- sns_build_actuators.m and sns_test_actuators.m: generate/test BPA and biological muscle blocks. Existing reported matched-force error4.9e−10N is scoped to actuator block tests.
- sns_build_from_json.m: sns_build_from_json(jsonPath,mdlName,slxDst) creates editable Simulink networks from the supported JSON schema.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\export_network_json.py: graded-network exporter for that generator. Its source asserts NonSpikingNeuron and stores selected graded parameters. It is not evidence of general gated/spiking semantic round-trip.
- Demos and runners under D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos: sns_build_demo/sns_run_demo, sns_build_cpg_demo/sns_run_cpg_demo, sns_build_beer_demo/sns_run_beer_demo, sns_build_deng_cpg/sns_run_deng_demo.
- snsfig.m, sns_draw_circuit.m, sns_function_subnetworks.m: graphical/scientific display tools.
- GUI prototype: D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\gui\SNSGui.m. The currently assigned GUI session owns further frontend integration.

Evidence: README_SNS_Simscape.md and reports_20261024\DISSERTATION_EVIDENCE_INVENTORY_20261025.md. Selected neural twins agree with NumPy under matched inputs; reported410-neuron test maximum deviation4.201e−6mV is a particular historical export. Deng's1.938s rhythm is an implementation reproduction, not a newly invented CPG.

### 4. Simulink–MuJoCo command/sensor coupling

- D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge\matlab\laptop_build_bridge.m: build driver for the existing MathWorks blockset.
- make_simbridge_xml.py: creates sensor-enabled MJCF copy, preserving conversion source.
- run_bridge_tests.m: run_bridge_tests(), run_case(), parse_bus(); single-muscle, clock-rate, and sensor checks.
- e1_sns_reflex_loop.m and e2_cpg_open_loop.m: closed-loop reflex and all-channel neural command demonstrations.
- BRIDGE_REPORT.md: source-code patch to upstream initData() resets a valid keyframe; corrected sensor names jointpos/jointvel/actuatorpos/actuatorvel/actuatorfrc. Reported sensor readback was bit-exact over the scored80-time-step window. Later report chronology must be read past the historical PARKED section.

Status: bridge is third-party; local compatibility/init fixes, SNS wiring, and verification are the project integration. Demonstrated on laptop/EB475WS4. The extended run guide says no compiled mex found in this easteregg2 checkout, so do not imply installed everywhere. The bridge lets a Simulink controller operate MuJoCo mechanics; it does not automatically translate every Simscape mechanism or arbitrary Simulink model into standalone MuJoCo/Python.

### 5. Restricted mechanical interchange

- D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\osim_import\mjcf2urdf.py: selected MyoConverter MJCF→URDF geometry/kinematic subset, with anchor-offset chain construction.
- D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\sns_osim_import.m: URDF→Simscape route.

Critical boundary: mjcf2urdf.py explicitly ignores sites, tendons, muscles, and equality constraints; neural/actuator layers must be attached separately. CAD closed loops require deliberate reconstruction. Do not draw a double-ended “lossless import/export” arrow.

### 6. Custom actuator interfaces and selectable force laws

- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\bpa_muscle.py: festo4(), maxBPAforce(), tendon_spring_rate(), BPAMuscle. Port of Ben's characterized force law and series-compliance calculation.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\bpa_mujoco.py: BPAMuscleSystem, set_activation(), attach(), callback force application.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\add_bpa_to_mjcf.py: add_bpa(), tendon/site/zero-gain actuator construction.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\sns_bpa_demo.py: closed-loop neural/BPA/knee integration example.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\actuator_interface.py: force_at_pressure(), inverse_pressure(), pressure_step(), neural_step(); pressure-state and static inverse fixtures.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\series_elastic.py: SeriesElasticBank, step_bank(), solve_one(); early elastic-tendon prototype.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\thelen_dynamic.py: DynamicThelenBank/RigidThelenBank/ConstitutiveBackend, force_velocity(), balance_and_slope(); later source-aligned opt-in force backend.

Status distinction: biological tendon compliance, joint damping, and BPA series compliance are separate mechanisms. The selected current CAD262 stepping GIF uses54reference-tension actuators. The historical pressure-BPA CAD demonstration uses only2pressure-modeled vastus routes and52reference tensions. Its pressure/BPA pair is a real bounded integration demonstration, not a full54-BPA result. The loading-direction term in the neural/BPA fixture is explicitly assumed phenomenological behavior, not an experimentally identified hysteresis loop; pressure fill/vent constants are assumptions. Call it an experimental stateful actuator/uncertainty test interface, and state calibration scope.

### 7. Accelerated execution and persistent optimization workflows

- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\exact_fast_eval.py: ExactFastForward and compiled _edges(). Calculates active edges with preallocated work arrays while preserving the reference NumPy reduction order. Source explicitly rejects spiking/electrical channel types; separate event/spiking implementations require their own tests.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\virtual_walker_20261031\rustuna_runtime.py and rustuna_myo_worker.py: persistent worker interface to the scientific environment; protocol/resume capability checks nearby.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\hunt_trial_adapter.py: preflight(), evaluate(), handle_request(), dispatch to owned experiment outputs.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\design_gate.py: inspect_gate()/require_calculation_release(); validates provenance/design readiness before calculations.
- Spiking interface work: fast_event_port.py, event_conductance_port.py, and spiking_family_evidence.py in that same fsa_design folder. Electrical equivalence is a different gate from body success.

Best full-pipeline benchmark to cite: D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\BACKEND_AND_VALIDATOR_COMPARISON.md. Full source process199.624s native vs29.284s accelerated =6.82×; neural-update portion179.079s vs9.780s =18.31×; saved trajectories match exactly. A separate validator optimization is1.31×. Do not relabel kernel-only~200× as complete-walker throughput. Elastic-tendon roots retain additional overhead, and its isolated92-unit bank benchmark is not whole-walker performance.

### 8. Reusable experiment, acceptance, and figure tools

- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py: experiment wrapper; build(), harness(), standing_solver(); stages and body/controller integration.
- controller.py, intrinsic.py, pf_targets.py, motor_fit.py, train_fourcell_hunt.py, and fourcell_curriculum.py in the same folder: staged design/fitting components. Hunt's training stages are sourced methods; implementing executable fitting/protocol machinery is the tool contribution.
- regional_contacts.py: RegionalContacts; contact_events.py: HeelStrikeEncoder; one event per armed landing with hysteresis/refractory handling rather than sustained reset throughout stance.
- verification_capture.py: AlignedCapture and sparse_current_snapshot(); aligns neural/mechanical/control signals for meaningful comparisons.
- verification_figures.py: read_run(), choose_window(), plot(), render(); produces short-cycle figures, recorded-state animations, and print stills.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\support_mass_accounting.py: selected_robot_mass()/renormalize_saved_record(); prevents prescribed-lift labels from replacing measured assistance.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\how_to_run_code_extended.md and ZCODE_SIMULATION_PROMPTS.md: reusable operator/agent workflow, machine guards, model identity, and owned output conventions.

These are scientific tools when they make tests repeatable and expose failure criteria: protected knees, actual whole-mesh foot clearance, alternating unloading, signed harness forces, bilateral parameter ties, fresh-process repeats, and local perturbation checks. They are not a mathematical global-optimality or universal walking guarantee.

## GUI workflow that can honestly be presented

1. Import a source model: OpenSim via MyoConverter, CAD via its audited export/assembly route, or a prepared model manifest.
2. Inspect and repair: frames/units, anatomical landmarks, joint axes and travel, inertias, loop closure, tendon attachments, moving/conditional points, and single-actuator signs.
3. Select/controller-edit: known SNS family and supported neuron implementation; export exact parameters and readable subnetworks.
4. Select actuators: native biological muscle, explicitly opt-in compliant-tendon backend, reference tension, or a characterized BPA/pressure-state variant. Expose units and calibration scope.
5. Choose an experiment: intrinsic/air/standing/supported-ground/free trials, stimulus layer/time/amplitude, perturbations, objective, and mirrored parameter bounds.
6. Run the headless backend: verified fast calculation mode, Rustuna where supported, bounded/resumable workers, explicit concurrency and output ownership.
7. Validate and inspect: measured support, force balance, clearance/unloading, knee limits, neural activity, fresh repeat, parameter neighborhood, and failed-trial records.
8. Export: figure/GIF, model, parameter/actuator/controller manifest, source hashes, machine/environment identity, and exact replay command.

This frontend would be an orchestration and inspection app around demonstrated but partly model-specific tools. Its contribution can be workflow accessibility and reuse. Label unfinished end-to-end automation as proposed/in development. GUI editing and headless trials complement each other; Python backends can be controlled by a desktop or web GUI, so absence of a universal Python block language is not a fundamental barrier.

Recommended manifest fields: source and generated model IDs; anatomical frames and units; constraints/independent DoF; actuator roster/force backend and calibration; neural family and supported features; timestep/backend; protocol and support settings; mirrored optimization parameters/cost terms; acceptance thresholds; dependencies/machine; source hashes; exact replay; saved results/assets. Preview alone must never count as a passed experiment.

## Strong defense examples, with scoped evidence

1. SNS-Simscape custom blocks and matched-input reference tests demonstrate graphical neural construction and executable parity.
2. Existing MathWorks bridge plus local fixes demonstrates command/sensor coupling under tested rates and initial conditions.
3. Source-derived dynamic muscle-route transfer demonstrates handling of moving/conditional attachment mechanics, with virtual-work derivative checks.
4. Gait2392 paired fast/reference trajectory checks demonstrate acceleration without changing the accepted trajectory under the tested backend.
5. Current Gait biw2lvar-Hunt-PT and CAD262 records demonstrate staged experiments, measured assistance, reproducible outputs, and visible failure boundaries.

Current result indices: D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\RESULTS_AND_FILES.md and D:\Users\Ben Bolen\CAD3DOF_20261102\stage4\RESULTS_AND_FILES.md. Do not conflate their different neural families, actuators, support levels, or model physics.