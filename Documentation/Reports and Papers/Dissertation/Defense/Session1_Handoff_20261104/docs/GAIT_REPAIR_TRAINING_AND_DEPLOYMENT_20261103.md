# Gait2392 source fidelity, fast execution, and the training workflow

Session 1, easteregg2, 3 November 2026. Read-only implementation review and proposed defense/GUI workflow. CAD simulations are stopped pending Ben's routing revision. No converter, production controller, physical muscle/damping law, optimization, or hardware settings were changed in this review.

## The proposed defense demonstration

Select a subject-scaled Gait2392 source, convert and verify its mechanics and muscle functions, bind a named controller, select an explicitly identified actuator/physics profile, derive the applicable functional subnetworks, calibrate the controller with Hunt's staged method, and validate it through the five physical environments. Show a held-out gait/disturbance and aligned neural/mechanical outputs. Export the same controller interface for hardware testing after timing and actuator checks.

The current code supplies many of these components. A single guided application that executes, validates, and versions every branch is proposed orchestration. A GUI should report which steps ran, the actual model/profile identity, and which gates remain unresolved. It must not turn the existence of a helper script into a passed experiment.

## 1. Repair the conversion against the scaled source

There are different repair layers, and their roles should be explicit:

- D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\repair_robotbody_muscles.py constructs the special robot-route OpenSim copy. It mirrors all 46 right-side GeometryPaths and re-enables source appliesForce=false muscles. It is not a general human-muscle calibration program. Its current script also performs work at module scope; a reusable version needs a guarded main and explicit source/output/symmetry options.
- D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\convert_to_mujoco.py runs MyoConverter. Its current geometry folder is a laptop path, and validation=False/speedy=True select a quick first pass. A new CLI should require source/geometry/output paths, expose conversion stages, and use full validation for a fidelity conversion. Quick sampling remains a separately named preview option.
- D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\runner.py, function patch_xml, applies legacy runtime repairs and apparatus settings. It mixes axis adaptation, numerical regularization, muscle pruning, trunk force correction, and support. These belong to separate declared profiles. Body support or biological muscle deletion must not silently enter a source-fidelity conversion.

The source of truth must be the exact **scaled .osim file**, not the generic unscaled Gait2392 model or a previously cached source. Record source, geometry, converter, coordinate adapter, and actuator-law hashes. Cache keys must also include pose/velocity grids, source parameters, unit conventions, equilibrium/activation settings, and schema/version. Source extraction must equilibrate muscles when appropriate, distinguish equilibrium from a transient, and preserve failures rather than inserting zeros.

### Repair gates in order

Before muscle fitting, audit anatomical mass, COM, inertia, the generalized mass matrix/gravity, and the source-to-native coordinate/frame transformation. Legacy follower masses/armatures must be declared separately. A future kinematic-path or reduced-coordinate implementation can avoid fictitious inertial loading only if its path derivatives, generalized forces, and constrained dynamics are independently verified. Do not silently remove numerical mass from an old result or merely change its body-weight denominator.

1. **Kinematics and routing:** match body poses, coordinate signs, muscle-tendon lengths, moving/conditional path points, wraps, and complete constrained moment arms. Check r=-dLmt/dq using finite differences through equality followers. Use full-ROM/multijoint samples plus actual gait trajectories; do not interpret near-zero moment-arm relative errors as large mechanical errors.
2. **Source parameters:** preserve names, actuator type, Fmax, optimal fibre length, tendon slack length, pennation, activation constants, passive parameters, force-velocity parameters, and disabled-force flags. Specify how geometry and strength were scaled. Standard geometric scaling is not a justified universal strength scaling rule.
3. **Force components:** compare passive force and active increments at multiple activations, including zero/low/intermediate/full. For compliant tendons, solve each activation's equilibrium; do not assume force is globally linear in activation. Preserve Fmax as a source parameter; distinguish it from a surrogate force scale fitted to another force law.
4. **Dynamic response:** compare shortening, lengthening, activation/deactivation transients, and, for compliant profiles, fibre/tendon states and equilibrium errors. Check positive and negative velocity branches and discontinuities. Static agreement alone cannot certify these.
5. **Generalized torques:** compare tau=R(q)F and passive joint effects, not only per-muscle peaks. Hold out poses and complete trajectories; report failures by muscle, joint, speed, activation, and loading direction.

MyoConverter's intended sequence is XML/kinematics conversion, moment-arm optimization, then force-length optimization, with separate validation at every step. The project should add the activation/component/dynamic/source-cache contracts above. Official pipeline: https://myoconverter.readthedocs.io/en/latest/pipeline.html.

For a fast native-MuJoCo surrogate, use a sequential, normalized fit: first geometry/length/moment arms, then active/passive force components, then dynamic validation. A candidate objective combines squared length, moment-arm, active-force, passive-force, and torque errors, each normalized to a declared physical scale, plus bounded deviation from source parameters. Keep moment-arm errors in absolute units near zero. Tie homologous parameters when the intended source is symmetric; preserve an explicitly asymmetric subject when it is not. Bilateral source fitting is preferable to arbitrary right-to-left copying. New caches and profiles must be applied before standing solutions, Fmax normalizations, and sensory encoders are built.

### Existing evidence to build on

D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\MOMENT_FORCE_AUDIT.md documents an important failure mode. At neutral, legacy BFSH total force at full activation was similar right/left (736.19/732.60 N), while passive force was 215.57/127.07 N. Matching total peak tension hid a large low-drive difference. The joint bilateral fit reduced held-out BFSH RMSE from 10.37 to 3.07 N and soleus from 81.79 to 29.27 N. These are static, limited-pose six-family results, not all-muscle/dynamic certification.

Relevant implementations are in the same physics_fidelity directory: bilateral_moment_arm_audit.py, shared_force_surrogate_audit.py, export_shared_force_law_v3.py, conversion_cache_static_check.py, conversion_activation_check.py, opensim_static_check.py, and validate_thelen_dynamic.py. Read their manifests before selecting a version; these files do not constitute a universal conversion button.

### Selecting fidelity and damping

Native MuJoCo 2.3.7 muscle force is a normalized FLV surrogate with a constant tendon-length offset. Its actuator lengthrange plus normalized range imply fibre optimum and tendon offset; equal Fmax alone does not make it equal to OpenSim. See https://mujoco.readthedocs.io/en/2.3.7/modeling.html#muscle-actuators.

A source Gait2392 Thelen model should first be compared to the same Thelen law. Selecting Millard is a declared alternative model, not an automatic exact repair. Existing thelen_dynamic.py and series_elastic.py are opt-in components with their own validation scopes. Compliant biological tendon, muscle/fibre damping, joint losses, BPA compliance, and pneumatic hysteresis are independent choices.

A useful passive joint-loss option is tau_d=-b1*v-b2*abs(v)*v-tau_c*tanh(v/epsilon), with nonnegative coefficients and epsilon>0. It dissipates energy because tau_d*v<=0. That property does not identify its physical coefficients. Preserve imported losses and choose additional values from measurements or declared sensitivity studies; do not tune losses solely to obtain a desired walking GIF. Nonlinear losses must participate consistently in the solver's implicit/Jacobian treatment or receive a separate timestep check.

## 2. Acceleration beyond the current backend

Correction (4 November): the current Hunt path already uses sparse_pairwise_eval.py and preallocated I/O with Numba edge updates; exact_fast_eval.py is an older reference layer. The additional compact implementation and measured incremental gains are recorded in D:\Users\Ben Bolen\NeuroAcceleration_20261104\EXPERIMENT_RECORD.md. The proposals below describe remaining/general pipeline work, not features absent from every existing runner.

The next implementation should expose two backends:

- **Reference-compatible exact mode:** retain the validated reduction strategy while reusing voltage/gate/current/output buffers; precompute static inverses/indices; replace eligible allocation-heavy expressions with out buffers or verified compiled kernels. Check output aliasing and keep old/new state separation.
- **Sparse compiled mode:** fixed CSR/edge tables, per-neuron current accumulators, active ionic-channel tables, sparse input/output bindings, and fused C++ or Numba stepping. Storage and synaptic work become proportional to neurons plus actual edges rather than dense N squared. Different summation order, fused multiply-add, or fast-math flags require separate numerical/closed-loop acceptance, not an assertion of bit identity.

Remove permanently absent channels and structurally zero weights at build time. Do not prune a pathway because its current activation is zero, a conditional muscle point is inactive at one pose, or no spike occurred this step. Spiking conductance tails and ionic/adaptation states can still affect future activity. Rebuild/version topology when enabled pathways change. Preserve independent output-observed states and the supported electrical/event laws.

Preallocate encoder, muscle-bank, Jacobian, solver, and diagnostic workspaces too. Resolve body/joint/actuator names once. Reuse geometric calculations within a timestep. Keep numeric factorizations that depend on q current; cache their structure/workspace rather than stale values. Use a warm-started, bounded, fixed-workspace recruitment solver where suitable. Store routine telemetry in bounded buffers and flush it outside the control loop. Profile each component and end-to-end latency before deciding whether C++ work or another optimization is worth doing. MuJoCo's physics is already native C.

Existing paired Gait evidence is D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\figures\BACKEND_AND_VALIDATOR_COMPARISON.md: 6.82x full-process and 18.31x neural calls for its declared comparison. Additional speed gains are proposed until measured. Kernel-only figures must remain distinct from complete simulator and real-time control timings.

## 3. Migration to Jetson Orin and physical actuators

These changes help deployment by making memory, operations, interfaces, and timing predictable. Export a fixed controller contract: initialize, reset, step(sensor_packet, dt, output_packet), plus complete state save/restore and configuration hashes. In simulation, the sensor packet comes from MuJoCo; on hardware it comes from encoders, pressure/force sensors, IMU, and contacts. The actuator backend changes from simulated excitation/pressure dynamics to measured pressure/valve control, retaining explicit units and delays.

D:\GitHub\Bipedal_Robot\Code\Hardware\walker_io\README.md labels the existing Orin I/O code as a **skeleton**. It has sensor/valve/watchdog interfaces and a compiled-SNS hook; it is not completed hardware bring-up or demonstrated 1 kHz performance. The real robot does not necessarily need a full musculoskeletal simulation each control tick.

For one circuit of hundreds to a few thousand neurons, a compiled ARM CPU loop is the first deployment baseline I would benchmark. A GPU may help large populations, batched rollouts, optimization, perception, or a sufficiently large estimator. Kernel launches, synchronization, precision, and coupled serial work can outweigh GPU benefits for a small single network. Integrated Jetson memory removes some copying needs, not synchronization or scheduling costs. NVIDIA recommends profiling and minimizing unnecessary transfers: https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/index.html.

Keep physical time constants separate from computation speed. For example, a 1 kHz sensor/valve task can hold commands between 500 Hz neural updates; migrating the I/O rate does not automatically justify halving the neural timestep or retiming the circuit. Validate any new multi-rate schedule.

Measure complete sensor-to-command latency, maximum and high-percentile jitter, missed deadlines, startup/JIT costs, thermal behavior, and failure handling. Faster-than-real-time simulation alone does not prove a hardware deadline. Native ARM builds, configured real-time scheduling, fixed buffers, and validated watchdog behavior are separate work. NVIDIA's supported RT-kernel options are release-specific: https://docs.nvidia.com/jetson/archives/r36.4.4/DeveloperGuide/SD/Kernel/RealTimeKernel.html. GPU physics would be a separately validated backend; do not assume the current MuJoCo 2.3.7 callbacks and features automatically migrate.

## 4. FSA, Hunt calibration, and the five physical stages

These are three complementary levels:

1. **Functional design:** specify represented quantities, units, gains, saturation, temporal response, and pathway sign. For compatible GLIF neurons, Szczecinski 2020 Table 2 selects the encoding domain, transient type and duration, bias, membrane timing, synaptic decay, and functional conductance. It is an analytical design sequence, not seven physical training stages. Graded circuits use the appropriate graded FSA; NaP/Ca ionic oscillators need their own dynamics. Source: https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2020.577804/full and graded FSA https://www.frontiersin.org/journals/neurorobotics/articles/10.3389/fnbot.2017.00037/full.
2. **Hunt A-D calibration:** A measures/selects an intrinsically slow free rhythm, approximately twice target period without drive; B calibrates expected-afferent timing; C freezes timing and fits feedforward motor peaks; D freezes feedforward and refines sensory-driven waveforms. Cite A. J. Hunt, Real-time control of dynamic walking in a dog-like robot, INRC2 Meeting, 15 February 2019. Ia/Ib refinement is the lecture step; adding II pathways is an explicit project extension where applicable.
3. **Five curriculum environments:** calibration plus recurring closed-loop acceptance gates. An unsuccessful stage is not repaired by merely renaming it or changing all network weights simultaneously.

| Physical stage | What it tests | How calibration complements it |
|---|---|---|
| 1. Deafferented suspended motion | Calibrate intrinsic rhythm; verify suspended motion | Includes Hunt A and relevant feedforward calibration. Current physical_eval only evaluates saved fits; intrinsic.py is a separate fitter. |
| 2. Afferented suspended motion | Calibrate afferent timing, feedforward, and sensory waveforms in ordered blocks | Hunt B/C/D fitting complements closed-loop air motion. Feedforward calibration may also precede stage 1. Current fitters are separate from physical_eval. |
| 3. Ground standing | Load-bearing posture and balance feedback | Calibrate a separate postural/descending module and verify muscle authority. This is additional control design, not a fifth Hunt step. |
| 4. Harness-supported ground stepping | Contact transitions, stance load, unloading, and real mesh clearance | Reuse the saved A-D solution; make bounded, documented changes to implicated blocks and recheck preceding gates. |
| 5. Ground contact refinement and weaning | Contact-PF integration, decreasing assistance, broader gait/disturbance acceptance | Freeze candidate blocks for tests, and revise only identified functions if needed. Held-out robustness tests also apply at earlier stages. |

There is no strict one-to-one correspondence. FSA establishes component contracts; Hunt supplies an ordered calibration loop; the five stages progressively expose the controller to mechanics. A body/actuator change invalidates relevant component and physical checks, so the workflow loops back to those checks. One result must not inherit a different controller's pass.

Exact local implementations: D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\table2_glif.py, graded_profile.py, design_gate.py, and DEFENSE_CONTROLLER_MATH_20261102.md; hunt_training\intrinsic.py, sensory_timing.py, motor_fit.py, train_fourcell_hunt.py, and physical_eval.py. Existing electrical/model-specific limitations remain in ALL_CONTROLLER_DISPOSITION_20261102.md. The presence of the seven-step module is not evidence that every current spiking family completed it.

## 5. What the IK/gait database should do

The reference bank is D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\gait_refs. It supplies candidate cadences, joint trajectories, phase/duty targets, and design bandwidth contexts. The earlier audit used a shortest usable cached stride of0.559 s and a fourth-harmonic7.1556 Hz context for one encoder calculation. These are cached design contexts, not universal physiological limits.

The active Hunt path is narrower: hunt_training\prepare_targets.py builds model-consistent targets from its selected IK/GRF pair; physical_eval.py loads the resulting targets.npz. kine_ref.py also names a specific subject01 IK/GRF source. The code supports reference comparisons, but the latest nominal success and small parameter checks do not demonstrate an exhaustive held-out test of the entire library.

Use the bank in four places:

- **Design/scaling:** select the intended subject/speed domain and derive sensor ranges, fibre lengths/velocities, cadence, and signal bandwidth from it. Retarget anatomy/coordinates explicitly. Missing MTP, transverse, or lumbar channels remain missing rather than becoming zero-valued training evidence.
- **Training targets:** combine q, filtered qdot/qddot, geometry, inverse dynamics, and GRF/CoP where available to derive feasible muscle/afferent targets. IK alone does not uniquely determine activation, force, or EMG.
- **Held-out transfer:** separate participants and complete trials, not random adjacent frames or cycles from one recording. Include slower/faster walking and, only when in scope, running. Recompute targets for each mechanics/actuator profile.
- **Scoring:** phase-normalize shape while scoring cadence, phase, offsets, amplitude, duty, contact loading, support forces, clearance, separation, and slip separately. Otherwise time normalization can hide incorrect speed or dragging.

## 6. Perturbation and neural noise

Recent CAD local tests sample static parameter factors once per trial; the Gait eight-case result also has a limited selected-scalar scope. They are not ongoing stochastic noise, physical push tests, a full gait-library generalization study, or a state-space region of attraction. The inspected W2L spiking launcher similarly exposes a parameter-perturbation seed; its sensory populations use deterministic initial phase staggering.

Proposed robustness layers, each with recorded seeds, units, timing, and independent held-out tests:

| Layer | Examples | Required distinction |
|---|---|---|
| State | Initial joint/COM velocity, posture, neural state, oscillator phase | Different from changing a fixed conductance |
| Parameters | Muscle strength, delays, gains, source uncertainty | Fit homologues together; declare independently asymmetric tests separately |
| Sensors and timing | Force/position noise, bias, delay, dropped packets, sample jitter | Keep causal filtering and actual physical units |
| Neural | Current fluctuations, spike reset/threshold variability, synaptic variation | Record exact injection location, correlation, and time scale |
| Mechanical/environment | Phase-specific pelvis impulses, surface height/friction, belt-speed changes, support changes | Log applied wrench and test actual recovery and slip |
| Actuators | Pressure delay, valve deadband, force-map/loading-direction errors | Keep biological and BPA uncertainty models distinct |
| Reference/task | Held-out participants, speeds, amplitudes, cadence changes | New reference domains are not synthetic Gaussian noise |

Tiny neural noise is not universally beneficial. It can desynchronize rate-coded populations or reduce artificial phase locking, and it can also corrupt timing or destabilize a circuit. The 2020 FSA paper's neuromechanical comparison explicitly uses reset noise and compares individual cells with populations; that is a sound model-specific motivation for a matched test. It does not establish that all controllers improve when random current is added.

Start with clean runs plus several declared noise severities and independent seeds. Define severity relative to an encoded physical range, not the instantaneous signal (which can be zero). Compare independent, common-mode, and correlated disturbances separately. A useful colored-current model is eta_next=a*eta+sigma*sqrt(1-a^2)*xi, a=exp(-dt/tau), so stationary standard deviation and correlation time remain explicit when dt changes. Do not use arbitrary per-step Gaussian amplitudes and then assume equivalent noise after changing timestep. Reset-noise models need their own rate/bias calibration.

Measure task performance and recovery distributions rather than only mean score. Keep deterministic source/fast equivalence fixtures; use common recorded noise sequences or per-cell seeded streams for comparisons. Training with bounded disturbances can regularize a solution, but withheld disturbances and clean conditions must still pass. Noise should not be used to conceal a structural wiring or support-control failure.

## Suggested 3-minute defense walkthrough

1. **Select:** scaled source model, geometry, and a reference trial; show source/frame/unit identity.
2. **Convert and verify:** show one muscle's OpenSim/MuJoCo moment-arm and active/passive-force overlays plus the acceptance ledger. Mention that raw conversion and scientific validation are different steps.
3. **Bind the controller:** match names, side, function, and units, not just actuator count. Show one functional subnetwork and a named actuator profile. Source strength, neural encoding, and delays are recalibrated for the new mechanics.
4. **Design and calibrate:** show the FSA/Hunt/physical-curriculum graphic. Point to one analytical gain/time-constant example and the ordered fitting blocks.
5. **Run and challenge:** launch a short prepared headless trial, display a few aligned RG/PF/MN/joint/load cycles, and select a held-out gait or disturbance. Use saved validated results for the expensive fitting parts rather than promise a full live optimization.
6. **Export:** show the model/controller/configuration manifest, GUI-editable network, trace/GIF outputs, and the proposed compiled sensor-to-actuator interface for Orin.

Same-topology, same-muscle Gait2392 models are structurally compatible and can reuse controller bindings. Matching counts alone between unrelated models is insufficient: names/order, muscle actions, normalization, strength, mechanical DoF, sensor directions, and actuator states must agree. New tendon or damping choices require relevant revalidation before controller training. The demonstration should use a supported branch that actually exists and clearly mark orchestration or hardware components that are still proposed.
