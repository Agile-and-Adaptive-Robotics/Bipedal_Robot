# Saved follow-up: neural BPA coordination and actuator-model robustness

Saved by Session 1, 30 October 2026, at Ben's request. **Activated by Ben on 31 October 2026.** The first bounded tests and CAD audit are recorded in `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\EXPERIMENT_RECORD.md`; read that record before repeating work. The prompt below remains the broader continuation scope.

## Prompt to give the next session or a physics/controller subagent

Work on **easteregg2**, in **`D:\GitHub\Bipedal_Robot`**. Check hostname before doing anything. Use **`D:\Anaconda\envs\myo\python.exe`**, not bare Python. Read AGENTS.md and applicable Python, MuJoCo, SNS, neuroscience, and BPA modeling guidance before changing or running models. Use Numba, preallocation, and existing exact neural speed-ups wherever appropriate. Do not commit or push.

Implement and test **neural coordination of the CAD robot with a BPA actuator interface**, starting from Ben's static characterization and explicitly testing robustness to imperfect actuator dynamics. A neuron-to-actuator conversion is expected: AnimatLab polygon maps, SNS–Simscape conversion blocks, and muscle excitation/activation interfaces do not invalidate neural coordination. Neural output need not represent a literal force target.

The intended loop is:

**SNS recruitment/timing → actuator command conversion → actual pressure response → BPA force at actual contraction/loading state → robot mechanics → sensory feedback to SNS.**

The conversion may produce pressure references directly, or invert the static characterization to approximate a requested tension. Assess achievable forces at current contraction, pressure limits, and slack conditions. If an inverse is infeasible, saturate explicitly and record the force shortfall. Avoid treating desired pressure as actual chamber pressure when testing pressure dynamics. Shortening, lengthening under tension, and approximately isometric states are determined from actual BPA kinematics, not from pressure direction alone. Do not add biological muscle eccentric force enhancement to a BPA without supporting data.

### Start by auditing the existing implementation

- Robot runner: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\run_neural_robot.py`.
- BPA forward model: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\bpa_muscle.py`. Trace it to the current MATLAB characterization sources; do not assume all historical parameter files are current.
- Historical GIF configuration/source: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\slide_media\session1_20261030\SIMULATION_GIF_EDITOR_HANDOFF.md`.
- Biological compliant-tendon and damping prototype/report: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\README.md` and `series_elastic.py` in that folder.
- Related static/contact initialization findings: `MOMENT_FORCE_AUDIT.md`, `STANDING_INITIALIZER_AUDIT.md`, and `CONTACT_INITIALIZER_AND_ALIGNMENT.md` in that same folder.

The historical robot controller uses a joint-torque allocator and `motor_target_current`, which cancels background MN synaptic currents to track the allocator's targets. Preserve that case as a clearly labeled diagnostic comparator. In the neural-coordination test, remove that cancellation and independent runtime force allocation. Preserve meaningful RG/PF/MN and sensory contributions. An allocator may be an explicitly labeled offline teacher or initialization aid, but its own performance must not count as neural coordination. A prescribed harness remains acceptable if its forces and torques are measured and disclosed.

### Run bounded experiments before broad optimization

1. Check forward BPA force curves and any inverse map against Ben's static data across the relevant pressure/contraction range. Test pressure bounds, force direction, units, slack, and numerical stability. Do not silently extrapolate beyond characterized ranges.
2. Establish a simple nominal pressure-response model. Add separately configurable filling/venting response, delay, force-map mismatch, and loading/unloading effects only to the extent supported by data. Label assumed parameters as assumptions and sweep plausible ranges rather than presenting them as measured values.
3. Evaluate neural feedback on versus selected feedback ablations under identical actuator conditions. Keep baseline parameters and disturbances matched, and record actual sensory inputs, RG/PF/MN states, pressure commands, chamber pressure, and delivered forces.
4. Fit neural feedback on a specified training set of actuator variations. Evaluate held-out variations and disturbance seeds without further tuning. Demonstrate compensation for effects actually included in the evaluated model. Do not claim hardware robustness from nominal-model success alone.
5. Progress from suspended movement and supported standing toward actual harness-supported ground steps. Use the established stage definitions and Hunt training procedure; do not call suspended stepping stage-4 ground walking. If testing a Gait2392 branch instead, identify the model/actuator configuration separately and do not transfer robot outcomes to it.

Preserve bilateral parameter tying while fitting the nominal controller, including shared mirrored connection strengths. Deliberate unilateral actuator disturbances can be held-out robustness tests and must be explicitly identified. Retain separate sartorius/gracilis pathways and the selected RF/hamstring PF architecture; do not silently change the circuit to obtain a pass.

### Mechanical and success requirements

- Robot anatomical axes: x anterior, y superior, z right; correct gravity and ground orientation. Keep the confirmed ball-and-socket lumbar joint and 0.25 m rigid spine unless Ben changes them.
- Enforce no knee hyperextension for new runs. Check ground placement before loading; do not accept an initial interpenetration impulse as push-off.
- Log harness forces/torques and actual supported body-weight fraction, not just the configured support setting.
- Measure stance/swing loads, sustained whole-foot clearance, alternating contacts, crossings, progression, joint limits, and neural activity. A dragging, frozen, crossed-leg, or NaN trial is not successful stepping.
- Distinguish feedback compensation from success caused by increased harness assistance, changed initial placement, greater actuator capacity, or a hidden force allocator.
- Preserve MuJoCo 2.3.7, implicitfast, and the established 2 ms setting for comparable trials unless an isolated convergence study justifies a separately labeled variant.

### Timing, records, and deliverables

Benchmark warmed and cold timing separately, and report both neural-kernel and complete simulation throughput against matched reference calculations. Validate numerical agreement where the accelerated path is intended to be exact. Do not transfer isolated muscle-bank speed ratios to complete-walker claims.

Put new raw JSON/NPZ/log outputs in a dedicated gitignored run folder. Keep source, small curated configurations, and readable reports reviewable. Use a **separate BPA-neural-coordination experiment record**, and add only a short handoff entry to CHATGPT_REPORT.md. Preserve the existing Gait2392 training record and continuation prompt.

Deliver: actuator-interface description; actual active circuit diagram; exact runnable commands with machine and absolute paths; bounded nominal and held-out results; matched feedback ablations; actual assistance measurements; runtime benchmarks; and a short animation with a caption stating the controller and actuator model actually used. Report failures and remaining limitations plainly.

## Existing physics subagent status

The previous `tendon_damping_physics` subagent completed its work. It built an isolated, accelerated compliant biological-tendon prototype, checked static OpenSim agreement and damping/contact behavior, and investigated contact-aware initialization. It did **not** implement BPA hysteresis/pressure dynamics or validate a complete compliant-tendon SNS walker. Read its reports before reusing it; biological tendon compliance and BPA actuator dynamics are separate modeling questions.
