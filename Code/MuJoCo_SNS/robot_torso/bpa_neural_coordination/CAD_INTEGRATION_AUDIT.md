# CAD BPA integration audit and force-direction correction

Session 1 BPA subagent, **easteregg2**, 2026-10-31. This extends the separate EXPERIMENT_RECORD.md. No full-CAD simulation is claimed.

## What exists, and what remains unassigned

The current `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\ah001_neural.xml` contains **54 active reference-tension actuators, 27 per side**, selected from 92 source routes. The model mass is 75 kg and gravity is `[0,-9.81,0]`. Geometry and moving/conditional path functions exist. Those are usable interface inputs.

`robotbody_routes.json` explicitly calls the force model **unassigned**, and says source human Fmax values are provenance rather than robot actuator capacity. The reference actuator caps come from the master's MIF targets. The available human optimal fiber and tendon slack lengths are biological source parameters; they are not BPA sleeve/cable dimensions.

No complete route-specific BPA parameter assignment is present for any of the 54 active CAD actuators. The missing assignment comprises:

- sleeve resting length and usable contraction fraction;
- nominal BPA diameter and number in parallel;
- nonmuscle path length (cable, two fittings, and any separately reconciled constant offsets);
- effective series compliance/stiffness and its modeling scope;
- intended pressure ceiling;
- identified or explicitly assumed fill/vent response and delay;
- provenance tying those values to this retargeted CAD route.

The old **Routing27_v2** result folder does contain historical optimized designs for **15 of the 27 master actuator IDs**. Their MAT files have readable resting length, tendon length, 20 mm diameter, 620 kPa operating pressure, fitting length, and selected single/parallel/tackle configuration. They are useful **candidate design sources**, not nonexistent data. However, they were solved on older robotbody primary-DOF routes, with moved attachment points and their own compliance and pulley conventions. Some selected candidates have large target-torque deficits. They have not been checked against the present AH001 CAD path lengths and three-dimensional motion.

Missing historical design IDs in that folder are **1, 2, 3, 15, 17, 18, 19, 20, 21, 22, 23, and 25**. The readable IDs are 4–14, 16, 24, 26, and 27. A historical selected design is not automatically a feasible mounted-CAD design.

The extraction is reproducible with `inspect_routing_candidates.py`; detailed values and exact MAT paths are in `runs\cad_parameter_audit_20261031\HISTORICAL_ROUTING_CANDIDATES.json`. No historical Xi values have been adopted or refitted. The design-map/readme values are historical, not an override of the dissertation's current adopted compliance parameters.

For comparison, the old W2L BPA test **constructs hypothetical actuator sizes** using the largest sampled route length minus a cable length plus a bias, with supplied diameter/count/contraction settings. It does not provide measured CAD actuator dimensions. See `Code\MuJoCo_SNS\spinal\w2l_mujoco\test_bpa_stepping.py` around the BPA construction loop.

## Implemented opt-in plumbing

`cad_interface.py` provides `RoutedBPABank` with explicit per-route parameters, actual pressure state, separate fill/vent response, command delay, positive tensile force, and either direct-pressure or inverse-force command conversion. It never edits neuronal voltages or cancels synaptic currents. It rejects incomplete configurations. A two-unit direct-pressure test passed, as did the refusal of the unconfigured 54-actuator CAD template.

`CAD_BPA_PARAMETER_TEMPLATE.json` inventories all 54 active actuators, their current path lengths, and their reference metadata, with required hardware/assumption fields left **null**. No sleeve lengths, counts, or dynamic properties were invented. The template is deliberately not runnable until completed or populated with clearly identified design assumptions.

The caller must explicitly disable the old reference-force gains on an isolated in-memory model before applying these BPA forces. `assert_native_force_disabled` checks that condition. The caller can then obtain ordinary SNS motor-neuron outputs, update `DynamicRoutes`, compute the pressure/BPA response, and add the bank's generalized force to its existing support-force calculation. Native reference tensions must not run simultaneously.

The bank currently supports the Festo static surface, finite **total assembly** axial series stiffness, and a constant nonmuscle-length offset. Its stiffness input is not automatically multiplied by the parallel count. Bend-dependent Xi3 losses, directional bracket compliance, tackle transmission, and any identified dynamic hysteresis still require explicit reconciliation/implementation. The single neutral-pose route length in the template is an inventory measurement, not the complete RoM envelope or a basis for silently inferring sleeve length.

For the CAD moving-point routes use **`DynamicRoutes.length_jacobian`**, including the moving-site contributions, and apply:

`Q_BPA = -length_jacobian.T @ positive_tensions`

The bank does not clear an external harness's existing generalized forces or install a global MuJoCo callback. The caller owns force-buffer resets and physics scheduling. No heavily supported full-CAD trial was run, because it would currently require unrecorded actuator-sizing choices rather than merely interface plumbing.

## Confirmed defect and active correction

The existing shared helper `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\bpa_mujoco.py` previously added **`+J_length.T @ positive_tension`**. This reverses tensile action. Its comments also inferred matrix orientation from the first dimension alone, which is ambiguous when the number of actuators equals the number of DoFs.

The helper is now corrected to **subtract** the positive tension contribution and prefer MuJoCo 2.3.7's documented `(nu,nv)` actuator-moment layout. A nonsquare transposed compatibility fallback remains. Its original ownership of `qfrc_applied` is unchanged; callers must still combine other generalized forces explicitly.

The pre-fix file is preserved at:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\runs\cad_parameter_audit_20261031\bpa_mujoco_before_sign_fix.py`

Pre-fix SHA256: `dbd10111a7a8b2fff3ca6ea0a9252ab3206649ca3011d7fc10b74eec540010a9`.

### Meaningful regression checks

`check_bpa_glue_parity.py` compares the **actual corrected helper** against native MuJoCo tendon actuators applying negative scalar force for positive tension. Seven cases cover positive and negative tendon-length gradients, shortening and lengthening velocities, square and nonsquare moment matrices, and repeated callback calls. Joint forces and mechanical power agree to the asserted numerical tolerance. Power independently satisfies `Q·qdot = -T·Ldot`.

The simplest flexor example has `dL/dq=-0.03 m/rad` and T=100 N: native MuJoCo produces **+3 N·m**; the old positive-sign mapping gives **−3 N·m**. This is a force-direction issue, not an argument about whether a neuron-to-actuator conversion is legitimate.

### Affected callers and historical evidence

Direct consumers found in the repository:

- `Code\MuJoCo_SNS\sns_bpa_demo.py`;
- `Code\MuJoCo_SNS\spinal\w2l_mujoco\test_bpa_stepping.py`;
- `Code\MuJoCo_SNS\spinal\w2l_mujoco\_curriculum_bpa.py`;
- the archived campaign copy `Code\MuJoCo_SNS\spinal\campaigns\20260930\bpa_run\test_bpa_stepping.py`.

Indirect consumers include `_curriculum_bpa_spiking.py`, `_bpa_runbest_20261003.py`, and `Code\MuJoCo_SNS\run_w2l_bpa_extended.py`. Historical BPA trials through these paths need a **sign-corrected replay** before they can support physical-control claims. Their existing results have not been overwritten, and no dissertation numbers or figures have been silently revised.

**`Code\MuJoCo_SNS\spinal\runner.py` does not import or call this helper.** It supplies native actuator controls and calls MuJoCo stepping. The CAD robot's historical `run_neural_robot.py` also uses native reference tensions and its allocator rather than this BPA helper. Those simulations are not invalidated by this particular glue defect, though their own documented limitations remain.

### Specific Appendix D BPA example: dependency confirmed

The historical **1.542 s period, 63.5° left-knee range, peak foot loads 880/1372 N, and 372 kPa pressure** all come from the gate documented in `Code\MuJoCo_SNS\spinal\campaigns\20260930\bpa_actuation.md`, section 1.3, and its `bpa_run\bpa_gate.log`. The documented command ran the saved `bpa_run\test_bpa_stepping.py`. That driver imports the shared `BPAMuscleSystem` and directly installs **`mujoco.set_mjcb_control(system.cb)`**; there is no compensating sign reversal in the driver. Therefore this dissertation example follows the affected force-application path and requires explicit sign-corrected replay before it can support a claim of physically correct tensile BPA actuation. Its old numerical values and source logs remain preserved.

The report's discussion of the knee-extensor route producing a flexion-sign torque treats the **length gradient** as the direction of a positive tensile force. Its statement that Hill-muscle runs carry that same torque sign is not established by sharing tendon geometry: native Hill muscle force is negative in MuJoCo's transmission convention. This is precisely the distinction checked by the new native-actuator regression.

`CROSS_PLATFORM_REPORT.md` repeats those same BPA figures and gate verdict; it is downstream reporting, not an independent replication. The shared helper has one prior Git commit, **854a7923 (2026-09-08)**, whose source contains the same positive-sign force application and predates the reported gate. No source hash embedded in the original gate log was found, so the dependency assessment rests on the preserved driver, shared-helper history, and the documented positive-force callback convention rather than a recovered executable-state hash.

Keep the later **native Hill-muscle W2L feedback/wean trials** separate from this BPA gate. Sharing a body/route family does not make them consumers of the defective BPA helper. Likewise, the current Gait2392 spinal runner is not a consumer. Do not replace or invalidate those separate results because of this dependency finding.

## Commands on easteregg2

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\check_bpa_glue_parity.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\audit_cad_integration.py'
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\inspect_routing_candidates.py'
```

First reconcile candidate/design assumptions with the current CAD routes, then test an explicitly configured neural-only heavily supported robot. Do not bypass missing dimensions by reusing human fiber lengths or MIFs as if they were identified BPA hardware.
