# Standing initializer audit and isolated corrected callable

easteregg2, 30 October 2026. No runner/Hunt edits, no force-law JSON changes.
Evidence: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\runs\audit_20261030\STANDING_INITIALIZER_AUDIT.json`.
Reproduction: `standing_initializer_audit.py` in this directory.

## Confirmed implementation problems

`runner.solve_standing_activations`, lines 185–244, constructs its design from
positive `data.actuator_moment.T * Fmax`. Native muscle force is negative along
the tendon transmission, so this has the wrong sign for positive activation.
It also assumes maximal force rather than the actual pose-dependent active
capacity and does not subtract the native muscle passive bias force. That
passive muscle force resides in `qfrc_actuator`, not `qfrc_passive`.

It fits 99 of 105 generalized rows, including equality-dependent pathpoint
slides, without reducing through E. The raw objective consequently mixes slide
forces (N) and hinge torques (N·m); its norm is not an anatomical torque norm.
Equality-aware reduction converts the dependent force effects into the correct
anatomical generalized torques.

After cropping rows, the preload indexes the cropped vector using original
joint qpos addresses. The current all-scalar model has nq=nv, so the principal
observed error is cropping without remapping, rather than a quaternion offset.

| Intended preload | Actual destination in this model |
| --- | --- |
| Right hip flexion, −40 | Right ankle angle |
| Right knee angle, −60 | `semiten_r_semiten_r-P2_x` follower slide |
| Right ankle angle, +40 | `semiten_r_semiten_r-P2_y` follower slide |
| Left hip flexion, −40 | Left ankle angle |
| Left knee angle, −60 | `semiten_l_semiten_l-P2_x` follower slide |
| Left ankle angle, +40 | `semiten_l_semiten_l-P2_y` follower slide |
| Lumbar extension, −20 | `psoas_r_psoas_r-P3_y` follower slide |

These targets are imposed net torques/forces, not a general implementation of
antagonist coactivation. The isolated corrected callable does not reproduce them.

## Pose and comparison definition

Two ground-on poses use the current repaired model and a solve-style support
rig (`kxy=1500`, other standing rig defaults). Both have the lowest foot mesh
1 mm below the actual ground plane, qvel=0, and pelvis orientation zero.
Neutral has the selected anatomical joint coordinates zero. The bent pose
uses native bilateral hip −7°, knee −10°, ankle 0°, and hip-adduction 0.04 rad.
No 70% support force is invented in the main comparison.

At zero activation, the audit captures muscle passive forces and the existing
constraint forces. It then measures actual active force capacity from the
activation-1 minus activation-0 native force response at that same pose.
The corrected matrix is `E.T @ (-actuator_moment.T * active_capacity)`.
The corrected target subtracts native muscle passive generalized force.
The linear optimization freezes zero-activation contact forces; it is therefore
a static initialization approximation, not a proof of equilibrium.

## Measured anatomical torque residuals

All entries below are L2 norms over the seventeen selected rotational
coordinates in N·m, not norms over mixed follower/root rows.

| Variant | Neutral, frozen constraints | Neutral, recomputed dynamics | Bent, recomputed dynamics |
| --- | --- | --- | --- |
| Legacy with preload | 45.23 | 44.37 | 45.27 |
| Legacy without preload | 18.34 | 20.51 | 22.42 |
| Legacy algebra, force sign only corrected, no preload | 12.38 | 15.51 | 11.30 |
| Corrected capacity/passive/projection, free activations | 6.06 | 9.53 | 6.33 |
| Corrected capacity/passive/projection, bilateral ties | 6.36 | 9.61 | 6.69 |

The corrected tied solve has exactly equal L/R activations, 46 free parameters
for 92 actuators, and activation bounds enforced within optimization. At
neutral its activation range is approximately 0.00079–0.34763.
It uses the source-style posture prior with ridge 0.25 and torque scale 0.01;
retaining that regularizer leaves a nonzero residual rather than demanding
exact torque cancellation at any activation cost.

**A smaller torque norm does not guarantee smaller acceleration everywhere.**
At neutral, maximum selected-coordinate acceleration is about 660 rad/s² for
the corrected tied result, versus 438 without legacy preload and 2644 with
legacy preload. At the bent pose the corresponding values are 397, 400, and
2556 rad/s². Small distal inertias and regularization matter. No dynamic
rollout or causal walking improvement is claimed by this audit.

Historical measured curriculum passes retain their original initializer and
provenance. A corrected-initializer variant must rerun its own stages; these
static results neither invalidate those observations nor establish a new pass.

## Callable for an explicit opt-in variant

File:
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\standing_solver.py`

Import within the current spinal/functional_gait2392 environment:

```python
from physics_fidelity.standing_solver import solve_standing_activations

x = solve_standing_activations(
    model, data, Fmax=None, preload=False, bilateral=True
)
```

Signature:

`solve_standing_activations(model, data, Fmax=None, preload=False, *, bilateral=True, return_diagnostics=False)`

Default return is the 92-element activation vector. With diagnostics enabled,
return is `(x, diagnostics)`. `Fmax` is accepted for runner-call compatibility
but ignored because the actual force response is measured. `preload=True`
raises: the caller must explicitly disable the historical artificial preload.
No production function is replaced by importing this module.

The function includes existing `data.qfrc_applied` and maps any already-present
body wrench `data.xfrc_applied` using `mj_applyFT`; it does not assume or invent
support from a trial configuration. An already-applied 70% BW pelvis-force
test passes. Incoming activation and velocity are restored, and forward caches
are recomputed. The reusable callable matches the isolated audit result to
within 1e−8 activation. Bilateral activation equality does not imply identical
delivered forces under the current asymmetric force surrogates.

The intended next test is an explicitly labeled initializer variant with the
same controller/model/support, followed by causal trajectory and curriculum
checks. Do not silently substitute it in the accepted baseline.
