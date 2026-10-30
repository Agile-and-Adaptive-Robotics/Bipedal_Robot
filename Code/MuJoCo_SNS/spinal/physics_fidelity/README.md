# Isolated compliant-tendon and damping investigation

Machine: easteregg2. Python: `D:\Anaconda\envs\myo\python.exe`.
MuJoCo: 2.3.7, implicitfast, 2 ms. OpenSim Python: 4.4.1.
No production runner, neural controller, model, timestep, or trained winner was changed.

The follow-up [matched-coordinate route and shared-force-law audit](MOMENT_FORCE_AUDIT.md)
finds largely symmetric routes but substantial low-activation force asymmetry,
and supplies a jointly fitted, opt-in six-family static surrogate parameter file.
The [standing initializer audit](STANDING_INITIALIZER_AUDIT.md) also supplies an
isolated corrected callable with sign, passive-force, equality-projection,
bounded-solve, and bilateral-tie corrections. Its static improvements do not
establish dynamic balance or walking.

## Findings that matter for the current model

1. The OpenSim source uses **Thelen2003Muscle**, not Millard2012. Its tendon
   strain at maximum isometric force is 0.033. The native MuJoCo muscle has
   an inelastic biological tendon, no explicit pennation state, and different
   force–length–velocity curves. MuJoCo's spatial *tendon* is the routed total
   muscle–tendon length, not an elastic biological tendon by itself.
2. The parent session traced bilateral strength differences to independent
   Step-3 surrogate fits. Source soleus Fmax is 3549 N on both sides; fitted
   scale factors 0.9567481/0.89782008 produce approximately 3395.5/3186.36 N.
   These are jointly fitted surrogate curve parameters, not evidence that the
   original anatomical muscles differ in strength. Equalizing that one number
   without revalidating the complete curves is not a faithful conversion fix.
3. Equality-follower armature does add anatomical effective inertia, but the
   current `bsolve_ik.py` docstring's approximately 36 kg·m² knee claim is not
   supported by the measured projection below. Do not repeat that number.
4. `runner._apply_contact_damp('nonlinear')` changes `solimp[3]`, which is the
   **midpoint**, while calling it width. Width is index 2. Its index-1 change
   changes impedance at width, not midpoint. With default d0=0.9 and the new
   d_width=0.9, impedance becomes constant, so width/midpoint changes have no
   effect. The existing name and historical interpretation are misleading.
   This investigation does not alter that historical option or its results.

## Implemented prototype

`series_elastic.py` adds a projected fiber-length state and an explicit elastic
biological tendon, with constant fiber thickness to represent pennation. It
uses source-informed Thelen tendon and static fiber force–length curves.
The force–velocity curve is still an approximation (`1+tanh(v_norm/10)`), and
the added normalized fiber damping beta=0.05 is an experimental regularizer,
not a property recovered from the original Thelen source. Activation is an
input; this module does not replace activation dynamics.

For fiber projection s, fiber length l=sqrt(s²+h²), and cos(alpha)=s/l:

`Lt = Lmt - s`

`0 = [a*fL(l/l0)*fV(v/l0) + fP(l/l0) + beta*v/l0]*cos(alpha) - fT(Lt/Lslack)`

`v = cos(alpha)*(s_new-s_old)/dt`

A bracketed backward-Euler root updates s. The tendon has the Thelen
exponential toe region and linear continuation, normalized to force F0 at
the source's strain-at-F0 parameter. The implementation uses Numba and
preallocated state/force/residual arrays for the full 92-muscle bank.
Bound hits and force residuals are exposed rather than silently accepted.

This is a tested implementation scaffold, **not a full Thelen/Millard/Hyfydy
replication or a validated drop-in replacement for the fitted walker**.

## Measured checks, 30 October 2026

Active raw results: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\runs\audit_20261030\`.
The folder-local `.gitignore` excludes `runs/`. Earlier scratch files disappeared
from repo `tmp/` for an unknown reason; `REGENERATED_EVIDENCE.json` identifies
the reproduced checks. Historical curated parameter files remain unchanged;
regeneration does not claim to recover identical historical raw-file hashes.
Timing figures below describe the original measured run; current rerun timing
is recorded in RESULTS.json and can vary with machine load.

| Check | Result |
| --- | --- |
| 92 compliant units, 16 simulated seconds, 2 ms | 2.47 s warmed wall time; muscle bank only, about 6.47× real time |
| Same approximate fiber law, rigid tendon | 0.0776 s; compliant bank about 31.9× slower in this isolated workload |
| Largest equilibrium residual | 2.29e−8 N; zero bound hits |
| Independent SciPy-root trajectory | Maximum projected fiber-length error 3.29e−14 m |
| Actual OpenSim static checks, six activation/pennation cases | Maximum force difference 2.16e−6 N; fiber-length difference 8.32e−12 m |
| OpenSim tendon at identical fiber state | Force difference 0.0 N |
| OpenSim static equilibrium at prototype state | Residual at most 1.12e−9 N |
| 2 ms force RMS difference from 0.5 ms reference | 0.1376 N for the prescribed soleus test |
| 1 ms force RMS difference from 0.5 ms reference | 0.0463 N, demonstrating refinement convergence |
| Actual one-DoF MuJoCo spatial-tendon coupling | Four seconds finite at 2 ms; force residual below 7.3e−10 N |

Timing is workload- and machine-specific, excludes compilation, and is not
whole-walker throughput. The compliant root overhead remains real despite
acceleration. The Python-vs-compiled outer-loop timing in RESULTS.json shares
the same compiled scalar solver; it must not be presented as a full interpreted
versus Numba speedup. The single-unit coupling uses MuJoCo's supported
`mj_step1/mj_step2` path, which preserves implicitfast in version 2.3.7.
External muscle/body coupling is partitioned; fiber-state implicit integration
does not make the full neuromechanical system monolithically implicit.

The static OpenSim check explicitly sets KshapeActive=0.5, KshapePassive=4,
FmaxMuscleStrain=0.6, and FmaxTendonStrain=0.033 to match the source. Constructor
defaults do not all match; the initial unaligned-default comparison differed
by up to 0.257 N. The reported near-exact checks use matched parameters.
No dynamic OpenSim force–velocity equivalence has been established.

## Damping and contact checks

`benchmark.py` verifies a 1 kg/100 N/m oscillator against its analytic solution.
At viscous damping 10 N·s/m, maximum position error is 0.0644 mm, energy decays
to 4.14e−9 of its initial value, and sampled energy never increases. This
checks the sign/implementation; it does not identify anatomical damping.

`contact_benchmark.py` compares a 1 kg sphere dropped 0.10 m onto a plane:

| Contact variant | Peak force | Maximum penetration |
| --- | --- | --- |
| Baseline solref [0.02,1] | 142.9 N | 8.73 mm |
| Existing less-viscous solref [0.01,0.7] | 281.0 N | 3.14 mm |
| Actual legacy “nonlinear” parameters | 148.9 N | 8.68 mm |
| Correct width index, retaining depth variation | 141.6 N | 8.75 mm |

The legacy variant and a differently shaped constant-impedance control have
**identical traces**, confirming the indexing/constant-impedance interpretation.
Lower penetration does not by itself mean better contact physics; contact force
and dissipation must be fitted/validated against appropriate data. These toy
impact results do not select a production setting or emulate Hunt–Crossley.

Fiber damping, joint viscous damping, contact damping ratio, and harness
damping are distinct quantities with different units and physical roles.
Do not increase them interchangeably to hide a controller/geometry error.

## Equality-follower inertia audit

`follower_inertia_audit.py` computes `E.T @ diag(follower_armature) @ E`, where E
is the full polynomial equality projection used by `CoupledArms`. Its derivative
is checked by finite differences. It also computes `E.T @ M @ E`, with pelvis
coordinates fixed, to express the follower contribution as a fraction.

| Pose | Right-knee follower contribution | Fraction of reduced knee inertia |
| --- | --- | --- |
| Neutral anatomical coordinates | 0.03256 kg·m² | 6.19% |
| Right hip 30°, knee −60°, ankle 10° | 0.02570 kg·m² | 5.35% |

There are 82 follower DOFs in this model. Forward follower armature is 1.0;
the actual ID mode sets it to **0.005**, not zero, and disables equalities.
At neutral, ID-mode reflected follower contribution is 0.0001628 kg·m².
No forward-model follower values or flags were changed. The current feasible
Hunt C fitting uses live forward-model teacher trajectories, so the difference
must not automatically be blamed on training against the ID model.

## Recommended integration sequence

1. Preserve the original model and correct the shared bilateral neural fit.
2. Compare source and converted total muscle–tendon lengths, pennation, full
   force curves, and effective moment arms at matched poses. Independently
   fitted length/force normalization must be accounted for before applying
   physical source parameters to the converted routes.
3. Replace the approximate FV/activation handling with the selected source law
   and validate force transients against OpenSim, including low activation,
   slack tendon, high pennation, and rapid lengthening. Keep damping explicit.
4. Add the external fiber states as an opt-in force backend. Preserve activation
   dynamics, sign conventions, and equality-aware routing; avoid double counting
   native active/passive forces. Restoring a simulation snapshot must restore
   fiber states too. Rollback/substepping without those states is invalid.
5. Run a fixed-control replay and timestep/bound/energy checks before fitting
   the neural controller to this new physics. Do not compare independently
   retuned controllers as if that isolated the engine difference.

## Commands

From `D:\GitHub\Bipedal_Robot`, run each with
`D:\Anaconda\envs\myo\python.exe`:

- `Code\MuJoCo_SNS\spinal\physics_fidelity\benchmark.py`
- `Code\MuJoCo_SNS\spinal\physics_fidelity\opensim_static_check.py`
- `Code\MuJoCo_SNS\spinal\physics_fidelity\contact_benchmark.py`
- `Code\MuJoCo_SNS\spinal\physics_fidelity\follower_inertia_audit.py`

## Primary sources

- [MuJoCo 2.3.7 muscle modeling and differences from OpenSim](https://mujoco.readthedocs.io/en/2.3.7/modeling.html#muscle-actuators).
- [MuJoCo 2.3.7 solver impedance and damping parameters](https://mujoco.readthedocs.io/en/2.3.7/modeling.html#solver-parameters).
- [OpenSim Thelen2003 source equations](https://github.com/opensim-org/opensim-core/blob/main/OpenSim/Actuators/Thelen2003Muscle.cpp).
- [OpenSim damped equilibrium formulation](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53090555/Millard%2B2012%2BMuscle%2BModels).
- [MuJoCo 2.3.7 split-step integrator implementation](https://github.com/google-deepmind/mujoco/blob/2.3.7/src/engine/engine_forward.c).
