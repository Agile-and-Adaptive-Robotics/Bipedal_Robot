# Isolated compliant-tendon and damping investigation

Machine: easteregg2. Python: `D:\Anaconda\envs\myo\python.exe`.
MuJoCo: 2.3.7, implicitfast, 2 ms. OpenSim Python: 4.4.1.
No production runner, neural controller, model, timestep, or trained winner was changed.

The new [full-Gait2392 replay experiment record](EXPERIMENT_RECORD.md) advances
the prototype to a bounded full-model mechanics comparison. It records matched
saved neural excitation, rigid/compliant controls, bilateral physical ties,
strict knee guards, fresh-process parity, whole-mechanics runtime, timestep
refinement, and a recovered native implicit muscle-velocity tangent. It is
separate from the earlier isolated bank/one-joint results below. Its later
sections also report actual live-SNS supported standing and stepping checks.

The current section of that record includes the separate `thelen_dynamic.py`
source-matched force–velocity extension and constitutive velocity tangent.
Actual OpenSim dynamic checks pass, and an opt-in live SNS wrapper preserves
accepted constrained syn6 standing under rigid and compliant source laws,
with an exact fresh repeat. The unchanged Hunt controller initially loses
ground clearance with both source laws. Shared bilateral recruitment scales
now restore passing air1/air2 and supported-ground stepping with both laws.
The compliant case has an exact fresh-process repeat; its bounded parameter
basin is recorded separately. See external
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\hunt_motor_calibration_v1\HUNT_CALIBRATION_CHECKPOINT.json`.
This is an opt-in human Gait2392 example, with measured apparatus assistance,
not a production/default change. The earlier approximate prototype below is
retained unchanged.

The follow-up [matched-coordinate route and shared-force-law audit](MOMENT_FORCE_AUDIT.md)
finds largely symmetric routes but substantial low-activation force asymmetry,
and supplies a jointly fitted, opt-in six-family static surrogate parameter file.
The [standing initializer audit](STANDING_INITIALIZER_AUDIT.md) also supplies an
isolated corrected callable with sign, passive-force, equality-projection,
bounded-solve, and bilateral-tie corrections. Its static improvements do not
establish dynamic balance or walking.
The [contact initialization and ground-alignment follow-up](CONTACT_INITIALIZER_AND_ALIGNMENT.md)
finds severe initial foot penetration, provides a repaired-probe alignment
helper, and distinguishes successful fixed-activation holds from the still
unsuccessful SNS transfer of the new contact-aware initializer.

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

### Reference-constrained stance-load diagnosis

On easteregg2, the separate current record is
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\hunt_ib_load_v1\EXPERIMENT_RECORD.md`,
with numerical results in `CHECKPOINT.json`. The source-compliant engineering
clearance result does not pass the saved IK envelope; no lower-support gait
claim follows from it. The bounded new Ib study finds no fully accepted case.

Reusable additions are `biomechanical_envelope.py` (frozen reference and route
gates), `stance_ib.py` (opt-in source-normalized stance load pathway and bounded
Rustuna adapter), and `pt_current_chain.py` (shadow neural/physical transfer
audit without modifying the live standing trajectory). `balance_diagnostics.py`
provides current reconstruction, exact-array comparison, and checkpoint tools.

The subsequent source-afferent/hip–knee record is
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\hunt_spindle_v1\EXPERIMENT_RECORD.md`.
`source_force_comparison.py --actual-torque-audit --reference-model <captured-model-dir>
--folder <full-dynamics-replay> --out <external-output-dir>` reconstructs actual
source active/passive torques, separates contact/stop reactions, checks conditional
vasti capacity, and compares nine poses with original OpenSim. It changes no
controller or trajectory. `--reference-model ... --prepare-contact-reference`
declares offline pelvis-height diagnostics; these are not runtime clamps or an
accepted ground reference. The next requirement is contact-feasible ground
training data, not further ankle-only optimization.

### Native contact calibration checkpoint

`supported_step_teacher.py` now prepares explicitly declared in-place engineering
references and native soft-contact standing/keyframe calibrations. The current
external record is
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\supported_teacher_v1\EXPERIMENT_RECORD.md`.
Read `PLANTAR_CONTACT_AUDIT.md` before treating the original-source92 force-law
variant as equivalent to the legacy runner's six-muscle pruning policy.
Native inverse contact balance, source recruitment, and actual forward holds
are separate checks. The current successful single-support holds use 60% nominal
upward support and are not neural-controller or alternating-step results.

### Continuous source-muscle teacher

Record and gates:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_transition_v1\EXPERIMENT_RECORD.md`
and its sibling `CHECKPOINT.json`. Both swing-initiation references and the
complete left lift/return pass fresh/timestep checks; the opposite complete
return remains rejected. This is not an SNS or alternating-gait result.

On **easteregg2**, reproduce the accepted left reference into a new output
folder with PowerShell:

```powershell
$physics = 'D:\Users\Ben Bolen\NeuroCampaign_20261031\physics'
& 'D:\Anaconda\envs\myo\python.exe' `
  'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\physics_fidelity\source_force_comparison.py' `
  --reference-model "$physics\hunt_spindle_v1\reference_inverse_v1\actual_model_binary" `
  --transition-start "$physics\supported_teacher_v1\native_ds60_source_hold_v1" `
  --transition-via "$physics\native_transition_v1\unloaded_left_nearflat_source_hold" `
  --transition-end "$physics\supported_teacher_v1\native_right_ss60_fresh" `
  --spec "$physics\hunt_ib_load_v1\IK_ENVELOPE_SPEC_V1.json" `
  --transition-ramp 4 --transition-roundtrip --return-dwell .5 `
  --transition-feedback --command-filter .05 --feedback-warmup `
  --mtp-impedance-scaled `
  --out "$physics\native_transition_v1\new_left_return_01"
```

The default is the captured 2 ms physics step. `--hold-dt .001 --trace-stride 20`
tests 1 ms while preserving the 10 ms teacher clock. `--export-training-reference`
adds one full compressed physical reference; use it for selected final cases,
not exploratory trials. Completed output folders are protected against reruns.
The source92/pruned actuation policy must match every endpoint. Shared MTP
teacher impedance is computed from existing model stiffness and inertia;
physical joint damping, contacts, source anatomy, and production code are
unchanged. The accepted final reference exports actual excitation separately
from end-step activation and the split-step activation used for source force.

### Stationary-endpoint prerequisite

The two-second holds above are bounded mechanical examples, not demonstrated
stationary equilibria. They do not authorize bilateral Hunt C training. See
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_equilibrium_v2\EXPERIMENT_RECORD.md`
for the newer full-native force balance and five-second validation. Every
native follower residual is retained. New endpoint formats are rejected by
transition training until `accepted_stationary_endpoint` is true.

The first double-support endpoint passes all static and five-second gates,
with a bit-identical fresh repeat. Its authoritative validation folder is
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_equilibrium_v2\double_support_source92_jacobian_scaled_5s`.
This remains an original-source92, 60%-supported mechanical reference.
It changes no source dimensions, native stiffness/damping, or contact geometry.

Use `source_force_comparison.py --full-native-equilibrium` to fit a new endpoint;
`--resume-equilibrium` restores saved variables while retaining the original
target. `--static-force-priority` changes numerical objective weighting only.
Run `native_equilibrium_validation.py --model <captured model folder> --endpoint
<static fit folder> --spec <frozen envelope JSON> --out <new validation folder>`
before transition use. `--physics-dt .001` is an isolated refinement check;
the default captured physics step remains 2 ms.

All three endpoint prerequisites now pass: double support, left stance, and
right stance, each with exact fresh repeats and passing 1 ms checks. See the
external `native_equilibrium_v2\CHECKPOINT.json` for authoritative folders.
`--static-forward-acceleration` conditions the inverse fit with all actual
native accelerations; `--feasibility-first` prioritizes declared physical goals
over a soft posture preference. Neither changes the model or acceptance gates.
The validator's `--initial-audit-only` records native forward/inverse consistency
without another forward trial.

Both continuous no-reset lift/return baselines clear the swing foot and recover
double support, but touchdown slip remains above the 5 mm/s limit. Dynamic
Hunt C fitting stays gated. Before using `NativeCartesianTeacher`, preserve
loaded follower offsets in its reference and trial FK; its current ideal
projection would discard the newly qualified geometry. No further physical or
neural variant has been promoted.

### Loaded Cartesian target repair

The opt-in `--loaded-cartesian-followers` now preserves actual keyframe and
live follower geometry/velocity. Add `--root-aware-cartesian` and optionally
`--stance-posture-cartesian`; `--cartesian-tangential-clearance .005` separates
horizontal/orientation settling from vertical touchdown without changing the
stationary endpoints. This adapter explicitly uses the legacy captured human
z-up frame, not a reinterpretation of the robot's anatomical y-up frame.

`validate_loaded_cartesian.py` checks legacy/default-off replay, exact loaded
FK, and independent world-velocity/profile derivatives. Full default-off
physical replay is also exact. Current gates and rejected trials are in
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_touchdown_v3\CHECKPOINT.json`.
No teacher passes all physical gates. `--native-impedance-acceleration` is an
explicitly rejected experimental recruitment variant, not an adopted controller.
No source dimensions, physical damping, masses, or contact settings changed.

`predictive_source_recruitment.py --trace <existing full trace> --model <captured
model folder> --out <new folder>` validates native activation and source-force
forecasts without another physical run. `--snapshots <actuation_windows.npz>`
instead checks complete augmented-state replay and its command Jacobian.
Capture at most three short windows with `--actuation-audit-time <seconds>`.
The experimental `--predictive-source-horizon .01` recruiter is rejected on
physical gates: accurate conditional force prediction does not replace the
full q/velocity/activation/fiber/filter transition. It remains default-off.

Add `--inverse-reference-jacobian <jacobian_windows.npz> --spec <frozen ROM
JSON>` to the snapshot audit for preallocated full-state derivatives,
conditioning, bounded command authority, and one small regularized correction
only when derivative checks pass. `--branch-derivatives` instead diagnoses
one-sided derivatives and always disables the inverse. Current equilibrium
directional checks fail; no live teacher uses this method. Read the separate
`native_touchdown_v3` external experiment record for the measured scope.

`--stationary-start <accepted DS> --stationary-target <accepted left stance>
--spec <frozen ROM JSON>` runs the separate bounded nonlinear right-to-left
load-transfer diagnostic. It uses shared functional-family corrections and
native augmented rollouts. The initial 10 ms task fails; one analytically derived
302 ms horizon passes a 1% BW transfer, fresh-process and 1 ms checks included.
`--verify-stationary-inverse <saved inverse folder>` replays saved corrections
without fitting; `--physics-dt .001` performs the isolated refinement check.
This is an offline mechanical primitive, not a supported-step or neural pass.

Add `--continuation-seed <preserved1% inverse>` for the state-continuous
1→3→5% BW sequence. Both directions pass through 3% and stop at 5% under the
unchanged diagnostic reference hold. `--verify-continuation <saved sequence>`
replays only its accepted prefix with no fitting; an optional `--physics-dt
.001` checks refinement. `--audit-transfer-geometry <saved sequence>` separately
decomposes the loaded foot-body motion without another simulation search.
The stopped 5% result is a reference-position limitation, not a demonstrated
muscle-strength deficit. Read the separate external record before continuing.

`--rolling-stance` selects the distinct contact-consistent stance task. It
allows body-origin motion through rolling while retaining actual contact-slip,
native friction/unilateral-force, mesh, source-envelope, root, and knee checks.
The existing inverse-teacher subtalar/MTP +/-10 degree posture envelope is
retained. Independent contact-kinematic checks and exact default-off replay
precede the new runs. Both rolling sequences pass 1→3→5% BW transfer, fresh
exact replay, and unchanged-command 1 ms gate checks. Right-toe acceleration
remains resolution-sensitive. Both feet are still loaded; this is not a swing
or full stepping-teacher pass. Historical body-position-hold failures remain.

`--extend-transfer-from <accepted5% rolling sequence>` freezes that exact
continuous prefix and attempts 7.5% then 10% BW with the unchanged 302 ms
horizon and constraints. Both directions pass through 10%, fresh exact replay,
and unchanged-command 1 ms gate checks. Both feet remain loaded, and detailed
toe acceleration remains timestep-sensitive. Use the separate external record
before treating these load-transfer primitives as a stepping reference.

Extending an accepted 10% sequence with `--extend-transfer-from` attempts
12.5% then 15% with the same policy. Both current directions pass 12.5% and
stop at 15%: right force residual, left vertical-root drift. Accepted 12.5%
prefixes pass fresh/1 ms checks but are not stationary endpoints. Separate
actual-force and whole-foot-mesh unloading scores remain false.
`--audit-transfer-authority <saved attempt>` diagnoses group command caps,
delivered recruitment, root displacement, and contact ownership without a search.


## Stationary loaded-reference validation
The source92/60%-support12.5% load-split references now pass full-native105-coordinate static and5-second/fresh/1-ms checks. See external native_stationary_loadsplit_v1/EXPERIMENT_RECORD.md. Simple full-source activation interpolation fails dynamic gates; no swing or Hunt C claim. native_equilibrium_validation.py --campaign-summary curates repeatability and failed-motion evidence; loaded_contact_min_BW is an optional native_transition audit threshold with historical default .025 preserved.


## Full-family derivative and branch audit
On easteregg2, predictive_source_recruitment.py --full-family-audit --terminal-window --physics-dt .001 --finite-difference-step .001 --paired-rollout-budget 53 retains all 46 task-normalized family directions and deliberately performs no fit. Supply the frozen model,start,two stationary targets,spec,and external output directory through the existing flags. native_equilibrium_validation.py supplies --native-branch-audit,--plot-native-branches,and a zero-effect --predictor-branch-target/--predictor-branch-start shadow audit.
The completed example is external physics/native_manifold_v1/full46_1ms_audit_v1. Right derivatives are consistent; left force/velocity blocks fail. Native constraint states remain fixed in the baseline,while a source-predictor FV branch mismatch is localized. The physical source law and production defaults are unchanged; no movement is promoted.
