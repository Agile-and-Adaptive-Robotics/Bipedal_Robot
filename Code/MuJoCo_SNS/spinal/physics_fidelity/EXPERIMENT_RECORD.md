# Gait2392 tendon and damping experiments

Machine: easteregg2. Python `D:\Anaconda\envs\myo\python.exe`, MuJoCo 2.3.7.
This record is separate from the handoff report and manuscript. Active raw
output is external: `D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\`.
One reusable driver, `tendon_replay_campaign.py`, owns this new campaign.

**Current:** source-matched dynamic Thelen forces and the constitutive velocity
tangent are validated. Live constrained syn6 standing passes with all three
force backends, including a bit-identical compliant repeat. Hunt ground
transfer loses clearance with both new source laws; see the current section
at the end. Earlier approximate-law results remain historical comparisons.

**Current:** source-matched dynamic Thelen forces and the constitutive velocity
tangent are validated. Live constrained syn6 standing passes with all three
force backends, including a bit-identical compliant repeat. Hunt ground
transfer loses clearance with both new source laws; see the current section
at the end. Earlier approximate-law results remain historical comparisons.

## Protocol and baseline

The first full-model comparison replays the same saved 92-channel neural
excitation sequence through retained MuJoCo activation dynamics. It does
**not** recompute SNS feedback after changing mechanics. Activations are
initialized identically from the first excitation sample. The saved source is
`hunt_training/runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/contact_feasible_start_v1/ground_guide.npz`.
The recorded 20-channel `act` array is a selected-muscle diagnostic, not the
full input; the full `ctrl` array is used instead.

All comparisons start from the same recorded pose and use the same repaired
converted Gait2392 geometry, equality followers, contacts, frontal restraints,
and 70% body-weight upward support. Vertical harness spring/damping are off;
horizontal/orientation restraints remain. At the primary 2 ms timestep,
`implicitfast` is unchanged. Knees have a zero-degree extension limit and a
2-degree protective margin; every physics state is checked, with immediate
rejection at any positive native knee angle. There is no position/velocity
clipping or hidden force allocation.

The external rigid/compliant pair uses the same Thelen-informed approximate
fiber force law and all source-derived physical parameters, explicitly pooled
between bilateral homologues. Only biological tendon compliance differs.
Normalized fiber damping is beta=0.05 or 0.10; this is an assumed parameter,
not an identified tissue property. Hip/knee/ankle joint damping is tied
bilaterally and tested at multipliers 1 and 2. The native fitted MuJoCo law is
a separate reference, not an identical constitutive law.

## Eight-second full-model screen

`screen_8s_v2/SUMMARY.json` and `VALIDATION.json` contain the nine cases.
All completed eight seconds with finite states, zero fiber-bound failures,
and no knee hyperextension. The maximum elastic equilibrium residual was
3.223e-8 N. Bilateral physical-parameter disagreement is exactly zero.

| Force backend | Joint damping | Fiber beta | Warm wall time | Peak single-foot normal load | Leg angular velocity RMS | Peak absolute pelvis tilt |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Native fitted reference | 1 | Not used | 2.689 s | 1.773 BW | 1.761 rad/s | 23.73 degrees |
| Source-informed rigid | 1 | 0.05 | 2.311 s | 0.750 BW | 0.837 rad/s | 24.03 degrees |
| Source-informed compliant | 1 | 0.05 | 3.895 s | 0.799 BW | 0.862 rad/s | 23.48 degrees |
| Source-informed compliant | 2 | 0.05 | 3.784 s | 0.911 BW | 0.859 rad/s | 23.47 degrees |
| Source-informed compliant | 1 | 0.10 | 3.659 s | 0.617 BW | 0.507 rad/s | 20.75 degrees |

Loads in this screen count calcn/toe contacts against the ground plane.
The mean combined measured load was approximately 0.298 BW for the external
rigid/compliant cases. Peak values in the table are the largest **single-foot**
normal load, while the plot shows the sum across feet.

Compliance alone did not reduce the baseline peak load: 0.750 to 0.799 BW.
Doubling joint damping increased that peak in the compliant case, whereas
doubling fiber damping reduced it to 0.617 BW. These quantities therefore
cannot be treated as interchangeable tuning parameters. The large pelvis
excursion remains; completing this assisted command replay is not an accepted
standing or stepping result.

The baseline compliant replay ran at 2.05 simulated seconds per wall second,
versus 3.46 for the source-informed rigid replay: **1.686 times the complete
mechanical replay cost**, not the earlier muscle-bank-only ratio of about 35.
These times include mechanics, muscle updates, contacts, guard checks, and
logging buffers; they exclude live neural computation, imports, compilation,
and NPZ saving. The native beta-change negative control is bit-identical.
A fresh-process eight-second compliant repeat reproduced qpos, force,
activation, and equilibrium-residual arrays **bit-identically**.

Plot: `screen_8s_v2/PHYSICS_COMPARISON.png`, visually checked. It compares the
right soleus force, combined normal foot load, and pelvis tilt under the same
command sequence. Individual NPZ files retain all positions, forces,
activations, and diagnostics.

## An integration difference found and isolated

Simply moving native muscle forces into `qfrc_applied` preserves instantaneous
generalized force to about 6.82e-13 but does **not** preserve native implicit
integration. Native muscle velocity derivatives contribute to the
`implicitfast` matrix; arbitrary externally supplied forces do not. This is
consistent with the [MuJoCo 2.3.7 integration equations](https://mujoco.readthedocs.io/en/2.3.7/computation.html#numerical-integration).

Three controls separate force mapping from integration:

| Control | Maximum full-qpos disagreement |
| --- | ---: |
| Isolated Euler audit, native versus force-only external replay, 0.2 s | 2.37e-15 |
| Pinned implicitfast, force-only external replay, 2 s | 0.01389 mixed qpos units |
| Pinned implicitfast, native-law force plus matching local affine velocity tangent, 2 s | 1.39e-11 mixed qpos units |

The affine control reconstructs the force and velocity tangent of the native
law through supported actuator parameters. Its small mismatch establishes a
practical way to retain implicit muscle damping when developing a custom
backend. It is presently a **native-law diagnostic**, not an implemented
implicit compliant-tendon tangent. The rigid/compliant screen above uses the
same partitioned external integration on both sides; it does not conflate
that comparison with native/external trajectory equivalence.

Evidence: `preflight_euler_v3/SUMMARY.json`,
`preflight_affine_2s_v5/SUMMARY.json`. The earlier 0.2 s implicit difference
was 0.01354, retained in `preflight_v2` and `preflight_affine_v4`. The first
attempt, `tendon_replay_v1`, stopped during setup because it incorrectly
treated selected `act` channels as all 92 activations; its manifest is retained,
and no simulation result is attributed to it.

## Current boundary and next checkpoint

The existing compliant prototype still uses approximate force–velocity
behavior and experimental fiber damping. Source Thelen static agreement does
not validate a complete dynamic Thelen/Millard/Hyfydy model. This campaign
advances it from a muscle bank and one-joint fixture to a tested full-Gait
mechanical replay, with source/bilateral/provenance checks and measured
whole-mechanics cost. It has not yet replaced the live SNS force backend.

Two isolated two-second compliant timestep-refinement controls at 1 and 0.5 ms
completed with no hyperextension or fiber bounds. Across the 17 independent
anatomical joints, RMS disagreement against the 0.5 ms control decreased from
**0.03554 degrees at 2 ms to 0.01142 degrees at 1 ms**; maximum disagreement
decreased from 0.4804 to 0.1986 degrees. Across all 92 muscle tensions, RMS
disagreement decreased from **0.7326 to 0.2395 N**. Coordinates are aligned at
physical state times after each step; applied forces are aligned at the
command times before each step. These are full mechanical-model convergence
checks, distinct from the earlier prescribed-muscle trajectory checks.
Raw data: `refinement_1ms`, `refinement_0p5ms`; summary in
`screen_8s_v2/VALIDATION.json`. Production defaults stay at 2 ms.

Next, validate dynamic muscle responses against OpenSim, and carry the
constitutive velocity tangent into an opt-in live SNS backend before any
controller retuning. Rustuna's isolated runtime exists, but no parameter
optimization was invoked without a validated physical target/objective.

The driver checks an immutable per-run manifest and skips completed cases,
allowing bounded batches to resume. Manifests record source trace, model,
parameter source, and script hashes; script snapshots are retained with the
main comparison. No manuscript, production runner, trained winner, or curated
shared-force profile was changed.

## Source-dynamic and live-SNS physics validation

`thelen_dynamic.py` separately implements the source Thelen2003 force–velocity
relation, including activation dependence and extrapolation. Explicit settings
are Af=0.3, Flen=1.8, maximum contraction velocity10 optimal fiber lengths/s,
and minimum activation0.01. Source static-curve parameters and tendon strain
atF0=0.033 are checked during loading. Additional fiber damping is zero for
these results. The earlier approximate module is unchanged.

### Independent physics checks

- Actual installed OpenSim4.4.1, 270 force/velocity cases: maximum active-force
  disagreement **5.19e-8N**.
- Three actual OpenSim prescribed-motion transients (2mm,1Hz,1s): forceRMS
  errors **0.0442/0.1820/0.1675N** at2ms, decreasing to
  **0.0221/0.0911/0.0838N** at1ms. No fiber bounds.
- Analytic length tangent:82 independent finite-difference cases, maximum
  relative error **4.21e-8**. Predicted-length velocity tangent:138 cases,
  same maximum error. Repeated evaluations never commit state. Two-sided
  smoothness is not claimed at source curve kinks.
- Matched source rigid/compliant mechanical replays both complete8s without
  hyperextension. Times: **3.255s rigid /6.658s compliant**. Elastic residual
  <=5.69e-9N. Local linearization defect reaches0.57%F0 compliant and10.7%F0
  rigid at impacts; the adapter is one linearization, not a monolithic solve.
- Full-model2s refinement against0.5ms: jointRMS disagreement decreases
  **0.02801→0.01206degrees** (2→1ms), and forceRMS **1.677→0.540N**.

The adapter predicts end-step routed length, supplies force and its velocity
tangent through affine actuator parameters, then commits fiber state once at
actual new length. Native MuJoCo activation dynamics are retained. Production
defaults stay2ms/implicitfast. Evidence: `dynamic_source_v2/VALIDATION.json`,
`source_dynamic_full_{rigid,compliant}_8s`, and
`source_dynamic_compliant_{1ms,0p5ms}`. Validation entry point:
`validate_thelen_dynamic.py`.

### Accepted live neural standing

`live_tendon_trial.py` installs the backend only in its own subprocess. It
reuses the exact three-parameter constrained syn6 candidate from the other
agent's `gait_bpa/STANDING_CHECKPOINT.json`. Neural parameters, native static
initializer, corrected pose, sensory normalization, and apparatus stay fixed.
Actual committed muscle forces feed the existing force sensors. There is no
new neural fit, MN-current cancellation, or force allocator.

| Runtime force law | Accepted | Peak tilt | Sway | Mean foot load | Up/down apparatus force |
| --- | --- | ---: | ---: | ---: | ---: |
| Native fitted reference | Yes |2.646degrees|43.85mm|0.9930BW|0.00606/0.00032BW|
| Source-rigid | Yes |2.704degrees|31.79mm|1.0267BW|0.00045/0.02823BW|
| Source-compliant | Yes |2.368degrees|32.67mm|1.0240BW|0.00072/0.02578BW|

Each requested8s neural protocol includes the inherited0.6s pinned
initialization, followed by7.4s free physics. Native numeric arrays reproduce
the other agent's saved candidate exactly. The compliant fresh process
reproduces time, qpos, force, activation, and excitation **bit-identically**.
Its knee maxima are **-1.906/-1.976degrees**, checked every physical step;
equilibrium residual <=2.64e-9N; prediction defect <=0.171%F0.

Observed whole-protocol wall times are29.70s native,33.37s source-rigid,
37.66s source-compliant (fresh repeat37.11s), including initialization and
neural/mechanical measurements. These are not kernel timings or
faster-than-real-time liveSNS claims. Lateral/rotational apparatus and soft
ground penetration remain: this is supported standing on the human Gait2392
model, not independent standing or robot gait. Inner runner force-profile
metadata describes native initialization; wrapper metadata identifies the
actual runtime law.

Compact evidence: **`SOURCE_DYNAMIC_CHECKPOINT.json`**; raw conditions:
`live_syn6_{native,rigid,compliant}_v1` and
`live_syn6_compliant_fresh_v2`.

### Hunt ground transfer and next checkpoint

The unchanged16s corrected-pose Hunt example was also run with all three
laws. Native physics reproduces accepted ground stepping. Source-rigid and
source-compliant both complete safely with5/4 knee cycles but **0/0 sustained
visual-mesh-cleared swings**. Source-compliant knee phase error is0.00814
cycles, measured assistance0.6104BW, and whole-foot clearance maxima about
-0.024/+0.115mm. All knee/force/clearance gates remain unchanged.

The source-rigid failure shows that clearance loss cannot be attributed to
tendon compliance alone. The saved recruitment was fitted for the native
surrogate. Next: backend-consistent motor/recruitment calibration, followed
by the original ground gates. This Hunt result and the accepted constrained
syn6 standing result remain separate controller families.

Dynamic force-law agreement is now measured, but complete OpenSim activation,
extreme pennation/fiber-limit behavior, contact calibration, and end-to-end
OpenSim/Millard/Hyfydy equivalence remain unestablished. No manuscript,
production source, trained winner, or curated force profile was changed.

## Backend-consistent Hunt motor calibration

Fresh native and source-compliant16 s air simulations pass both original stages
with the same saved shared controller (5/4 knee cycles, negative knee angles
throughout). In the 5–13 s actual-physics window, source-compliant stage 2 knee
excursions are 61.63/59.44 deg versus native86.81/85.27 deg, despite retained
alternation. This supports recruitment calibration before ground transfer.

The isolated `hunt_motor_calibration.py` adapter runs a bounded 12-trial Rustuna
study (30min wall budget), varying three **bilaterally tied** existing PF-to-MN
group scales: hip flexors, knee flexors (including the hamstrings), and ankle
dorsiflexors. The objective compares actual free air joint trajectories with
the accepted native example, allowing one common ±0.2 s time shift, normalized
by native excursion. It never imposes motion. Each candidate must pass both
unchanged air gates before any ground evaluation. Intrinsic/timing/reflex
weights, source muscle dimensions, activation dynamics, and apparatus remain
fixed. The existing per-muscle length midpoint/half-range derives from routed
lengths along the source IK trajectory; it is retained because replacing it
would invalidate the separately fitted Hunt B input mapping.

Raw results and restart manifest:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\hunt_motor_calibration_v1\`.
Fresh air controls and delivered-force breakdown:
`physics\calibration_air_baselines\BASELINE_DIAGNOSIS.json` under that campaign
root. These compare the same controller after forward physics, not forces at
identical prescribed kinematics.

**Diagnostic correction:** the native `SwingDiagnostics` zero-activation probe
cannot estimate passive force when an external stateful law is represented by
affine actuator bias. Previous source-backend `.mechanics.npz` passive/active
columns are unsuitable for that interpretation. Actual backend-trace force,
activation, motion, and acceptance measurements remain valid. New source runs
disable that probe and log active fiber, passive fiber, and damping terms at
the actual committed state. Their sum matches delivered tendon force to
4.92e-9 N in these air controls. Source damping beta remains zero.

### First source-compliant supported-ground passer

The shared scales **hip_flex=1.5, knee_flex=1.5, ankle_df=1.5** pass both air
stages and the unchanged ground protocol. Ground results: 5/3 sustained
whole-visual-mesh-cleared swings, peaks53.85/35.57 mm, measured contact
onsets4/4, contact phase error 0.02641, and maximum knee angles
−1.618/−1.539 deg. Measured vertical apparatus assistance is 0.68774 BW.
An exact fresh process reproduces time, complete physical state, activation,
excitation, delivered force, fiber state, and force components **bit-for-bit**.
This is supported in-place stepping with source-Thelen compliance. Basin,
standing-stage, and reduced-support validation remain separate requirements.

Compact external checkpoint:
`physics\hunt_motor_calibration_v1\HUNT_CALIBRATION_CHECKPOINT.json`.
The bounded air-only Rustuna study continues; its best tracking score is not
automatically a ground winner. The first verified ground passer is retained.

`source_force_comparison.py` compares force laws at identical recorded ground
motion and activation, reconstructing/checking source fiber state to 4.38e-9 N.
In the failed uncalibrated ground example, mean soleus passive force is
171/132 N source versus 0/0 N native; rectus femoris passive force 18.3/13.5 N
versus 6.38/4.66 N. BFSH active force during knee flexion is 116/116 N source
versus 136/133 N native. These are actual constitutive differences that motivate
recruitment calibration without changing anatomical tendon/fiber dimensions.
The query of zero activation is performed only on a separate native-law model.

One attempted rigid air comparison accidentally reused ground apparatus
settings and was excluded (`rigid_candidate_1_air1/INVALID_PROTOCOL.json`).
Its corrected original-air counterpart is saved under the `_v2` name.

The matched source-rigid law also passes both air stages and ground with the
same 1.5/1.5/1.5 gains:4/3 sustained cleared swings, actual contact4/4, contact
phase error 0.01953, assistance0.66135 BW, and maxknees−1.700/−1.553 deg.
For the compliant ground passer,7700 physics steps have no fiber-bound hits,
max equilibrium residual7.02e-9 N, and max local force-prediction defect0.00520F0.
Observed complete16 s protocol wall times are 64.07 s compliant and 56.61 s rigid
(including initialization, logging, and scoring, during concurrent jobs;
these are not an isolated speed benchmark).

The 12-trial Rustuna fit completed. Best air tracking objective 0.0143764 versus
identity0.0184482 (22.1% lower) uses hip_flex 1.4153501, knee_flex 1.5273282,
ankle_df 1.4543317. Its separate ground test passes with 5/3 cleared swings,
4/4 contact onsets, phase error 0.02637, and assistance 0.69958 BW. Its fresh
repeat is also bit-identical. Better air tracking therefore does not imply
lower ground assistance; the first passer remains the primary checkpoint.

### Completed recruitment basin

All **8/8** corners of the three tuned gains at ±1% (each gain 1.485 or 1.515)
pass both original air stages and then the unchanged ground protocol. Every
ground case has 5/3 sustained cleared swings. Measured assistance ranges
0.68317–0.69271 BW and contact phase error 0.02573–0.02650 cycles. The largest
knee extension over every saved physical state is −1.55898 deg; there are zero
fiber-bound evaluations and the largest force-balance residual is 7.49e-9 N.
The basin does not test physical-parameter, initial-state, disturbance, or
long-duration generalization. All runs are complete; no background job remains.

The next useful checkpoint is the same fixed controller on a measured-support
ladder, alongside a separate stage-3 balance check. A future Hunt B refit could
use shared bilateral source/fiber-normalized afferents: existing IK-derived
midpoints/half-ranges retain small sagittal and larger gluteal asymmetries
(`LENGTH_NORMALIZATION_AUDIT.json`). This campaign deliberately preserves that
input definition and all source anatomical dimensions.

`SOURCE_STABILITY_AUDIT.json` verifies one unchanged trained-weight set and one
hash per invoked controller/helper, runner, and force-backend file across all
53 valid calibration/validation cases. The three declared recruitment scales
are the only fitted quantities. Curated force-profile v1/v2/v3 hashes remain
unchanged. The external folder retains source snapshots and restart manifests.

## Source-normalized stance Ib and reference-envelope checkpoint

Machine: easteregg2. Complete evidence and the separate experiment record are
in `D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\hunt_ib_load_v1\`.
`CHECKPOINT.json` is the compact numerical index. The earlier 50% support
clearance result remains an engineering result: its joint trajectories violate
the saved retargeted subject01 IK envelope, especially at the ankles.

`biomechanical_envelope.py` freezes per-side hip, knee, and ankle bounds before
the next fit, intersects knee extension with the user's strict zero-degree
limit, and scores the logical 5–13 s window. These are reference-specific bounds,
not universal human ranges. Equality-aware route checks exclude per_tert from
new plantarflexor load recruitment because its ankle moment arm is dorsiflexing.
Source-normalization-only retains legacy routing; force-pool and stance-load
variants filter the new load pathway. Source muscle dimensions remain fixed.

`stance_ib.py` reuses existing excitatory pools, removes their additive RG_E
drive, and instead permits their outputs through endogenous phase, measured
whole-foot load, and descending drive. Force inputs use original source F0.
Air stages 1/2 and zero-vertical-support standing pass; both new conductances
are exactly inactive throughout those prerequisite traces. Eight bounded
Rustuna trials fit only common bilateral knee and plantarflexor load-output
conductances. No trial passes the frozen IK envelope, so support is not reduced.
The best diagnostic case improves the envelope penalty but still reaches
42.16°/35.91° ankle dorsiflexion. Its accelerated fresh repeat is bit-identical.
Regular SNS differs by at most 3.02e-13 rad in joint position and 3.32e-9 N in
force, with the same engineering outcome.

For this 1033-neuron source-compliant simulation, matched 16 s trials take
68.38 s accelerated and 233.43 s regular: 3.41× for the complete simulation.
The measured neural/controller portion takes 15.67 s and 185.02 s: 11.81×.
This includes sensory adapters and differs from a neural-only kernel benchmark.

`pt_current_chain.py` independently copies the neural state and probes source
steady forces at the actual standing pose. The live trajectory remains
bit-identical. At current scaling, 0.05-unit commands produce forward response
0.947 and lateral responses +0.902/−0.899 per unit, but backward response is
−0.415 with vertical cross-effect −0.597. Increasing all PT excitation worsens
the combined static error. This is a local transfer audit, not recovery
validation. Source-consistent motor/length-afferent calibration and backward
PT cross-coupling are separate next targets. No production or manuscript edits.

### Source-fiber and loaded hip/knee checkpoint

The current external record is
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\hunt_spindle_v1\EXPERIMENT_RECORD.md`.
Six shared ankle motor trials, local Ia removal, and source-fiber ankle encoding
do not pass the frozen anatomical envelope. The transient PF candidate passes
the air/standing prerequisites but fails ground; its fresh repeat is exact.
No support promotion or source-length change follows.

`source_force_comparison.py` now audits archived reference contact feasibility
and actual equality-aware hip/knee dynamics. The fixed-pelvis reference has
12/17 contact-free sampled poses and cannot serve as a loaded ground teacher
unchanged. A height-only correction leaves large stance slip. The unchanged
incumbent's full-force diagnostic replay is bit-identical to its earlier trace;
force reconstruction and full dynamics close within 5.42e-9 N and 5.41e-12,
respectively. Source quadriceps force-velocity depression and impact/stop
transients are material; sampled ordinary loaded states retain static source
capacity. Median passive hamstring opposition is zero, so it does not explain
the general failure. Direct OpenSim checks at nine actual full-coordinate poses
confirm the short TA/long soleus fibers are postural consequences. Original
source dimensions remain fixed. Next work needs a contact-consistent supported
reference and dynamic Hunt C recruitment targets, rather than another ankle-only
gain sweep. Exact evidence and explicit diagnostic scope are in the external
`reference_inverse_v1` folder. No production or manuscript edits.

### Native-contact/source-muscle keyframes

Current external record:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\supported_teacher_v1\EXPERIMENT_RECORD.md`.
Native inverse-contact calibration plus source recruitment now produces
verified double-support and separate single-support mechanical holds.
At declared60% upward support, original-source92 right/left holds keep the
airborne foot above10.95/11.26mm for2s, with zero measured normal force,
no knee extension, and less than1.83mm loaded-foot penetration. Fresh repeats
are bit-identical;1ms checks retain the gates. These are static muscle-driven
calibrations, not an SNS or alternating-step result.

`PLANTAR_CONTACT_AUDIT.md` records the missing distinction between distal-toe
and metatarsal-ball contact, the native soft-contact/root gate, and the important
actuation-policy difference: the source-law research variant includes all92
original muscle forces; legacy runner pruning disables six short rotators.
Explicit pruning-preserved holds and a5% activation-floor probe remain separate,
with their failures recorded. No production pruning, geometry, or source strength
was changed. Hunt C remains gated on a smooth, native-contact-consistent dynamic
reference with an explicitly selected actuation policy.

### Continuous native-compliant transition checkpoint

Current compact record:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_transition_v1\EXPERIMENT_RECORD.md`.
Both supported swing-initiation references now pass unchanged source-ROM,
zero-hyperextension, native-contact, clearance, and source-rate gates. Fresh
repeats are exact and1ms checks pass with the same10ms teacher clock. The
specific fix is shared model-scaled MTP feedback: the generic180/12 joint
gains demanded toe torques far above source capacity. The isolated teacher
uses10Nm/rad and0.225536Nms/rad at both MTP joints, without changing physical
damping, source strength, or geometry.

A continuous left-foot lift/return also passes, including recovered two-foot
loading:13.698mm peak clearance,1.928s unloaded,4.838mm/s maximum loaded slip,
1.894mm maximum penetration, and maximum knee angle-4.343deg. Its fresh
repeat is exact and1ms check passes. Actual excitation, fibers, forces, and
motion are exported in one19.25MB compressed final reference for later motor
calibration. This remains an offline source-muscle teacher, not an SNS or
alternating-gait result.

The opposite complete return fails6.733mm/s slip and2.217mm penetration.
Slower lowering alone does not solve it. A root/world-foot tracking mismatch
is now the next causal target; no bilateral HuntC fit or support reduction
has been promoted. Original-source92 actuation, the1% source floor, prior
pruning-preserved negatives, and the rejected5% floor remain explicitly
separate. Production and manuscript files were not edited.

Endpoint interpretation correction: those two-second holds are not verified
stationary equilibria. The left single-support knee drifts from -5 deg to
-9.843 deg, including -2.070 deg in the last half-second. The root/clearance
measurements remain valid under their original bounded-hold gates. A full
105-DOF native equilibrium solve, retaining soft-follower forces, and five-second
constant-control/fresh validation now precede any bilateral dynamic teacher or
Hunt C fit. See the external record's correction and
`STATIC_HOLD_TARGET_ALIGNMENT.json`. No source dimensions, contact geometry,
native stops, or production physics were changed.

### Full-native endpoint calibration

The active compact record is
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_equilibrium_v2\EXPERIMENT_RECORD.md`.
The revised solve retains every native force component and solves loaded soft
follower offsets rather than projecting away their residuals. Numerically
verified translation invariance permits analytic elimination of the two
horizontal static rig coordinates at their spring references. This adds no
live constraint and preserves the original spring/damping/contact model.

All static gates now pass for one original-source92 double-support candidate:
maximum full native residual 0.00227 N or Nm, foot-position error 31.5 micrometres,
foot-rotation error 0.707 mrad, and maximum penetration 0.964 mm. Removing the
competing exact 8% toe-load preference allowed the existing 6–10% admissible
interval to operate as intended. Acceptance limits were unchanged.

Its five-second source-compliant repeat is bit-identical in every saved field,
with 29.6 micrometres of root drift and 0.0397 deg maximum joint drift. It still
fails the stationary-endpoint gate: final static-source-force residual is
0.0681 N versus the frozen 0.05 N limit. Source activation is unchanged, fiber
bounds and knee hyperextension are absent, and no dynamic Hunt C fit is promoted.
The next bounded refinement tests native force-balance numerical precision.
`native_equilibrium_validation.py` is the reusable forward validator; new
full-native endpoints cannot enter transition training without its acceptance.

Subsequent inverse-Jacobian numerical scaling resolves that remaining hold
failure without changing physics or gates. The accepted double-support folder
is `native_equilibrium_v2\double_support_source92_jacobian_scaled_5s` under the
external physics root. Root drift is 12.1 micrometres, joint drift 0.01595 deg,
late joint speed 0.00776 deg/s, and full final static residual 0.02727 N. Fresh
saved trajectories are bit-identical; the same controls also pass a five-second
1 ms check. No knee hyperextension or fiber-bound events occur. The separate
single-support endpoint still requires qualification before a bilateral dynamic
teacher or Hunt C fit; see the external `CHECKPOINT.json` and current record.

Both single-support endpoints now qualify, including exact fresh repeats and
1 ms checks. The right-side limiting issue was numerical conditioning of the
small swing-foot inertias. Adding all 105 actual native forward accelerations
to the static objective converges in five evaluations, without changing source
physics. Right five-second drift is 0.000025 deg and final static residual
0.000198 N; left drift is 0.00103 deg and residual 0.01023 N. The separate
validator's initial audit confirms native forward/inverse force agreement near
1e-12 N or Nm.

Matched continuous lift/return examples now use these stationary endpoints
without any state reset. Both clear and unload the swing foot and recover
double support, but both fail the unchanged touchdown-slip gate: 37.82 and
17.89 mm/s versus 5 mm/s. No dynamic Hunt C fit is promoted. Existing Cartesian
teacher FK projects away loaded follower offsets and needs correction before
it can supply a physically consistent no-slip target. All curated paths,
numerical objectives, and rejected cases are in the external checkpoint.

The source mass audit finds 3.2 kg in 32 auxiliary path-point bodies beyond
the original OpenSim 75.1646 kg; original anatomical body masses are unchanged.
Those auxiliary masses remain in the captured model. Its fixed ground body
also contributes 0.1 kg to historical BW normalization. These are documented
alignment/accounting limitations, not silently changed physical parameters.

### Loaded-FK touchdown follow-up

Current record:
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_touchdown_v3\EXPERIMENT_RECORD.md`.
The optional Cartesian adapter now preserves full loaded keyframe geometry,
live follower displacement, and its measured velocity. Default-off class
replay and the complete five-second physical replay are bit-identical to the
frozen references. Independent world-foot velocity errors are below 2e-10 m/s.

The old Cartesian path requested 11.10 mm/s horizontal speed at touchdown.
A separate, validated C2 tangential/orientation phase finishes above 5 mm
planned clearance, reducing its requested low-clearance tangential speed to
0.030 mm/s while preserving endpoints and vertical motion. Actual touchdown
still fails. Native root-wrench recruitment gives the closest matched pair:
13.32/15.75 mm/s slip, with left peak clearance 9.773 mm below the 10 mm gate.
Full XYZ targeting exposes pelvis sag and an ankle-bound failure; instantaneous
mass-consistent impedance is rejected on large toe oscillation and source-rate
failure. Activation/fiber prediction remains a hypothesis to validate, not a
proven sole cause. No Hunt C fit, mass sensitivity, or gate relaxation follows.

The follow-up predictor now reproduces native activation to 1.11e-16 and
source-force replay below 2.7e-10 N on two immutable historical trajectories.
At 10 ms, causal torque forecasts improve on held force and instantaneous
command-gain baselines. Live conditional recruitment nevertheless fails both
sides after preserving the original warmup. No horizon/amplitude sweep follows.

The augmented-state diagnosis resolves that apparent contradiction: independent
native window replay is exact, but the conditional-force/instant-acceleration
Jacobian differs from the full command transition by 153%/264% in the mass
metric. Activation, fiber, filter, q, and velocity states must propagate together.
Current acceleration has no direct command feedthrough. All exact matrices,
source replay checks, and final default-off bit-identical checks are recorded
under `native_touchdown_v3\full_actuation_jacobian` and its checkpoint. No
dynamic teacher or neural fit is promoted; physical parameters remain frozen.

## Full augmented transition, conditioning, and bounded inverse

`AugmentedTransition` now restores and propagates all 486 q, velocity,
activation, fiber, and excitation-filter states together, with native contact
at every 2 ms step. Its preallocated 10 ms rollout reproduces saved q within
1e-16, and its command-to-acceleration columns match the previous independent
full-state audit exactly. No reference, gain, source dimension, force floor,
mass, or physical damping/contact parameter changes.

`native_touchdown_v3/augmented_inverse_audit/AUGMENTED_INVERSE_AUDIT.json`
records both frozen windows. At 0.5 s, directional errors of 7.97–11.31%
reject the inverse despite exact replay. At 0.7 s, half-step and independent
directional errors are approximately 2.1e-4 and 1.3e-4, respectively. The
23-row objective Jacobian has full numerical row rank but condition number
about 100,000. A single regularized correction capped at 0.02 excitation
reduces the actual nonlinear objective norm from 342.345 to 334.993; the
linear prediction is 335.005. Source ROM excess is zero. Every desired row
lies outside the independent linear correction range of this small command
box at that already failing state. This is a local 10 ms authority result,
not a claim of globally insufficient muscle strength or a successful teacher.
The complete audit uses 3,725 native steps in 6.541 s, including restoration,
finite differences, checks, and the inverse; this is not a live-loop speed.

The separate `augmented_branch_derivatives/AUGMENTED_INVERSE_AUDIT.json`
tests one-sided columns at 1e-4/5e-5 excitation and independent directions
at 1e-5. At 0.5 s, 117 of 460 native activation-step distances and 25 fiber
velocities are within 1e-8 of zero. The tiny column perturbations cumulatively
switch 400 activation branches and 368 fiber-velocity signs. Signed one-sided
column superposition still fails independent checks by 11–37%; merely
shrinking central differences or choosing each muscle's activation branch
does not establish a valid coupled derivative. At 0.7 s no branches switch,
and refined independent errors fall to 1.3–1.6e-5. These observations implicate
the coupled nonsmooth transition, without isolating activation, fiber
velocity, or native contact as its sole cause. No inverse is attempted in
this diagnostic mode. It uses 6,310 native steps in 11.440 s.

The equilibrium derivative gate remains closed. No live augmented inverse,
physical teacher, Hunt C fit, gain/horizon sweep, support promotion, or mass
sensitivity follows this checkpoint. The protected constitutive engine hash
remains a4033ef11f6cf03bf7d35e896fc4690e8500c183eaeae5dfc6131751d776ba69.

## Nonlinear shared-family load-transfer checkpoint

`stationary_group_inverse` starts from the unchanged accepted native double-
support endpoint, not the failed predictor trajectory. Static source-force
reproduction is exact. Eleven common functional-family amplitudes multiply
the saved bilateral differences toward the accepted left-stance endpoint.
These are endpoint-state command corrections, not fitted independent leg gains.
RF, HAM, and GAS remain separate; sartorius and gracilis stay unchanged.
Every muscle command correction is at most 0.02 excitation. All other commands,
source parameters, native contacts, guides, and 60% prescribed support remain.

The explicit diagnostic task transfers 1% protocol BW from right to left while
retaining loaded foot geometry. Component residuals have declared physical
scales: 1% BW normal force, 0.2 mm foot position, 5 mm/s loaded-contact slip,
5 mm root translation, and 2 degrees root rotation. The nonlinear trust-region
search evaluates actual full native augmented rollouts and reports each
component separately. It assumes no differentiable command Jacobian at the
source-law/contact branch crossings.

At the initial fixed 10 ms horizon, the best bounded correction changes foot
load by only about 0.00029 BW and fails the 1% BW task. A three-pole small-signal
fixed-MTU model derives the response time from the unchanged 50 ms engineering
filter, native activation/deactivation constant, and original source fiber
velocity tangent. The slowest selected channel is right medial gastrocnemius:
77.66 ms activation/deactivation time, 40.41 ms conservative fiber time, and
301.98 ms to 90% force response. Its predicted force response at 10 ms is only
0.092%. This derives one next fixed horizon of 302 ms; there is no horizon sweep.
These are isolated constitutive response estimates, not whole-model settling
times. Physics remains at 2 ms and the existing neural/control clocks do not
change.

At 302 ms, actual right/left normals are 0.170806/0.227965 BW versus targets
0.170752/0.228020 BW. The maximum error is 0.0000549 BW. Root translation is
0.0512 mm and rotation 0.0350 degrees. The bounded solve uses 160 evaluations
across the two declared horizons and 13,074 native steps in 19.47 s, excluding
setup; this is inverse-audit throughput, not a live-controller speed.

Fresh saved-command replay reproduces final q, velocity, activation, and fiber
state bit-identically. All-window maximum loaded slip is 0.433 mm/s, minimum
visual-mesh gap -0.989 mm, source speed/acceleration ratios 0.00366/0.00968,
and both knees remain below -15.95 degrees. The same saved correction at 1 ms,
without refitting, retains the task and window gates: force error below
0.0000743 BW, slip below 0.438 mm/s, and mesh gap above -0.989 mm. Maximum
final q difference is 1.44e-5 native units.

Evidence is in `native_touchdown_v3/stationary_group_inverse`,
`stationary_group_inverse_fresh`, and `stationary_group_inverse_1ms`.
The exact original inverse source version has SHA256
6a765558fddd82894427330c8acd7b0f8c2366d618ef24c7812a9528b410d488 and is saved
once in external `code_snapshots`, alongside the reusable validation version.
This verifies a small load-transfer primitive only. It is not complete swing
unloading, a dynamic stepping teacher, a neural controller, or a new stationary
endpoint. No Hunt C fit or live-controller promotion follows.

## Continuous bilateral 1/3/5% BW continuation

Both directions now pass continuous 1% then 3% BW transfer with the same eleven
functional families, fixed 302 ms segments, component weights, and local
0.02-excitation trust regions. Direction-specific saved single-support
activation differences specify phase-state references; there are no learned
independent left/right gains. The first right-unload 1% command and its final
q, velocity, activation, and fiber states remain exactly the preserved result.
All physical/latent boundary-state errors are zero. Native contact, excitation
filtering, activation, fiber state, velocity, and loaded followers continue
through each accepted segment without a physical-state reset.

At the 3% targets, maximum normal-force errors are 0.000455 BW for right
unloading and 0.000501 BW for left unloading, within the unchanged 0.002 BW
diagnostic tolerance. Through the accepted prefixes, maximum loaded slip is
0.702 mm/s, minimum visual-mesh gap -1.065 mm, and both knees remain flexed.
Fresh saved-command replays are bit-identical for every saved q, velocity,
acceleration, activation, fiber, filtered excitation, command, source force,
and normal load. At 1 ms, both saved prefixes retain their gates; force errors
remain below 0.000515 BW and loaded slip below 0.715 mm/s.

Both searches stop at their first failed 5% target. Right unloading achieves
normal errors +0.002513/-0.002399 BW; left unloading gives
-0.003888/+0.003816 BW. All retained best-path slip, penetration, source-rate,
root, and knee checks pass, but the load task fails. Their absolute calcaneus-
body position changes reach the declared 0.2 mm bound. Thirty-nine of 107
right-side and 28 of 104 left-side recorded 5% trials terminate specifically
at that position constraint. The right best also reaches the hip-abductor
local command limit; the left best reaches no upper group limit. This bounded
search does not establish global force or task infeasibility, and no convergence
certificate was retained from its 100-call trust-region budget.

Read-only loaded FK identifies the limiting coordinate as lateral displacement
of the left calcaneus body in both directions: +0.199989 mm or -0.199844 mm.
The bodies also rotate by small fractions of a degree. A body-origin position
hold therefore constrains foot rocking and compliant geometry as well as
translation; it is distinct from actual loaded-contact slip, which stays below
1.070 mm/s on the retained failed paths. The 0.2 mm value is this diagnostic
reference hold, not a physiological slip definition. Neither constraint is
relaxed here. The next decision concerns a contact-consistent stance reference,
not another gain or horizon sweep.

Source snapshots for the searches have SHA256
f677cab37ee8734d7016a4d3d5609cb893f5d325506884fe8ef74c046efe95a1.
The searches take 59.03/76.42 s for 33,557/48,264 aggregate native steps;
these are inverse-search workloads, not live-controller timings. Only one
compressed final continuous trace per direction is saved; other evaluations
retain compact metrics. Earlier search metrics named `fiber_speed_l0_s` measure
the tendon-aligned fiber-state speed per optimal fiber length. The reusable
code now labels that quantity explicitly and additionally reports the actual
pennation-corrected fiber-length speed; neither diagnostic is a fitted limit.

Evidence: `native_touchdown_v3/continuous_transfer_{right,left}_unload`, their
`geometry_audit` subfolders, and `continuous_transfer_{right,left}_{fresh,1ms}`.
No accepted stationary endpoint, source physics, support, or original 1%
checkpoint changes. No actual swing, stepping teacher, live controller, or
Hunt C fit is promoted.

## Explicit contact-consistent rolling-stance task

`--rolling-stance` registers a distinct default-off task. It replaces the
calcaneus-body 0.2 mm position hold and its soft cost with actual loaded-contact
relative tangential slip, unilateral normal force, and the existing native
friction cone. The old position-hold searches and their failed 5% results are
unchanged. Source joint/rate envelopes, knee protection, mesh penetration,
root bounds, force laws, support, and contact geometry are unchanged. Calcn-body
motion remains a measured diagnostic.

Foot posture stays inside the existing full-native equilibrium teacher's
subtalar/MTP +/-10 degree envelope, with the frozen six-axis ankle/hip/knee
reference bounds retained. The original OpenSim subtalar/MTP mathematical
ranges are +/-90 degrees and unclamped; they are not interpreted as
physiological ranges. No physical joint limit is enlarged. The new stance
posture bound comes from the existing inverse-teacher contract, not a new
claim about anatomical RoM.

Before either new search, forty saved contact fixtures pass independent
material-point FK/Jacobian velocity checks (3.37e-12 m/s), contact virtual
work (7.03e-17 W), and native contact generalized-force projection
(5.68e-14 maximum coordinate-force error). Native friction-cone violation
is zero. Both historical default-off accepted prefixes remain bit-identical
for every saved state and force. The contact-coordinate and friction-cone
definitions follow the [MuJoCo 2.3.7 computation documentation](https://mujoco.readthedocs.io/en/2.3.7/computation.html#contact).

Both new directions pass continuous 1/3/5% BW transfer using the same eleven
functional-family policy and fixed 302 ms segments, without physical/latent
boundary resets. At 5%, right-unload normal errors are +0.000223/-0.000048 BW;
left-unload errors are -0.000379/+0.000227 BW. Maximum loaded slip is
1.146 mm/s, minimum visual-mesh gap -1.167 mm, and knees stay flexed.
Native cone violation is zero, maximum cone utilization 0.134, and all active
contact normals stay positive. Maximum subtalar/MTP posture is 2.075 degrees.
Root translation stays below 1.21 mm and rotation below 0.160 degrees.

Actual calcn-body movement reaches 0.242/0.216 mm, modestly beyond the old
diagnostic position hold while satisfying the new contact constraints. This
is evidence that body-origin immobilization can block an otherwise accepted
rolling/contact task, not grounds to relabel the old failures.

Fresh saved-command replays reproduce every saved q, velocity, acceleration,
activation, fiber state, filtered excitation, command, source force, and normal
load exactly. Both full sequences retain all declared gates at 1 ms without
refitting. Terminal force errors remain below 0.000416 BW and slip below
1.17 mm/s. Pointwise detailed dynamics are not fully converged: right MTP
acceleration near 0.806 s is 3.371 rad/s2 at 2 ms versus -7.312 rad/s2 at
1 ms. The maximum toe-speed difference is 0.01597 rad/s; the largest position
difference is 0.0001194 rad at the left knee (about 0.00684 degrees).
Right normal-load transient difference peaks at 0.001118 BW. This localizes
a toe/contact resolution sensitivity; no damping, inertia, source curve, or
geometry adjustment is inferred or made.

Search-source SHA256 is
c7ce37745af7f2d2e6077f274f85aaa901f2de5de018d0bd85171254caba6692.
The bounded searches take 69.44/92.06 s for 34,579/49,981 aggregate native
steps. These are inverse-search workloads, not live-controller speed tests.
One compressed final trace per direction is preserved; fresh/refinement checks
save summaries only. Results live in `rolling_transfer_{right,left}_unload`,
`rolling_transfer_{right,left}_{fresh,1ms}`, and
`rolling_transfer_right_1ms_localized`, under external `native_touchdown_v3`.
Prechecks are in `rolling_stance_validation` and `rolling_defaultoff_{right,left}`.

This is a verified 5% BW load-transfer task under 60% prescribed support and
the retained guide. Both feet remain loaded: normal force on the unloading
foot is still about 0.131/0.168 protocol BW at the two 5% endpoints. No swing
clearance, complete supported-step teacher, neural controller, or Hunt C fit
is promoted. Old stationary endpoints and the protected constitutive engine
remain unchanged.

## Rolling-stance continuation to 10% BW

The accepted 1/3/5% rolling prefixes are frozen and extended through 7.5%
then 10% protocol BW in both directions. Saved q, velocity, activation, fiber,
filtered excitation, and source forces at every old endpoint match exactly.
All accepted boundary-state errors are zero. The same eleven functional-family
policy, 302 ms segments, 0.02-excitation local trust regions, source physics,
support, and rolling-contact constraints remain. No horizon or gain sweep runs.

Both 7.5% and 10% milestones pass. At 10%, right-unload normals are
0.081299/0.317646 BW with maximum target error 0.000548 BW; left-unload
normals are 0.280283/0.118043 BW with maximum error 0.000469 BW. Across
the full prefixes, loaded slip stays below 1.876 mm/s, visual-mesh gap above
-1.370 mm, root translation below 2.146 mm, and rotation below 0.334 degrees.
Subtalar/MTP posture stays below 3.865 degrees within the unchanged teacher
envelope. Knees remain flexed, with maximum extension below -12.37 degrees.
Native friction-cone violation is zero and normal forces stay unilateral.

Fresh 2 ms replays match every saved state, command, source force, and normal
load exactly. Unchanged-command 1 ms replays pass all milestones and physical
gates: terminal target error below 0.000600 BW, loaded slip below 1.904 mm/s,
and mesh gap above -1.371 mm. This does not establish pointwise detailed
convergence. The right MTP acceleration difference remains 10.683 rad/s2
near 0.806 s. During left unloading, the left MTP near 1.050 s differs by
22.216 rad/s2 (-15.729 at 2 ms versus +6.487 at 1 ms). These localized toe
acceleration sensitivities are preserved; no damping, inertia, or source-force
adjustment is inferred or applied.

Search-source SHA256:
e6450eee3fa18bffbf3de8adac8202db33f273492923fa19d742f8489af9db17.
The two bounded extensions use 32,918/35,787 aggregate native steps in
57.56/69.47 s, excluding setup. These are inverse-search workloads, not
live-controller speed comparisons. One compressed final continuous prefix
per direction is retained; exploration and validation save compact summaries.

Evidence: external `native_touchdown_v3/rolling_transfer10_{right,left}_unload`
and `rolling_transfer10_{right,left}_{fresh,1ms}`. The initial fit reports
retain inherited base-stage wording in `task_policy`; their exact `stages`
arrays and `extended_from` fields specify the 7.5/10% continuation. The writer
now uses declared-stage wording and saves the requested stage list. Raw
search reports are preserved.

The unloading feet still carry 0.0813/0.1180 protocol BW. These milestones
remain double support under the retained guide and 60% prescribed support.
No single-support unloading, swing clearance, complete stepping teacher,
neural fit, or Hunt C promotion follows. Historical body-hold failures and
all stationary/source-physics checkpoints remain unchanged.

## 12.5% transfer and first 15% failures

The exact accepted rolling 10% prefixes are extended toward 12.5% then 15%
BW in both directions, with unchanged source parameters, eleven-family basis,
302 ms segments, local command bounds, and physical gates. Both directions
pass the 12.5% dynamic load-transfer milestone. Old-prefix state/force errors
and accepted boundary-state errors are zero. Donor-foot normals at 12.5% are
0.056380 BW for right unloading and 0.092584 BW for left unloading.

Both stop at the first 15% failure, with no further search or bound changes.
Right unloading completes the segment but misses the force task: donor normal
0.035374 BW versus target 0.030752 BW, with maximum force error 0.004623 BW
above the 0.002 BW tolerance. All recorded right-side physical trials remain
inside the contact/mesh/root/source-rate/knee gates. Hip-extensor and
knee-extensor family corrections reach their allowed upper bounds.

The right knee-extensor endpoint-derived basis cap is 0.0175923 excitation;
the maximum feasible command correction along that same direction before a
muscle reaches its native command floor is 0.0189329. Donor-right vasti
commands are already near 0.01005, while stance-left vasti are approximately
0.183/0.191/0.252. This identifies a constrained recruitment-direction limit,
not global muscle-strength insufficiency or a reason to change Fmax.

Left unloading stops physically at 1.998 s, before the planned 2.114 s end:
the vertical pelvis_ty coordinate drops 5.0065 mm, crossing the unchanged
5 mm root bound. All 94 recorded left-side 15% trials meet that root guard.
Critically, the first trial uses exactly the previous accepted 12.5% commands
unchanged, yet crosses the same vertical bound after 206 ms (at 2.018 s,
height drop 5.0152 mm). Thus the 12.5% result is a dynamic transfer milestone,
not a stationary equilibrium. Continued stance drift already exists without
requesting a larger load transfer. Stance-right vasti in the selected failed
attempt deliver about 0.030 activation; the loaded recruitment path needs
stabilization before progression. No single limiting muscle or global
infeasibility claim follows from the bounded search.

Read-only authority audits reconstruct endpoint native normal loads to
1.2e-16 BW. The right failed donor load is on the calcn body, with zero toe-
body force; calcn-body ownership is not a heel-only anatomical classification.
The left stopped attempt has 48.89 N on calcn_l and 6.48 N on toes_l.
Native friction cones and unilateral contact remain satisfied. No contact
geometry, damping, inertia, muscle dimension, or anatomical bound changes.

Fresh 2 ms replays of accepted prefixes through 12.5% match every saved state,
command, source force, and normal load exactly. Only those accepted prefixes
receive 1 ms checks; both retain all gates. Terminal force error stays below
0.000652 BW, loaded slip below 2.410 mm/s, and mesh gap above -1.464 mm.
Detailed MTP acceleration remains resolution-sensitive, as previously
localized; no pointwise convergence claim is added.

Actual unloading is scored separately from the target percentage. At a
0.02 BW force threshold sustained for 80 ms, both accepted prefixes and
retained failed attempts have zero qualifying duration. Actual whole-foot
mesh clearance also fails the independent 1 mm/80 ms and 10 mm peak checks.
The fit reports score the retained final attempt including its failed last
segment; the fresh/refinement reports score accepted prefixes only. Neither
kind of low-force/clearance event is inferred from a target label.

Search-source SHA256:
7add4265b1fe60f6d6372b6b3db7a37629532bfd4c1e65efa31ae4c989ba8b70.
Search workloads are 31,559/27,629 aggregate native steps in 52.73/55.03 s.
One compressed final attempt per direction is retained; exploration saves
compact metrics. Evidence is external `rolling_transfer15_{right,left}_unload`,
their `authority_audit` subfolders, and `rolling_transfer15_{right,left}_{fresh,1ms}`
under `native_touchdown_v3`. `--audit-transfer-authority` reuses the saved
data without another forward search.

Both feet remain loaded under the guide and 60% prescribed support. No 15%
pass, actual unloading, swing, complete motor primitive, or Hunt C promotion
is claimed. The next causal task is loaded stance stabilization and functional
recruitment authority, preserving the height and contact gates.


## Full-native stationary12.5% load-split references
Both full105/source92 references at60% prescribed support pass static residuals,5-second holds,exact fresh repeats,and unchanged1-ms repeats. Maximum1-ms root drift is2.32nm/right and3.68um/left; final raw residuals1.75e-5/.00656N or Nm. Original DS default-off arrays remain exact. No source/engine/contact/bounds change.
The first full92 shared-C2 command interpolation fails: right root5mm at3.176s then left-ankle ROM at5.284s; left source acceleration at2.420s,slip14.35mm/s,and root5mm at3.836s. No Hunt C promotion. The proposed46-family bilateral basis is unoptimized.
External evidence: D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_stationary_loadsplit_v1\EXPERIMENT_RECORD.md, CHECKPOINT.json, STATIONARY_STANCE_BASIS.npz, and side-specific STATIONARY_MOTION_AUDIT.json.


## Small stationary manifold and augmented motion
External native_manifold_v1 now contains two accepted 1% stationary intermediates: 5-s/fresh/1-ms passes, exact fresh arrays. Full q/v/activation/fiber/filter shooting retained all physics gates, used 23 paired rollouts/168,176 native steps, and accepted no motion. A terminal-window sensitivity audit isolates inconsistent left velocity/force blocks; a bounded nonlinear candidate reduces the normalized pose/fiber residual by 20.1% but still fails stationarity. The six global torque singular modes omit ankle/foot directions, so full46-family authority remains untested. Next use a full or task-normalized shared functional basis; no repeated six-mode tuning or Hunt C promotion. Evidence: D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_manifold_v1\EXPERIMENT_RECORD.md and CHECKPOINT.json.


## Full46 one-millisecond derivative audit
53 paired full-family probes,775,072 native steps,and exact restore0. Full task-row coverage is1.0; no distal truncation. Right objective-block derivatives pass; left mean-state derivatives are locally consistent but force/velocity blocks remain non-smooth. Shadow-off/on full-state-force hashes match and predicted-force replay error is0 N; only right gluteus minimus1 has5 predicted/corrected FV disagreements near the terminal MTP spike train. Toe contact and native constraint states stay unchanged. All physics/source parameters and A403 constitutive backend remain unchanged. No fit,Hunt C,or load progression. External evidence: D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_manifold_v1\full46_1ms_audit_v1\EXPERIMENT_RECORD.md and CHECKPOINT.json.

## Source-consistent compliant integration causal test
The isolated `consistent_source_step.py` whole-body Newton adapter closes source tension below 1e-7 N with unchanged source curves, native physics, and bilateral parameters. Five frozen contexts pass directional finite differences and discrete work checks. Four matched 7.312-second replays retain the MTP peaks: original/consistent right peaks are 17.55760/17.55731 rad/s2 at 1 ms and 17.51299/17.51318 at 0.5 ms. Controls and activations match exactly; terminal stationarity still fails. The numerical correction does not demonstrate movement improvement. Full-105 force decomposition closes to 2.14e-12 and shows large canceling contact, equality, passive, and muscle terms; individual magnitudes do not establish causality. Original A403 remains byte-identical, and the optional backend constructor defaults to the original implementation. No production adoption or Hunt C promotion.
External record and compact checkpoint: D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_manifold_v1\consistent_integration_v1\EXPERIMENT_RECORD.md and CHECKPOINT.json. Exact executed source snapshots are in the external physics `code_snapshots` folder.

## Anatomical task-pattern reduction rejected
Two exact paired A403 baselines completed, with no new controller or accepted transition. The ten desired anatomical task directions failed actual distal row-projection coverage despite explicit ankle/foot names. The owned search was stopped during its first nonzero proposal, and a default-off preflight now rejects that insufficient reduction before any rollout. No restarted campaign or further simulation during review. Exact best zero correction, original baseline hashes, measured coverage, and the correction are preserved in D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\native_manifold_v1\task_pattern_1ms_v1\CHECKPOINT.json and EXPERIMENT_RECORD.md. Historical source-compliant 50%-support neural stepping remains an engineering example that fails the later saved-IK ankle envelope; stationary source92 references remain mechanical benchmarks.
