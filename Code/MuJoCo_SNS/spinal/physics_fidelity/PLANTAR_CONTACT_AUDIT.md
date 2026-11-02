# Gait2392 contact and source-actuation audit

Machine: **easteregg2**. Numerical evidence is external at
`D:\Users\Ben Bolen\NeuroCampaign_20261031\physics\supported_teacher_v1\`.
These are mechanics/reference-preparation experiments, not SNS curriculum passes.

## Current model and teacher mismatch

The source `gait2392_simbody.osim` used for conversion contains no contact
spheres, contact meshes, Hunt–Crossley forces, or coordinate-limit forces.
`hunt_training/prepare_targets.py:48` applies measured whole-foot GRF to
`calcn_r/l`. It does not load the articulated toe segment separately. The
official [OpenSim example](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53088704/Simulation-Based+Design+to+Reduce+Metabolic+Cost)
likewise describes external loads on the calcaneus segments.

The converted model has explicit mesh/ground collision pairs.
`runner.py:431` retains six foot/ground pairs, for talus, calcaneus, and toes.
At the audited standing poses, each foot actually contacts at a heel point on
the calcaneus mesh and a distal-toe point. There is no independently calibrated
plantar metatarsal-head/ball-of-foot contact patch. A nominal 50/50 split between
these two contacts is therefore **not** a generic heel/forefoot pressure split.

At one failed reference pose, original-source flexor digitorum/hallucis capacity
at MTP = 0 was 3.105 N m, with moment arms 6.47/5.79 mm. The distal contact was
57.49 mm beyond the MTP origin. That limits its load to approximately 54 N at
the tested pose; it does not support assigning half the stance load to the toes.
An 8% distal-toe load target was selected with capacity margin for an explicitly
heel-dominant engineering calibration. It is not measured human pressure data.

Rigidly fixing stance-foot height also conflicts with compliant contact during
load transfer. Tangential no-slip and native-law vertical/pitch deformation
must be distinguished. Optimizing arbitrary contact forces is only a necessary
feasibility check; the native inverse contact law and forward dynamics are
separate gates.

## Completed native-physics checks

`supported_step_teacher.py --native-stance` solves actual native inverse
contact equilibrium, foot microdeformation, pelvis pose, and a safe knee target.
It does not change geometry, original springs, source dimensions, or strength.
`source_force_comparison.py --native-stance-recruitment` then fits original-source
Thelen forces and runs an unclamped, constant-activation forward hold.

- Double support at 50% nominal upward support: 2 s compliant hold passes;
  pelvis drift 2.011 mm, rotation 0.231 degrees. The fresh repeat is bit-identical.
  At 1 ms, final pelvis translation differs by 23.6 micrometers and rotation by
  0.00323 degrees. The production timestep remains 2 ms.
- Right and left single-support keyframes at **60%** nominal upward support:
  both original-source 92-actuator holds pass for 2 s. The airborne feet remain
  at least 10.95/11.26 mm above the floor, with zero measured contact force.
  Both fresh repeats are bit-identical; both 1 ms repeats retain the clearance
  and knee gates. These are separate static holds, not alternating steps.
- The tested 50%-support right single-support keyframe fails the declared
  2 mm penetration limit at 2.176 mm. This is not a global infeasibility proof.
- Matched 50%-support double-support control takes 1.915 s compliant versus
  1.077 s rigid for 2 s physical time. This is physical stepping plus force
  logging, without SNS; the rigid case uses the same compliant-calibrated
  activations. It is not a complete MuJoCo–SNS speed benchmark.

## Actuation policy must be explicit

The source-law prototype reads all 92 original muscle forces. That includes
quad_fem, gem, and peri bilaterally, which `runner.py:287` intentionally disables
in the legacy model. Thus historical source-law experiments are an **original
source 92-actuator variant**, not an unchanged legacy-pruned model.

The isolated `--preserve-pruned` option explicitly masks those six force and
force-tangent outputs, retaining positive source constitutive parameters. Its
double-support hold passes, but right single support loses the airborne foot;
left single support misses 10 mm clearance slightly. A shared 5% activation-floor
probe fails the bounded inverse fit criterion; it is not a demonstrated remedy.
No pruning policy was changed in production. A feasible inverse activation
vector is not unique and is not an EMG prediction or a fitted neural controller.

## Literature and a possible later contact variant

[Falisse, Afschrift, and De Groote (2022)](https://doi.org/10.1371/journal.pone.0256311)
show that toe mechanics affect stance-knee mechanics. Their model uses six
contacts per foot and unactuated toes with added stiffness/damping; their
25 N m/rad and 2 N m s/rad values are not directly transferable to this model
with active toe muscles. Our current MTP surrogate is 10 N m/rad and
0.05 N m s/rad, added in `runner.py:506`.

[D'Hondt et al. (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11218950/)
separately represent heel, arch, metatarsal-ball, and toe contacts and calibrate
vertical placement using static equilibrium. This supports auditing an isolated
plantar contact variant attached to the correct foot segments. It does not
justify copying contact locations or stiffness without calibration.

No contact pad was added. A future variant needs anatomical placement,
static native-contact calibration, measured contact-force/CoP checks, original
versus modified geometry comparison, and matched controller/timestep repeats.
Hunt C remains gated on a smooth, full-root/native-contact-consistent dynamic
teacher; static keyframes and free-contact inverse fits do not satisfy it.
