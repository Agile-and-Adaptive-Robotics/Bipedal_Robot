# Stage 2: wiring, encoder parity, and causal isolation

Machine: easteregg2. This record concerns the opt-in biw2lvar-Hunt-PT branch.
Production s3k, w2lvar, and syn6 builders/winners are unchanged.

## What the failed test actually does

The full afferented-air test has no physical knee cycles with meaningful
excursion: late knee ranges are about 9.4/9.7 degrees, near the extension stop.
Its RGs still oscillate. Counting neural bursts or small coordinate ripples as
steps would therefore hide the failure. Ankle ordering also crosses.

Fresh-process four-way isolation, all other model/controller settings fixed:

| Timing feedback into RG | Motor reflex feedback | Physical knee cycles R/L | Phase error | Crossed ankle fraction | Original movement/frontal screen |
|---|---|---|---:|---:|---|
| Off | Off | 4/3 | 0.088 | 0 | Pass, but deafferented control, not a stage-2 success |
| On | Off | 7/7 | 0.333 | 0.313 | Fail |
| Off | On | 3/0 | 0.500 | 0 | Fail |
| On | On | 0/0 | 0.500 | 0.233 | Fail |

Thus motor reflexes cause an independent extension-lock problem, while timing
feedback independently disrupts bilateral phase and frontal separation. Neither
the loss of intrinsic rhythm nor one missing contact input explains everything.

## Verified encoder mismatch

The physical runner supplies Ib current as 0.5 times the phase-gated force
signal. The Hunt offline fit supplied the unhalved signal. The offline encoder
now includes that factor and derives modulation from the actual DRIVE neuron
voltage, matching the runner rather than assuming a constant drive fraction.

The original weights remain in runs/fitted. Only Hunt D was refitted into
runs/encoder_fixed_v1, with intrinsic, timing, and feedforward weights frozen.
Expected-afferent activation RMSE changes from 0.127291 to 0.127135. All four
initial corrected-encoder physical probes still fail. The mismatch is a real
implementation error but is not sufficient to explain or solve stage 2.

## Connection rules reviewed

Sources read include:

- SADb_audit/knowledge_base/pathways/ia-reciprocal-inhibition.md;
  type-ii-excitatory.md; ii-inhibitory.md; ib-excitatory.md;
  ib-inhibition.md; ia-or-ii-stance-to-swing.md;
  ib-contralateral-inhibition.md.
- Code/MuJoCo_SNS/spinal/WIRING_RULINGS_20261003.md.
- Code/MuJoCo_SNS/spinal/campaigns/20260930/lit_rules_conflict_audit.md.
- Neuromechanical_Models/Mujoco_SNS_models/Circuit_rules_CONNECTOME_md__connectome.json
  and the W2L builder/compiled graph used by the current branch.

These include curated interpretations and older translation audits, not just
primary experimental results. Their provenance must not be upgraded silently.

The apparent sign conflicts need endpoint, pathway, and state distinctions.
Ia homonymous excitation and IaIN-mediated antagonist inhibition are different
connections. Ib autogenic inhibition and locomotor extensor excitation can be
state-dependent pathways. Group-II excitatory and inhibitory relays have
different targets/context; deleting every II pathway or universally making it
inhibitory is not justified. A collapsed 'Ia -> MN' or 'Ib -> MN' family table
can report a conflict after erasing exactly those distinctions. Numerical gains
from different circuit/unit conventions are not interchangeable.

## Actual simplifications in this branch

The Hunt fit retained Ia homonymous excitation, Ib autogenic inhibition, and
Ib group excitation as fitted motor-reflex candidates. It silenced incoming
IaIN, IIX, IIIN, and Renshaw connections to MNs. The upstream cells still exist;
their presence in a neuron census does not mean those output pathways operate.
In particular, reciprocal antagonist inhibition was not active.

The initial Hunt build also deliberately disabled interleg wiring. Opposing
finite initialization plus shared mechanics is not a substitute for a built
commissural pathway, particularly after each leg receives differently fitted
timing feedback. The network visualization makes this absence explicit.

## Isolated corrective probes

Added default-off options for IaIN-to-antagonist-MN recruitment and a weak
disynaptic flexor-to-contralateral-flexor inhibitory path. The latter is a
C1-like functional motif, not a full Shevtsova/V0/V3 replication. Intrinsic
NaP parameters remain fixed. The original branch is unchanged when options
are zero. Group-II and Renshaw outputs have not been blindly restored.

Twelve initial wiring/timing/output-gain probes all failed the full stage-2
gate. Several recover bilateral physical cycles, but phase errors and crossing
remain. Subsequent bounded commissural/axis-feedback probes are recorded in
runs/stage2_coordination_v1. Axis feedback and a hip-abduction setpoint are
explicit engineering proprioceptive inputs; they must not be described as a
learned PF solution or runtime inverse-dynamics allocator.

No new downstream stage or ground-stepping trial is authorized by a failed
stage-2 screen. Candidates require fresh replay and perturbation checks before
progression, with all nondefault parameters and weight-file hashes retained.

## Intended training-stage interpretation

Ben clarified that the five stages were intended to emulate Hunt's four
functional training phases. Record that relationship explicitly: intrinsic
rhythm, afferent timing, standing/postural motor recruitment, loaded gait motor
patterns, and sensory/contact refinement. Hunt's four-step fitting order is
the organizing rationale; the five environment tests expose transfer failures.
Do not present the original intent as five unrelated validation environments,
or claim that a successful teacher replay has completed the corresponding
closed-body stage. Stage 4 still requires actual ground-supported stepping.

## Subsequent correction and limitation

The first passing candidate retains full timing/reflex inputs, has 0.1 times
the inherited IaIN-to-MN gain, 0.3 commissural conductance, 1.5 axis-feedback
gain, 0.08 rad native hip-abduction goal, knee-extensor output multiplier 0.3,
and knee-flexor multiplier 2. It passes fresh replay and four independent
one-percent perturbations of these configuration gains. Its 1.844 s physical
period remains slower than the reference 1.233 s; this is a movement-screen
repair, not a completed nominal-cadence training claim.

Removing reciprocal inhibition, axis feedback, or the motor-gain corrections
fails the screen. Removing commissural inhibition still passes for this case:
the evidence does not support calling that connection necessary here.

Standing passes with this revised candidate. Actual ground transfer remains
unsuccessful; the left leg generally remains loaded and does not alternate.
Tonic contact, contact-event, lateral-restraint, height, and finite initialization
probes are retained as failures rather than promoted to walking results.
