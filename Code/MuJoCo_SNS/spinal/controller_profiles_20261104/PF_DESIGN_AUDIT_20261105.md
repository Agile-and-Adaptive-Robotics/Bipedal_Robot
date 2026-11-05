# PF timing failure and diagram-routing correction

Machine: easteregg2. Session 1, 2026-11-05.

Ben identified two separate failures: hidden synapse paths in the drawings and
an inadequately verified PF design. Both were real. Prior C/D output-fit passes
did not establish correct PF timing and must not be described that way.

## Diagram repair

`diagram_routes.py` now routes all synapses orthogonally outside neuron boxes.
It checks every segment against every neuron interior before rendering. Both
the overview and the expanded flexor subnetwork have zero neuron-interior
intersections in `ROUTE_VALIDATION.json`. Crossings between synapses remain
possible and do not indicate a junction. The compiled circuit and recorded
neural traces are unchanged by this drawing repair.

Updated PNG/PDF/SVG assets retain their existing paths under:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\circuit_trace_20261105`

## Actual timing audit

The audit measures actual PF synaptic transmission `clip(V/5mV,0,1)` on the
original calibration clock. It reports onset/offset lists, duty, peak phase,
fundamental lag, excursion, activation errors, pairwise correlation, and singular
values. No individual PF channel is phase-shifted to improve the comparison.
The active threshold is 0.2. The provisional timing check requires absolute
fundamental lag <=30 degrees, duty error <=0.15, matching burst counts, and
excursion >=0.2. These are explicitly engineering checks, not biological laws.

Syn6 targets come from its source NMF primitive time profiles, not from a demand
that physical air stepping match grounded joint motion. The composed syn6
profile copies the right PF circuit/readout onto the left. Its left target
therefore uses the right primitive with the measured contralateral loading
phase offset. The independently factorized original left NMF columns are not
automatically the same channel identities. The first diagnostic input-fit run
used those original left columns; it is retained separately and superseded by
the mirrored-target experiment below. The right-side failure is present under
either interpretation.

| Finding | Measured result |
| --- | --- |
| Syn6 actual PF components needed for 99% centered variance | 3 per side |
| Intended syn6 primitive components needed for 99% variance | 6 per side |
| S5 duty, right / left | 1.000 / 1.000 |
| Intended mirrored S5 duty, right / left | approximately 0.325 / 0.330 |
| Right S1/S2/S4 fundamental lag versus target | -84.1 / -100.1 / -61.2 degrees |
| Syn6 channels passing the combined provisional timing check | 0/12 |
| s3k actual PF components needed for 99% centered variance | 2 per side |

Rank here is an SVD variance description, not a neuron count or proof of minimal
architecture. Nevertheless, it exposes the loss of intended independent timing
in syn6. Its outputs largely collapse into common E/F profiles, with S5 held
above activation threshold throughout the cycle. S3k's E1/E2 and F1/F2 overlap
strongly. Numerical intended onset/duration profiles were never declared for
s3k, so its phase labels alone cannot certify timing correctness.

## Code-level causes and missed checks

1. Syn6 reduces its intended primitive profiles to an E/F family assignment and
   a membrane-time-constant ordering. That is not a fit of their actual onset,
   duration, offset, and overlap. Shared family inhibition pushes cells toward
   common windows. S5 receives both complementary RG drives and is excluded
   from the normal PF inhibitory lamination, explaining its continuous drive.
2. In s3k `_build_pf`, `PF_SHAPE` is unpacked as `tau_m, _`. The second value,
   documented in `params.py` as an adaptation multiplier, is unused there.
   Different leaky membrane constants alone did not realize the claimed
   fast/short versus slow/long pattern windows in the measured circuit.
3. The recent timing fit checked period and bilateral phase; output fitting
   then targeted MN signals. It lacked an intervening PF timing gate and the
   required explicit phase-transition checks. That allowed MN/readout fitting
   to compensate for upstream errors rather than demonstrating a correct PF.

These are design/verification failures. They are not explained by damping,
tendon physics, or a plotting artifact.

## Bounded repair attempt

Two eight-evaluation experiments fit existing PF-input conductances using
discrete LIF current balance, with original RG parameters and MN weights frozen.
All gain multipliers are positive, bounded, and shared bilaterally. No target
voltage is clamped and no body physics is integrated. The first diagnostic used
independently factorized left targets; the second uses the actual mirrored
target contract.

The mirrored-target loss improved from 2.51171 to 2.15643, with an exact
fresh-process objective repeat, but **0/12 channels passed**. This is an
unsuccessful bounded gain-only repair, not a new accepted controller. The result
does not prove that all possible existing-edge gain fits must fail. It does show
that this attempt did not solve the design failure, and no body-stage advancement
or new walking result follows from it.

Experiments:

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_PF_input_repair_20261105`

`D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_PF_input_repair_mirrored_20261105`

Both have compact iteration records, a best trace, and fresh-process verification.
They are marked experimental, not accepted, and not installed in the body adapter.

## Downstream training hold

`pf_design_gate.py` now blocks new s3k/syn6 C/D, MN-readout, and body-optimization
launches through the current fitting entry points. Historical replay and
plotting remain available. Isolated PF-repair experiments remain available.
The hold must be cleared only after a revised implemented PF design has declared
timing targets, measured timing acceptance, and fresh-process evidence. The hold
is explicit; it does not pretend a repaired PF controller already exists.

Next design work must establish RG transition phase against its intended inputs,
fit or construct each PF channel's actual window before MN fitting, and test
the intended overlap and mixed-channel shutoff explicitly. S3k needs declared
PF windows first. Syn6 needs a principled mixed-channel and family-inhibition
revision, or an upstream functional subnetwork that can realize the six actual
primitive profiles. Repeat the full neural/physical curriculum only after that.

## Review files

Syn6 timing overlay and table:

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\circuit_trace_20261105\syn6\timing_audit\PF_TIMING_AUDIT.png`

`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\circuit_trace_20261105\syn6\timing_audit\PF_TIMING_AUDIT.json`

The corresponding s3k timing files are in the sibling `s3k\timing_audit` folder.
The audit and experimental fitter are `pf_timing_audit.py` and `fit_pf_inputs.py`
under `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104`.
