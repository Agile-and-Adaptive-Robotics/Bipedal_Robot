# s3k and syn6: PF readout and local reflex follow-up

Machine: **easteregg2**. Session 1. This continues the 20261104 calculation
batch. Raw outputs are under `D:\Users\Ben Bolen\ControllerTraining_20261104`.
The inverted-pendulum work was left unchanged while these controllers were tested.

## What changed

The previous feedforward fit used one multiplier per muscle. The new fit allows
different multipliers on each existing PF-to-MN projection, with homologous
left/right projections sharing the same multiplier. No PF cells or connections
were added. Discrete MN current balance proposes bounded conductances, then the
actual complete neural circuit is integrated to evaluate every proposal. The
original waveform, peak, rhythmicity, and bilateral-phase gates remain in place.

Both families now pass offline feedforward (C) and afferented waveform (D)
checks. These letters name calibration blocks, not physical curriculum stages.

| Family | C waveform RMSE | D waveform RMSE | D peak RMSE | Fresh-process D repeat |
| --- | ---: | ---: | ---: | --- |
| s3k | 0.168362 | 0.173783 | 0.246645 | Exact objective match |
| syn6 | 0.161525 | 0.158450 | 0.246249 | Exact objective match |

Errors are in normalized MN-output/target-activation units. Targets are replayed
expected muscle signals. Later validation cycles repeat the same reference;
they are not independent subjects or proof of physical tracking.

### s3k relay finding

Reducing raw Ia/Ib/II input gains alone did not repair the waveform. In the
existing builder, the stance-group `IBEXC` relay is driven by both Ib input and
`RG_E`. Its MN output can therefore remain strong when raw Ib input is weak.
An eight-case comparison isolated the downstream Ib-relay gain: with II-relay
output held at 1, reducing Ib-relay output from 1 to 0.1 changed waveform RMSE
from approximately 0.34 to 0.174 and passed the unchanged peak gate. Reducing
II-relay output alone did not achieve that. The selected parameters are shared
bilaterally: local Ia 0.1, Ib 0.1, II 0.01; Ib-relay output 0.1; II-relay output 1.
All tested gains were positive. This is gain calibration, not synapse deletion.

This distinction is now explicit in the replay and body adapters. Previously,
the optimizer exposed raw sensory-input gains while the separately driven relay
output remained fixed. Archived builders and old checkpoint files were retained.

### syn6 ankle finding

The earlier physical-air candidate had fast ankle oscillation. Attenuating all
ankle-crossing sensory sources also weakened the knee, because gastrocnemius
crosses both joints. The selected diagnostic keeps knee-crossing sources intact
and scales other ankle-group local Ia and II outputs by 0.15. It uses shared
PF-to-MN expression 0.3 and additional global local-II scale 0.2114065546411311.
Its underlying D checkpoint is the original accepted syn6 D, **not** the new
waveform-readout fit above. These are two separate candidates.

This ankle candidate passed both physical air stages with identical saved
weights/settings: stage 1 period about 2.466 s; stage 2 period 2.005 s, knee
ranges 20.450/22.840 degrees. Stage-2 maximum extension stayed below zero
(-1.910/-1.937 degrees). The stage-2 repeat matches the corresponding fine-ladder
case; see numeric repeat audit in the portable result folder.

Over 5–14 s, ankle-angle spectral power above 3 Hz fell from 50.9/55.0% to
9.0/12.3%. Absolute high-frequency RMS fell from 4.52/4.84 degrees to
2.70/3.13 degrees. The residual ripple remains visible. The fractional-power
improvement also reflects greater slow ankle motion; it is not complete removal
of the oscillation. The <=15% spectral screen is an engineering diagnostic,
not a literature-derived physiological threshold.

These are **suspended air-stepping** results: measured foot normal force is zero.
Vertical suspension force is apparatus bookkeeping, not a minimum ground-walking
support requirement. Target cadence and the remaining ankle ripple need further
work before this candidate advances to posture or supported ground stepping.

The newly fitted syn6 C/D candidate passes offline, but its separate physical
expression ladder fails stage 2. It is not substituted for the air candidate.

## Exact checkpoints and records

| Role | Full path on easteregg2 |
| --- | --- |
| s3k accepted C | `D:\Users\Ben Bolen\ControllerTraining_20261104\waveform_edge_feasible_20261104\s3k` |
| s3k accepted D / relay comparison | `D:\Users\Ben Bolen\ControllerTraining_20261104\s3k_relay_output_20261104` |
| s3k physical ladder | `D:\Users\Ben Bolen\ControllerTraining_20261104\s3k_repaired_body_20261104` |
| syn6 accepted new C | `D:\Users\Ben Bolen\ControllerTraining_20261104\waveform_edge_fit_20261104\syn6` |
| syn6 accepted new D | `D:\Users\Ben Bolen\ControllerTraining_20261104\waveform_edge_D_20261104\syn6` |
| syn6 new-readout physical failures | `D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_waveform_body_20261104` |
| syn6 physical ankle baseline | `D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_ankle_preserve_biarticular_20261104\baseline` |
| syn6 fine ladder / selected both015 | `D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_ankle_fine_20261104` |
| syn6 same-configuration stages 1/2 | `D:\Users\Ben Bolen\ControllerTraining_20261104\syn6_ankle_selected_repeat_20261104` |
| Portable GIF, plots, and numeric audit | `D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Defense\Session1_Handoff_20261104\experiments_20261104\s3k_syn6_followup` |
| Main training record | `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\TRAINING_RECORD.md` |

Every completed physical run owns `body.json`, `body.result.json`,
`capture\aligned_trace.npz`, and the actual compiled `capture\model.mjb`.
Use that configuration and compiled model together; the graph ledger alone is
not a runnable checkpoint.

## Numerical and scope checks

Sparse compiled neural updates, existing compact physics evaluation,
preallocated traces, persistent scientific workers, Rustuna for bounded D fits,
and one BLAS thread were used. No more than three heavy jobs were launched at
once. Offline replay takes roughly 2–3 s for 14.8 s simulated neural time;
physical timing is recorded per run and includes concurrent-machine contention.
No new matched baseline speedup factor is claimed.

Sampled chemical-conductance Euler load ratios stayed below 1 in the tested
offline and live-body probes (maxima approximately 0.31–0.33). This argues against
that specific overshoot mechanism; ionic and recurrent stability are not proved
by this diagnostic. It is not a replacement for timestep convergence testing.

Predictive knee protection remains enabled. Foot contact is disabled in these
air tests. No heel wiring, tendon compliance, joint damping, muscle force law,
CAD routing, manuscript, or Overleaf content was changed in this follow-up.
The separate failed experiments remain available and were not overwritten.

## Physical s3k completion

The five-expression physical ladder completed without a same-configuration
stage-1/2 passer. At shared expression 0.2, stage 1 passes, but stage 2 hip ranges
are 127.6/142.4 degrees, exceeding the unchanged 80-degree upper bound.
Knee ranges are 28.4/26.0 degrees and knees remain flexed. The problem is not
absence of all rhythmic leg motion; the mechanical output is too large at the
hip and bilateral phase drifts. Higher expression candidates fail earlier gates.

Two additional paired-stage checks reduce PF-to-MN gains for primary hip muscle
groups to 0.3 or 0.5 of that expression. Both fail stage 1: 0.3 largely loses
knee cycling, and 0.5 gives knee ranges 15.5/17.9 degrees, below the 20-degree
minimum. Hip-primary groups include knee-crossing muscles, so this also changes
biarticular action. These observations support a joint/functional readout
recalibration rather than treating one global hip gain as a complete repair.
No failed candidate advanced to ground stages.

Those tests are at
`D:\Users\Ben Bolen\ControllerTraining_20261104\s3k_hip_expression_20261104`.
All bounded jobs in this follow-up have completed. The final source snapshot is
`D:\Users\Ben Bolen\ControllerTraining_20261104\source_snapshot_s3k_syn6_followup`.

Visualization correction (2026-11-05): the air-stage GIF previously displayed
an inactive ground plane, so feet appeared to pass through a floor. The renderer
now hides that plane and states ground contact is disabled. Saved qpos, forces,
controller parameters, and acceptance results are unchanged. The earlier source
snapshot remains historical; the repository renderer contains this later fix.
