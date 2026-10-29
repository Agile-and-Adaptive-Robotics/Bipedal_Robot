# What the balance, time filters, and biarticular populations actually do

## Task-level feedback versus myostatic feedback

Ia muscle-length-rate and II muscle-length pathways are the local stretch/
myostatic branch. The implemented PT branch instead reads true model CoM and
torso orientation, computes a delayed body-level correction in Python, and
feeds four directional SNS populations plus four inhibitory relays. The SNS
does not infer CoM from spindle afferents. The inherited BAL controller also
contains engineering CoM/torso feedback, so additive PT feedback overlaps it.
No incremental harness-force benefit has been established for that addition.
This is not a physiological reconstruction of Peterka/Ting balance circuits.

## The misleading HUNT_BANK identifier

These cells were added as an engineering temporal basis; Hunt did not specify
them. Each leg has E/F low-pass filters at 20, 50, 100, 200, 400, and 700 ms,
plus corresponding inhibitory relay populations: 24 cells per leg. They let
conductance regression combine smoothed/delayed RG activity into PF waveforms.
They are neither new RG oscillators nor a required part of Hunt's four steps.
Diagrams now call them RG phase filters. The old IDs remain only to preserve
saved weight-file compatibility. Direct RG-to-PF fitting without this extra
basis remains a separate architectural comparison, not something already tested.

## Four-cell RF/HAM PF candidate

The current Deng-inspired implementation has separate single HAM and RECT
premotor integrators per leg. Adding one inhibitory relay from each to the
other makes a four-cell module: RF -> IN_RF -> inhibition of HAM, and
HAM -> IN_HAM -> inhibition of RF. The joint-PF inputs and MN output mappings
remain separate. Inhibition is tunable, not hard mutual exclusion; useful
coactivation must remain possible. The cells are graded/leaky in this variant,
not newly introduced intrinsic NaP oscillators.

This is a default-off candidate controlled by biarticular_pair_gain. It keeps
existing learned neuron IDs/weights compatible and preserves the passing
unpaired candidate. It has not been refitted through all Hunt phases yet.

## Actual Gait2392 function checks

Computed from the native equality-aware muscle Jacobian along the sampled
experimental gait trajectory, with the verified source/native sign adapter.
These are human Gait2392 results, not the CAD/BPA robot and not a full-RoM audit.
Median moment arms in mm, positive for the named column:

| Muscle | Hip flexion | Knee flexion | Hip adduction |
|---|---:|---:|---:|
| Rectus femoris | 44.45 | -45.29 | -11.99 |
| Biceps femoris long head | -60.29 | 31.54 | 16.08 |
| Gracilis | -6.43 | 25.65 | 56.79 |
| Sartorius | 53.49 | 3.35 | -28.14 |
| TFL | 38.40 | -12.82 | -48.47 |

Thus RF and hamstrings oppose each other across hip and knee in these poses.
Sartorius shares hip-flexion action with RF but has knee-flexion and abduction
actions. Gracilis is strongly adducting and knee-flexing, with a small median
hip-extension contribution here. Neither should be assigned wholesale to the
RF output merely because of one shared action. They can receive shared temporal
inputs while retaining function-specific PF inputs and independent MN outputs.
TFL is closer to RF in the sagittal action combination, but its strong abduction
action still needs separate control. Signs and leverage can depend on pose.

Raw sampled ranges: runs/biarticular_function_signs.json.
