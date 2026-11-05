# Fitted s3k and syn6: circuits beside neural activity

Machine: **easteregg2**. Session 1, 2026-11-05.

These figures reconstruct the saved **offline D fits**, not a generic template
or an idealized literature circuit. Each reconstructed circuit must reproduce
the saved objective and RG/PF/MN/target arrays exactly before rendering.
PF-to-MN conductances have matching left/right homologues: maximum difference
0.0 µS, with no missing mirrored projections. Figures display the right side.

## Expanded figures

- [s3k: schematic left, neural traces right](s3k/s3k_D_NEURAL_SCHEMATIC.png)
  ([PDF](s3k/s3k_D_NEURAL_SCHEMATIC.pdf), [editable SVG](s3k/s3k_D_NEURAL_SCHEMATIC.svg)).
- [syn6: schematic left, neural traces right](syn6/syn6_D_NEURAL_SCHEMATIC.png)
  ([PDF](syn6/syn6_D_NEURAL_SCHEMATIC.pdf), [editable SVG](syn6/syn6_D_NEURAL_SCHEMATIC.svg)).

The lower-left panels show summed PF-to-MN conductance by primary muscle group.
These sums are neither actuator force nor joint torque. The shared MN box is
a container; it does not assert all-to-all connections. Every individual muscle
projection is available in the maps and CSVs below. Inhibitory connections end
in filled circles; excitatory connections end in triangles.

## Which PF innervates which muscles?

| Controller | Direct RG-to-PF drive in these checkpoints |
| --- | --- |
| s3k | RG-E → E1/E2; RG-F → F1/F2 |
| syn6 | RG-E → S2/S3/S6; RG-F → S1/S4; both RG-E and RG-F → S5 |

The s3k PF cells are four **shared phase channels**, not a separate four-cell
PF oscillator for every joint. E/F names identify phase families. They do not
mean every recruited muscle extends/flexes every joint. For example, rectus
femoris receives all four s3k PF channels; its hip-flexion and knee-extension
roles cannot be summarized by one global flexor label. Iliacus receives F1/F2;
the vastus intermedius receives E1/E2; soleus receives E2. Tibialis anterior has
E1 as well as F1/F2 input in this fitted circuit.

Syn6 instead uses six overlapping primitive channels, with data-derived initial
projections and subsequent fitted per-edge gains. It is not six independent
joint-antagonist pairs. Iliacus receives S2/S3/S4/S5; rectus femoris receives
S2/S3; soleus receives S3/S5/S6. These facts describe the reconstructed wiring,
not a claim that every projection is an optimal or validated design.

- Full muscle maps: [s3k](s3k/s3k_PF_MUSCLE_MAP.png),
  [syn6](syn6/syn6_PF_MUSCLE_MAP.png).
- Exact channel-by-channel target lists and weights:
  [s3k](s3k/PF_ASSIGNMENTS.md), [syn6](syn6/PF_ASSIGNMENTS.md).
- Full-precision matrices: `s3k/PF_TO_MUSCLE.csv`,
  `syn6/PF_TO_MUSCLE.csv`; grouped sums: `PF_TO_GROUP.csv` in each folder.

All positive projections are retained, including small weights and zero-drive
muscle rows. Group labels come from the existing `muscle_map.py`; they are not
a complete anatomical moment-arm map. In particular, its historical sartorius
and gracilis group simplifications should not be interpreted as anatomy.

## Do flexors have Ib feedback?

**Yes locally; the checked pure flexors have no central Ib projection.**

| Checked source muscles | s3k local route | syn6 local route | Ib → RG/PF |
| --- | --- | --- | --- |
| Iliacus, short-head biceps femoris, tibialis anterior | Ib inhibits the homonymous MN directly | Ib excites IBIN; IBIN inhibits the homonymous MN | Absent for these muscles |

The fitted s3k direct inhibitory conductance is 0.035 µS for those muscles.
Syn6 has Ib→IBIN 0.0609156 µS and IBIN→MN 0.5 µS. These are pathway
conductances, not additive current or force contributions. This exposes a real
topological difference between the two checkpoints; an Ib relay was not drawn
into s3k merely to resemble the syn6 circuit.

Central Ib in the composed profile is restricted to selected stance sources:
gluteus maximus subdivisions, the vasti, and soleus. It reaches RG-E and the
E-family PF channels that carry those muscles, with a load-dependent gate.
For syn6, S5 receives both RG drives but is not assigned an E-only Ib route.
An absent flexor Ib central route is distinct from a temporarily zero gated
extensor Ib conductance. No connections were added or repaired in this drawing task.

- Actual iliacus flexor subnetworks: [s3k](s3k/s3k_FLEXOR_IB_SUBNETWORK.png),
  [syn6](syn6/syn6_FLEXOR_IB_SUBNETWORK.png).
- Direct sensory-muscle origin → PF maps:
  [s3k](s3k/s3k_PF_FEEDBACK_MAP.png),
  [syn6](syn6/syn6_PF_FEEDBACK_MAP.png).
- `AFFERENT_TO_PF.csv` in each folder records signs, fitted snapshot conductance,
  reference gate budget, shared gain, reversal potential, threshold, and saturation.

The flexor panels show homonymous reflex and carrying-PF routes. Antagonist
interneurons, recurrent feedback, descending commands, and commissurals are
omitted there rather than implicitly declared absent. The full external ledgers
include both sides and all reconstructed compiled edges.

## Scope and exact paths

The traces replay expected afferents; they do not integrate a mechanical body.
The syn6 checkpoint here is the **new offline D readout**, not the older-readout
syn6 candidate used for the improved air GIF. Their weights must not be mixed.

Fitted checkpoints:

`D:\Users\Ben Bolen\ControllerTraining_20261104\s3k_relay_output_20261104`

`D:\Users\Ben Bolen\ControllerTraining_20261104\waveform_edge_D_20261104\syn6`

Full bilateral ledgers:

`D:\Users\Ben Bolen\ControllerTraining_20261104\circuit_trace_audit_20261105\s3k\FITTED_EDGE_LEDGER.json`

`D:\Users\Ben Bolen\ControllerTraining_20261104\circuit_trace_audit_20261105\syn6\FITTED_EDGE_LEDGER.json`

Renderer/extractor:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\circuit_trace_figures.py`

Every figure has PNG, PDF, and editable SVG versions. `CIRCUIT_SUMMARY.json`
contains checkpoint identity, exact trace-reconstruction checks, mirror checks,
direct RG/PF routes, and per-muscle local/central Ib paths. No neural wiring,
muscle physics, curriculum result, or manuscript was changed to make a diagram.
