# Gait neural pathway composition review

Role: **Gait neural pathway composition reviewer**. Machine: easteregg2.
This review covers graded s3k, w2lvar, Syn6 K6, and the current Hunt/PT four-cell
architecture. It performs neural electrical calculations only. No body binding,
physics integration, tendon/compliance or damping application, optimization,
Hunt fitting, curriculum, or manuscript change was performed. Training remains
disabled until Ben explicitly approves it.

## Reviewed conditions

| Family | Saved standing settings retained | Separate expression fixture | Central current contract |
|---|---|---|---|
| s3k,410 neurons | PF-to-MN and IBEXC-to-MN gains remain zero; coincidence is inventory-only |162 active PF-to-MN slots and ten coincidence gates | No existing qualifying direct-central groups; inventory-only |
| w2lvar,888 neurons | Same saved zero PF/load-motor expression retained |106 active PF-to-MN slots and ten coincidence gates |24 existing groups,520 edges |
| Syn6 K6,794 neurons |318 PF-to-MN slots; saved gain multipliers retained, but right-canonical projection is an explicitly new condition |318 active PF-to-MN slots and ten coincidence gates | No existing qualifying direct-central groups; inventory-only |
| Current Hunt/PT,1033 neurons | Current four-cell RF/HAM architecture, shared existing paired gain parameters | Actual fitted and separate extracted-template RF/HAM operators reviewed | Existing local-reflex/timing/coincidence composition verified; no inherited canonical-family or body pass |

The new expression fixtures use canonical source templates with nonzero PF/load
paths. They do not alter saved standing configurations or inherit their physical
results. The source-default posture/II-reference flags are restored for these
separately named fixtures. Their original128-step seeded builder anchors are
verified before design transforms.

`central_afferent_budget.build_saved` now accepts an optional projection callback.
It first reconstructs and validates the **unchanged legacy saved graph**. With a
callback it rebuilds through the same extracted constructor, replacing only its
local projection function. Syn6's right-canonical policy therefore runs before
saved gain transforms, instead of averaging the source's unrelated left/right
NMF columns. The new candidate explicitly does not claim the legacy averaged
graph hash. No production builder or source K6 factor file is replaced.

## Complete electrical checks

Each of the six Gait conditions composes the final graded FSA profile, qualifying
central budgets, eligible coincidence circuit, and accelerated wrapper. Final
graded reapplication refreshes its whole-matrix/cache fingerprint. Graded terminal
and central conductance selections are disjoint; every active scoped verifier
passes after all transforms and actual neural steps.

All conditions check368 signed terminal pools, with maximum static error
4.44e-16 mV. Extracted SNS-Toolbox transients, two encoder-frequency probes,
actual compiled PF/current-budget probes, and eligible coincidence truth/release
tests pass. The per-encoder engineering phase budget is10 degrees at7.156 Hz;
it is not a whole-loop phase/stability bound or a rescaling of NaP/PF clocks.

The actual base-network step is exercised120 times through declared muscle-input
envelopes, with no encoder-envelope violation. An out-of-range adapter input is
rejected rather than clipped. Final dense/sparse128-step traces are exactly equal.
The source bounds are conditional runner bounds: Ia0..2.88, II0..2.92, and
Ib0..1.125 nA with gamma coactivation0, vestibular proprioceptive gain0, and
afferent multiplier1. Other input ports are zero in that test; this does not
claim body-measured ranges or simultaneous physical attainability of all maxima.

The load-gate0.225 nA reference uses the canonical full-stance/Fmax encoding.
It remains an explicit engineering reference until a live body's force/gain
normalization is bound. A masked standing output is not reactivated to make a
test pass. Existing inhibitory/composite and unqualified paths remain explicit
inventory entries; absence of a qualifying group is not relabeled as a repair.

## Current four-cell RF/HAM review

The primary source is **current1033-neuron Hunt/PT**, already containing four
cells per leg: RF and hamstring integrators plus their two inhibitory relays.
All eight actual cells and16 incoming/reciprocal edges are reviewed. Sartorius,
gracilis, and TFL remain in their separate groups.

Eight motor coefficients are active when the descending command enables them:
bilateral semimembranosus0.162680, semitendinosus0.155061, biceps femoris long
head0.203865, and rectus femoris0.0960335 uS. Command0 masks all eight; command0.5
halves them; command1 restores them. The authoritative PF/FF/bi buffers exactly
match the existing three-module composition evidence.

With actual fitted four-cell parameters and clamped named joint-PF inputs,
RF-pattern input gives RF2.698/HAM-0.720 mV; hamstring-pattern input gives
HAM4.390/RF-0.916 mV, identically on both sides. The separate canonical-template
operator also selects the requested antagonist. These are neural operator and
coefficient tests, not evidence of the contribution to body torque or stepping.

The older1027-neuron condition remains an untied historical comparison. Its
source PF/FF/reflex/bi paired mismatch counts are96/93/132/6. Timing mirroring
alone does not repair them. Current1033 has zero mismatches on existing paired
authoritative parameters. Its unmatched shared-lumbar values are only numerical
residue (maximum3.31e-14 uS); they are not treated as meaningful asymmetric drive.
No new midline drive is added automatically. The unbound PT mechanical pathway
is excluded from these electrical composition and four-cell tests.

## Control and mechanics scope

Each condition supplies a coverage ledger for six floating-base and17 internal
anatomical coordinates;105 native software coordinates are not105 anatomical
DoFs. Primary sagittal PF references, distributed muscle paths/local afferents,
foot/frontal restraint surrogates, lumbar/hip-rotation rig assistance, and
unsupported independent axis targets are distinguished. A path or a static
torque-map rank is not closed-loop control proof. Free3D walking is not a
prerequisite for this scoped design review.

Geometry, masses, inertias, source force laws, muscle/tendon dimensions, joint
limits, armature, damping, timestep, contact geometry, and apparatus settings
are unchanged by this work. Optional compliant-tendon and damping-repair code
is not imported or enabled. Existing mechanics reports are reference evidence
only. Zero vertical support would still not imply apparatus-free balance.

The rank8 review is in `SYN6_K6_K8_DESIGN_REVIEW.md`: dataset/factor matching to
production K6 is array-exact, but the improved inverse reconstruction is not a
forward production-RG/PF or walking comparison. K6 remains unchanged.

## Evidence and remaining work

Full external records:

- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\gait_full_composition_v1\{s3k,w2lvar,syn6}_COMPOSITION.json`.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\hunt_ibexc_composition\canonical_right_role_budget\CHECKPOINT.json`.
- The latter folder's `HUNT_FOURCELL_REVIEW.json` supplies registry-ready rows for
  Hunt/PT and four-cell RF/HAM, both pointing to the current1033 source.

Complete: scoped electrical composition, input encoding/domain checks, shared
parameter/cache consumption, independent pathway operators, backend equivalence,
source preservation, and declared control coverage. Remaining: all-family review
assembly by the parent session, any new operating-point/body calibration, neural
training, and curriculum/stepping validation. Ben's explicit approval is required
before training; these records do not automatically release any calculation gate.
