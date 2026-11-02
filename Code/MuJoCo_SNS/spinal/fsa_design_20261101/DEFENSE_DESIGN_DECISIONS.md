# Controller design: defense answers and corrections

Session 1, easteregg2. This is a review/handoff document, not manuscript text. Proposed questions reflect the expertise Ben described; they are not statements from committee members. Existing result/configuration records remain authoritative.

## Design choices that can be defended

| Likely question | Defensible answer | Evidence to show |
| --- | --- | --- |
| Turcic: Why a mixed spiking/nonspiking controller? | It is an incremental conversion experiment. Keeping downstream graded circuits fixed isolates the effect of a spiking encoder or RG/PF replacement. This supports controlled comparison; it does not establish a fully spiking controller. Mixed nodes are compatible with the functional approach when signal transformations are explicit. | Exact node census, input/output units, circuit schematic, matched graded controls, real events, and source-output-off tests. Micklem2026 instead keeps central circuitry graded and spikes at actuator outputs. |
| Turcic: What is controlling the body? | Show the full controller, actuator, body, and external-apparatus loop. The Gait standing example uses neural recruitment without an online torque allocator after initialization. The CAD demonstration uses a literal neural pressure pair inside a52-route hybrid context. | Full apparatus force/torque measurements; zero pair clamp-current proof; actual controller/configuration paths. |
| Reddy: Why these sensory pathways? | Pathways are hypotheses from specified literature/model lineages. Ia/II/Ib signal identity, target cell, sign, and phase dependence must be checked individually; wiring similarity alone is not biological validation. | Per-edge source/rule ledger, actual physical source binding, reversals and gain fields, single-path pulse and outgoing-off controls. |
| Reddy: Are six NMF components six biological spinal synergies? | No. They are a reduced nonnegative representation of a particular inverse-activation dataset. They provide a useful coordination basis, not identification of a unique biological circuit. | Dataset provenance, rank/VAF curve, normalization, reconstruction and held-out results. Keep DiRusso primitives, Stark clusters, and NMF factors distinct. |
| Jiang: What is analytical and what is optimized? | Single-path equilibrium/transmission parameters can be derived; interacting nonlinear channels still require forward tests. The NMF/inverse-force problems, selected conductances, intrinsic timing fits, and mechanical objectives are separate mathematical problems. | Actual equations, decision variables, constraints, units, solver, stopping criterion, and consumed parameter files. |
| Jiang: Have you proved stability? | No global walking-stability proof is claimed. Fresh repeats, perturbation tests, and passing finite trials are scoped empirical evidence. Local transfer functions apply at stated operating points; equilibrium ROA and periodic-orbit stability require different analyses. | Appropriate local Jacobian/transfer-function calculation if completed, and exact finite validation protocols; do not label +/-1%parameter tests a full state-space ROA. |
| Hunt: Was the four-part neural training completed? | Report each operation with the same saved configuration: free rhythm, expected-afferent timing, feedforward motor output, then sensory refinement. The physical curriculum is a separate set of tests. A different controller's standing pass cannot complete a stepping controller's failed stage. | One configuration lineage and a per-operation/per-environment status table. |

## Errors or omissions to acknowledge, with corrections

| Issue | Acknowledgment | Correction and status |
| --- | --- | --- |
| Full2020FSA claim | No whole current spiking controller was documented as completing all seven Table2operations. | Full-family audit completed. New design profiles must derive applicable quantities and verify pathway gain/transients before training. No retrospective relabeling. |
| Body-size transfer | Rat/cat model timing and gains were not consistently justified from actual human/robot inertia, stiffness, delays, and signal ranges. | Sutton-based size/time review and human cached running/harmonic audit completed. New profiles must declare behavior, bandwidth, physical ranges, and actuator response separately. No mass-only clock multiplier. |
| Units across tools | Like-named capacitance fields have different numerical time/unit conventions. | Explicit seconds/ms and microfarad/nanofarad adapters;30ms at1uS means30nF in the physical-unit Simulink block. Unit equivalence checks precede transfer. |
| Additive branch called a gate | Some IBEXC branches add RG and load excitation; they are not multiplicative stance/load gates. Other syn6 branches differ. | Name actual compiled topology per family; quantify input contributions and test intended gating behavior before making a mechanism claim. |
| Missing central afferent routes | A functional schematic may omit or a source model may lack specific sensory-to-RG/PF connections. A paper's wiring cannot be assumed present. | Edge-by-edge Ia/II/Ib central-target audit, with source IDs, fields, and missing/conflicting-path reports. Preserve actual and proposed graphs separately. |
| CAD stance | Historical initialization has loop gaps, poorly grounded feet, and toe-out; the pressure-interface example has distal overlap. | Mechanical repair is prerequisite to improved balance claims. Separate laptop Simscape handoff; imported constraints, inertia, ground contact, and near-closed loaded equilibrium must be verified. |
| CAD coactivation | The ideal-force controller's5%-of-MIF floor is not a calibrated physiological or pneumatic activation floor. | Document the actual force law and distinguish force floors, neural activation, pressure, and passive tension. Keep zero-floor diagnostic separate from neural results. |
| Template/current-model confusion | Old editor snapshots were labelled current/winners. | Remove legacy entries from the default menu, preserve them as explicitly historical data, and identify exact current runtime graph/configuration. User saved tabs are not deleted. |
| Spiking burst scoring | A local-peak detector split ripple into extra bursts. | Preserve old results and document the versioned common burst-envelope rule with missing/extra/irregular/isolated-event controls. Do not change raw waveform gates or actual neural traces. |

## Order of work mandated by Ben

1. Answer the design/CAD questions; deliver schematic and laptop prompts; curate the editor list.
2. Apply appropriate design rules, physical/signal scaling, and individual-pathway verification to ALL active controller families, beginning with Gait2392.
3. Only after the global design/pathway gate passes, run Hunt training and the corrected physical curriculum with verified speed-ups.

Training/curriculum jobs must not start merely because the first Gait family passes. Graded rules, GLIF rules, retained ion-channel dynamics, and external engineering components need explicitly appropriate checks rather than a misleading common completion label. Existing published/submitted manuscript, legacy winners, and ZCode's eight-synergy study remain separate.
