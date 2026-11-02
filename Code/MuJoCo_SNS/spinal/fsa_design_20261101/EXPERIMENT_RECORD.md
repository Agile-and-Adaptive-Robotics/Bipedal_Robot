# FSA and scaling design record

Session 1. Machine: **easteregg2**. Repository: `D:\GitHub\Bipedal_Robot`.

This record is separate from CHATGPT_REPORT.md and the historical training records. No manuscript or Overleaf edits belong to this campaign. ZCode owns the eight-synergy investigation; the existing K6 basis is preserved.

## Required sequence

Ben directed Gait2392 design first, followed by the other active families. Appropriate functional design, signal and physical-time scaling, individual pathways, and control coverage must be checked across **all** families before the next Hunt training or curriculum calculations. Design calculations and isolated electrical tests are permitted. No family receives a training release merely because another family passes.

`CAMPAIGN_GATE.json` remains closed. `design_gate.py` requires reviewed per-family evidence and current hashes for source, applied profile, design, and pathway results. Six prerequisite tests pass: reviewed fixtures release, boolean-only release fails, changed source invalidates evidence, missing control coverage fails, a second pending family prevents release, and unknown operations are rejected. The test fixtures are deleted from their owned temporary directories.

The existing historical runners are preserved. New campaign launchers must call `require_calculation_release()` before launching numerical training. This helper is not a claim that every old standalone script has been retrofitted with a global lock.

## Completed handoffs

- `DEFENSE_DESIGN_DECISIONS.md`: committee questions, defensible design choices, and correction record.
- `LAPTOP_CAD_SIMSCAPE_REPAIR_HANDOFF.md`: exact machine/file map, mechanical defects, and copyable laptop prompt.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\FUNCTIONAL_SUBNETWORK_OTHER_MACHINE_PROMPT.md`: actual-model graph/pathway audit and editable functional-subnetwork drawing prompt.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\CONNECTOME_TEMPLATE_CATALOG.md`: old templates hidden behind Legacy, preserved data, and current runtime graphs marked import-required. No user localStorage or live drawing was changed.

## Design evidence

### Dimension and torque coverage

`MODEL_DIMENSION_AND_CONTROL.md` explains six pelvis position/orientation coordinates plus 17 internal coordinates (seven per leg and three lumbar). It separates low-dimensional rhythmic coordination from the additional feedback, passive dynamics, ground contact, and apparatus stabilization.

The static `syn6_torque_span.py` diagnostic used four existing walking postures, actual native muscle active-force derivatives, and equality-aware anatomical moment arms. At each posture:

| Activation-to-internal-torque map | Numerical rank |
| --- | ---: |
| Full 92 muscles | 17 |
| 86 muscles represented in the current NMF basis | 17 |
| Fixed bilateral K6 activation basis, 12 primitive inputs | 12 |

Moment-arm finite-difference discrepancy at the checked posture was `1.8189912900576366e-11 m`. Source/model/basis/reference hashes were unchanged. The record is `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\SYN6_TORQUE_SPAN.json`.

These are static activation-map results, not a measured rank of the entire nonlinear SNS or a walking-stability result. Reflex/posture channels are additional inputs. Input bounds, unilateral muscle forces, torque margins, contact modes, and temporal dynamics require their own checks.

### Graded Gait2392 design

Optional profiles and default-off entry points are being checked for s3k, w2lvar, and syn6. The profile targets explicitly declared electrical functions; it preserves NaP rhythm dynamics, PF timing, MN/actuator dynamics, and existing default numerical behavior. A unit-current calibration amplitude is not called a universal physical sensor maximum. Full incoming-pool shunting and final compiled/mirrored/accelerated consumption must be checked separately from isolated endpoint algebra.

Source capacitance numbers use a seconds clock. Numerical C=0.03 with G=1 has a 30 ms leak time constant and corresponds to 30 nF if the physical conductance is 1 microSiemens. It must not become 0.03 nF in the Simulink physical-unit convention.

### Spiking functional design

`table2_glif.py` implements the appropriate seven-step pure-GLIF derivation, explicit physical units, and actual reset-and-decay synaptic conductances. NaP/Ca bursting states are explicitly outside that pure-GLIF derivation. They are retained or require their own ionic-spiking design.

The literal published conductance mapping and the selected spiking receiving-cell model produced an approximately 11–12% transmission excess in isolated tests. This negative result is retained; it is not repaired by an undocumented gain multiplier. `conductance_matched_design.py` derives a separate full-conductance firing-rate equation and solves it numerically by monotone scalar bisection. This is a design root solve, not Hunt training or an empirical controller sweep. It explicitly distinguishes rate matching in a spiking receiver from voltage matching in a graded receiver.

The corrected isolated transmission tests reported approximately unit gain across 25%, 50%, and 100% input and matched scalar/accelerated states within floating-point roundoff. Full results and applicability limits are in the agent's compact external design records. These module results do not themselves qualify a complete controller family or release the global gate.

### Non-Gait controller contracts

`NON_GAIT_DESIGN_INDEX.json` and `non_gait_ledgers` cover actual W2L, rat, CAD, SNS-Simscape, and reduced Deng sources. They distinguish verified electrical adapters from uncompleted family design/scaling/pathway work.

Clock factors must be interpreted from their actual implementations. In the current W2L spiking runtime, time and interval are divided by `clock_scale`; a factor 0.2 therefore gives five times the neural rate in physical seconds and one-fifth the neural time constants. It is not a cadence-only factor. Other controllers multiply time constants directly. A canonical profile whose source parameters already incorporate scaling must not be scaled twice. Physiological cadence, sensory encoding, spike carriers, synaptic/readout filters, and actuator dynamics need separately stated physical-time contracts.

## Current release status

The all-family gate is **closed**. Module-level checks are progressing; the full per-family prerequisite set has not yet passed. There have been no new Hunt parameter-training jobs, curriculum trials, whole-body optimization, source-model overwrite, or NMF refitting in this design phase.
