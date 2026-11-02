# Non-Gait functional contracts and electrical checks

**Current active-family/evidence entry:** `NON_GAIT_AUTHOR_CHOICE_RECORD.md`. The later unilateral134 prototype is audit-only and excluded from training. Original public120 and its explicit NEW122 bilateral-biarticular candidate are distinct; reduced Rybak is separately included. All training requires Ben's explicit OK after review, regardless of local electrical results. Earlier component tables below remain lineage, not blanket model approval.

Machine: **easteregg2**. These are opt-in engineering profiles and isolated electrical checks. They do not change a source model, fit Hunt parameters, run a body, or pass the global curriculum gate. Gait2392 remains the first priority, and **all required families must finish the appropriate design/scaling/pathway review before training**. `CAMPAIGN_GATE.json` remains authoritative and unchanged by this module.

The reusable implementation is `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\non_gait_contracts.py`. The eight ledgers under `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\non_gait_ledgers` contain exact source paths/hashes, actual units, controlled subspaces, explicit component targets, proposed profiles, consumed-field proof links, and remaining limitations. `NON_GAIT_DESIGN_INDEX.json` maps them to all seven non-Gait families in the global gate; the CAD family has separate original-hybrid and pressure-pair ledgers.

## Implemented and checked components

### Current appropriate-design checkpoint

The current source-specific result is `NON_GAIT_APPROPRIATE_DESIGN_RECORD.md` in this folder. Its evidence supersedes the earlier isolated-component inventory below. The actual W2L108, W2L704, and curated rat134 runtime constructors now pass their complete scoped electrical-role checks, actual physical-input dispatch checks, nonlinear ionic-mechanism checks, and explicit first-entry feedback-current budgets. Independent fresh processes reproduce all three complete result JSON files exactly.

The curated rat has **29**, rather than the earlier 17, eligible independent input encoders: heterogeneous outgoing presynaptic windows are preserved and checked. W2L704's sixteen real sensory-spiking populations are tested through their actual event/filter destination, not through unused graded shadows. Native reference snapshots and fresh-proof records are outside the repository under `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical`.

Coordinate coverage means each actual coordinate is assigned a neural, passive, constraint, or apparatus role. The rat's mid-spine hinge is explicitly passive; W2L's toes are passive and its sliders are apparatus coordinates. Free 3D walking is outside those source models' declared scope, **not a prerequisite for Hunt training**. CAD mechanical defects remain gates for new body trials, not failures of its isolated electrical contracts. Reduced Deng14 is a neural mechanism with no instantiated body coordinates.

All eight ledgers expose `appropriate_local_design_checks_verified` separately from the parent-owned global design/scaling/pathway flags. Those global flags remain false pending review and release; this helper does not start training or modify the campaign gate. New feedback budgets are explicit engineering candidates, not a relabeling of historical controller behavior or body success.

### Earlier isolated-component checks, preserved for lineage

| Family/component | Explicit contract and implementation | Electrical result and limits |
| --- | --- | --- |
| W2L108 graded | 44 actual input-only physical-sensor relays preserve their own source current/voltage domains and `1/Size` DC gains, with a declared 10 ms physical pole. Twelve MNs preserve their source physical membrane timing; their admissible motor-voltage domains come from the actual SI adapter and 5–95% of the positive active sigmoid plateau. | Source-kernel endpoint/midpoint and motor-domain transient checks pass, including a fivefold biological-rate case with physical timing compensated. Errors ≤3.42e−13 mV. Size, ionic parameters, and all edges remain unchanged. This is not every pathway or a body test. |
| Curated rat134 | Seventeen eligible input-only relays and 16 actual MN outputs use the same explicit unit/clock adapters with their own source values. | Errors ≤1.57e−13 mV for motor transients; all selected sigmoid-domain points satisfy the actual adapter limits. Primary runtime is `variant_controller()`, not the base constructor alone. The earlier accepted intrinsic-A and rejected B studies are unchanged. |
| Original CAD262 hybrid | Preserve and expose the engineering command-to-MN function, rather than relabel it neural coordination. The original current-target function cancels background synaptic currents and implements `Vnext = V + dt/motor_tau * (5*u − V)`. | The actual source function follows this equation in a changing-background electrical fixture. Consumed `motor_tau=0.005 s`, giving a 3.915 ms discrete pole at 2 ms steps. Its comment saying 30 ms is not the consumed tracking setting; passive MN capacitance and the downstream 30 ms actuator filter are separate. |
| CAD262 free pressure pair | New, separate endpoint-FSA candidate for the six existing PF_E1/PF_E2/IBEXC→vas_med MN edges. Keep their 0.2:0.3:0.6 ratios; derive from the **actual compiled 8 mV reversal**. Coactive full inputs give 0.25 mV MN voltage, normalized command 0.05, and 226.5 N reference tension. Derive capacitance for a 30 ms **discrete** full-input pole. | Isolated actual source-branch checks pass at five input levels, with no MN target current or force floor. The exact function is rational, not perfectly linear: maximum departure from a linear endpoint interpolation is 0.794% of full scale. An extra 0.5% linearity aspiration remains missed and recorded. Pressure conversion, limits, lag, and geometry are unchanged and were not simulated. |
| W2L704 sensory carrier | `PhysicalCarrierPort` reuses the unchanged source 32-cell LIF encoder on an independent **physical** clock: 200 µs steps, 5 ms membrane, 1 ms refractory interval, 10 ms count filter, requested 0–200 Hz, normalized 0–1 input. It rejects out-of-domain input rather than silently clipping it. | Actual mean rates for 0/50/100/200 Hz requests are 0/49.156/98.031/192.313 Hz. They satisfy the predeclared 10% numerical carrier bound. Causal states and events in common physical bins match exactly across 1×/5× caller partitions; complete state continuation is exact. This does not replace the production hook or derive Table 2 conductance synapses. |
| Reduced Deng14 | `run_deng14_passive_profile` executes the **original fixed-tau ODE code** with a private capacitance vector. Preserve all six NaP cells at 5 nF, h parameter 350 ms, all edges, and all six IN leak poles at 5 ms. A separate candidate gives the two feedforward MNs a 10 ms leak pole, retaining their DC maps. | In the 200 ms electrical pulse test, no-profile and empty-profile results are exact. All 12 upstream voltage traces and the h trace remain byte-identical; only the two MN readouts change. Reconstructed source DC conductances differ by ≤1.78e−15 µS. No 20 s mechanism result or mechanical walker is claimed. |
| SNS-Simscape graded | Actual source `NonSpikingNeuron` masks consume 30 and 15 nF at 1 µS, producing 30 and 15 ms physical RC responses. | MATLAB R2025a pure-Simulink block tests pass with maximum voltage error ≤1.78e−14 mV. No Simscape body or Multibody license was needed. |
| SNS-Simscape spiking | Preserve the actual fixed-threshold/reset-below-rest LIF law: rest −52 mV, threshold −46 mV, reset −58 mV, 0.4 µS, 2 nF. Derive its exact current/rate function instead of treating it as reset-at-rest adaptive GLIF. | Actual source-block rates at 2.9/5.3 nA are 84.715107289602/204.811096240665 Hz, matching the analytical source-reset law to approximately 1e−12 relative error. The current installed library lacked `SpikingNeuron`; this negative is preserved. The source-builder copy was tested in isolation, then its exact built library was promoted to the durable external folder below and tested directly again. Original library bytes are unchanged. |

The native 10 ms sensor target is an engineering conditioning choice, not a universal physiological value: 50 samples per pole at the source 0.2 ms clock, separated from locomotor cadence. Relative to the retained rat-A 0.292 s context it is approximately T/29.2. The motor contracts preserve source timing and the actual muscle-voltage map rather than applying this sensor time constant to MNs. The 5–95% motor domain refers to the **positive active sigmoid plateau before length/velocity/series/passive effects**, not total force or measured physiological activation.

The new pair conductances, shared left/right, are approximately 0.00586510, 0.00879765, and 0.01759531 µS. These are **proposed electrical profiles**, not replacements for the saved accepted pressure-pair weights. IBEXC includes additive RG feedforward; it is not a multiplicative stance-load gate. The endpoint 5% command is not a minimum force floor or a pressure-state assignment. Other neural inputs are held at rest only in the explicitly isolated branch test.

## Clock equations and scope

Use physical time, source biological time, and wall-clock computation time as three different quantities. In the current W2L central wrapper, `t_bio = t_phys / clock_scale` and the physical neural substep is `native_dt * clock_scale`. Therefore `tau_phys = clock_scale * tau_source`, and a configured biological firing rate becomes `rate_phys = rate_source / clock_scale`. A value of 0.2 speeds that whole neural clock fivefold; it does not slow it. The Hill/MuJoCo clock is unchanged. The wrapper affects sensory and MN membranes, carrier/refractory/filter timing, central adaptation, readout, and startup windows, so it is not an RG-only intervention or a computational speed-up.

Equivalently, let `lambda = biological_time / physical_time`. Then `tau_phys = C_source/(G*lambda)`. To achieve a selected physical passive pole, the source clock-coded capacitance is `C_source = tau_target*G*lambda`, with the appropriate seconds/ms conversion and any declared operating synaptic conductance. The SNS adapter can instead derive the exact discrete pole for its Euler step. The native adapter sets its millisecond TimeConst explicitly. This compensates only the named selected cells; it does not silently rescale the whole circuit.

The canonical Newer W2L profile records an authored neural-clock scale of two **already serialized in its actual source fields**. Its runtime `neural_time_scale=1` must not apply that factor again. Rat134's curated runtime consumes its actual source TimeConst/gate fields and uses `dt_s*1000` in the source kernel; no extra universal body-size clock was found in the frozen wrapper. Preserve intrinsic source fits and choose sensory, MN, and actuator physical timing separately. Do not derive human timing from mass alone.

`PhysicalCarrierPort` is an opt-in implementation of this separation for the existing sensory encoder. A future caller must disable the corresponding old encoder, advance this port with actual physical elapsed time, bind the same named output proxy, and hold it between physical ticks. Do not run both encoders or advance it again on the accelerated biological clock. The port exports/restores membrane, refractory, filter, count/event, pending-clock, identity/configuration, and explicit empty-RNG state. Its rate mismatch remains measured rather than hidden by a configured maximum.

## Appropriate FSA application

- Graded passive ports use declared current/voltage domains, DC functions, local poles, and endpoint conductance equations. A local endpoint match is not proof of linearity over the entire range or of simultaneous-input equivalence.
- NaP/Ca channels retain their original dynamics. Their m/h variables, h timing, and adaptation currents are not GLIF adaptive-threshold m or tau_theta. The new pure-GLIF module must not overwrite those cells by parameter-name analogy.
- The source Simulink spike synapses **add** ginc and cap conductance. Table 2's reset-to-peak synapse has a different event law. Do not install an Eq. 28/29 peak as ginc and claim equivalence. The existing count-filter sensory port is also not a reset-to-peak conductance synapse.
- The original CAD hybrid current-target path is a separately labeled engineering controller. Its verified equation does not make its background cancellation a biomimetic recruitment pathway.

## Mechanical degrees of freedom are part of the contract

W2L and rat source limb controllers address three hip/knee/ankle coordinates per side. WalkingPath/rails, frozen roots, and external apparatus remove or stabilize global motions. These are not the same control problem as a free 23-DoF human model. The reduced Deng14 circuit is a neural mechanism with only two hip MN outputs, not a complete walker.

CAD has 17 independent model DOFs: floating pelvis six, hips six, two knee motions after rank-four loop closure, and lumbar three. Its global support and much of its posture regulation are supplied by apparatus and the hybrid allocator. Welded feet remove joint leverage from 20 ankle/toe channels. Neither the pressure pair nor the sagittal source circuitry independently establishes hip frontal/transverse, lumbar three-axis, subtalar/MTP, or free-root balance coverage. Every ledger records this subspace and the remaining coverage; the laptop repair handoff remains the mechanical prerequisite.

## Use the opt-in profiles

Each family ledger is its profile file; select its `opt_in_profile` field. Existing runners do **not** load these files implicitly. Apply a profile after constructing the exact designated controller and before any warmup, preparation, or state capture.

```python
from non_gait_contracts import apply_native_passive_profile
# controller must come from the actual source constructor, e.g. rat
# rat134_receptor_runtime.variant_controller(), not a substituted base build.
apply_native_passive_profile(controller, ledger_path,
                             neural_rate_multiplier=1.0)
```

For compiled CAD SNS:

```python
from non_gait_contracts import apply_sns_passive_profile
apply_sns_passive_profile(compiled_forward, ordered_neuron_names,
                         pressure_pair_ledger_path, clock_unit="s",
                         neural_rate_multiplier=1.0)
```

Only existing selected cells/edges are changed. Sparse conductance caches are updated with the core array; ionic cells, unknown fields, missing/ambiguous paths, and topology additions are rejected. `None` is a literal no-op. Profiles are marked as engineering candidates, not adopted body winners. A production integration must retain that default-off behavior and its own complete preparation-state audit.

For Deng14, `run_deng14_passive_profile(None, ...)` delegates to the original function. The proposed profile selectively changes passive capacitances through its original vectorized voltage equation, without editing the source. Its primary reference is `simulate("fixed", t_end=20000, dt=0.1, i_pulse=10, t_pulse=100, pulse_ms=1)` in **milliseconds**. The wrapper's public arguments are seconds. Simulink and the legacy NumPy script remain separate comparators with different PF→MN numbers; they were not merged into this profile.

## Reproduction and records

Interpreter: `D:\Anaconda\envs\myo\python.exe`. Repository: `D:\GitHub\Bipedal_Robot`.

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\non_gait_contracts.py' --generate --electrical-probes --candidate-probes
& 'D:\Anaconda\envs\myo\python.exe' 'D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\non_gait_contracts.py' --electrical-probes --candidate-probes --output-suffix _fresh
```

These commands execute only the named electrical checks. They do not open MuJoCo bodies or train controllers. Existing completed evidence need not be rerun. Standard Simulink block check: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\non_gait_simulink_probe.m`, run using MATLAB R2025a `-batch`. It accepts explicit library and record-directory arguments and preserves the supplied library's bytes.

**Durable, usable neural library and launcher:** `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\simulink_neural_blocks\README.md`. That folder contains the actual tested `SNS_Library.slx`, byte-identical source `sns_build_spiking.m`, copied `non_gait_simulink_probe.m`, `open_verified_neural_library.m`, a direct promoted-library verification record, and `PROMOTION_MANIFEST.json` with source/stage/copy hashes. The promoted library is byte-identical to the tested build, its direct repeat passes, and it no longer depends on disposable `tmp`. A same-name library check preserves any other `SNS_Library` already loaded in the user's MATLAB session.

Reusable contracts and the design index remain in this repository folder. Raw electrical evidence now lives at `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical`; the nine original files were moved with exact SHA256 preservation, including independent fresh repeats and negative diagnostics. `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\MOVE_MANIFEST.json` records every original/current path, byte count, and hash. Historical helper hashes in those unchanged records remain lineage: only output locations and evidence references subsequently changed. No neural experiments were repeated for this storage pass. Future raw writers refuse to overwrite an existing record; use a new `--output-suffix` or an explicit new Simulink record directory.

- `NON_GAIT_DESIGN_INDEX.json` and `non_gait_ledgers\*.json`: contracts, source/anchor hashes, units, subspaces, and application recipes.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_ADAPTER_ELECTRICAL_CHECK.json` and `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_ADAPTER_ELECTRICAL_CHECK_fresh.json`: default no-ops, actual field/cache consumption, ionic rejection, clock compensation, and analytical current/rate inversion.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_CANDIDATE_COMPONENT_CHECKS.json` and `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_CANDIDATE_COMPONENT_CHECKS_fresh.json`: all scoped native, CAD, carrier, and Deng14 component checks.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_FRESH_PROCESS_CHECK.json`: independent fresh-process record hashes; complete JSON equality is checked, not only rounded scores.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_SIMULINK_BLOCK_CHECK.json`: actual block readback, response/reset events, source-library preservation, and staged-library paths.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_CAD_REVERSAL_UNIT_FAILURE.json`: the first new design test used an incorrect 20 mV reversal assumption; it failed and is preserved. The correction reads the actual compiled 8 mV field. Old accepted CAD results were not changed.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_ADAPTER_PRECISION_DIAGNOSTIC.json`: an overstrict absolute rate-roundtrip tolerance failed near the source LIF threshold. The recorded 10 Hz error is approximately 2.1e−9 Hz; the final analytical numerical check declares a relative 1e−9 tolerance. This was not a physical/held-out acceptance-gate change.
- `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design\non_gait_electrical\NON_GAIT_SIMULINK_MISSING_SPIKING_DIAGNOSTIC.json`: the initial installed-library missing-block result, retained alongside the staged-source test.

The primary rat anchor is the curated 134-cell receptor variant and its `variant_controller()` wrapper. The refinement12 Hunt-B result remains a **secondary rejected diagnostic**, despite a nominal RG pass. CAD's primary body anchor remains the historical mixed support085 recording with crossing/gap limitations, not the failed repair fixtures. Deng14 remains a mechanism reference with no mechanical winner. Root's `AUDIT_ANCHORS.md/.json` specifies the campaign anchors.

**Remaining gate:** selected components now have executable contracts and electrical proof. Full-family pathway signs/net effects, recurrent interactions, physical input coverage, production-hook preparation parity, mechanical axis coverage, and any new whole-body validation still require review. No ledger, unit test, or old body pass automatically grants the global all-family design/training gate. No source assembly, model XML/SLX, original controller weights, manuscript, automation, or ZCode eight-synergy study was changed.
