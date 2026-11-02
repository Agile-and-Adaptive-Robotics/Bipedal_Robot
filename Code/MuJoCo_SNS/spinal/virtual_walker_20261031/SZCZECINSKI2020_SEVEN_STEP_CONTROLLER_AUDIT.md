# Szczecinski2020 Table2: audit across controller families

Session1, easteregg2. Read-only review of the supplied PDF and existing source. No controller, source model, reference, physics, manuscript, or Zotero change; no simulations.

## Exact source and seven steps

Nicholas S. Szczecinski, Roger D. Quinn, and Alexander J. Hunt(2020), *Extending the Functional Subnetwork Approach to a Generalized Linear Integrate-and-Fire Neuron Model*, Frontiers in Neurorobotics14:577804. [Primary paper](https://doi.org/10.3389/fnbot.2020.577804).

User-supplied PDF read with PyMuPDF:
`D:\Zotero\storage\WUIWIT2B\Szczecinski et al. - 2020 - Extending the Functional Subnetwork Approach to a .pdf`.

**Table2 is on printed/PDF page11**, immediately after §3.4 on page10. Its title is “Expanded functional subnetwork approach for designing spiking pathways.” The exact goal wording is:

| Step | Exact Table2 goal | Parameter / equation |
| --- | --- | --- |
| 1 | Set network-wide activation parameters. | Fmax, R, θ0; Eq24 |
| 2 | For each neuron, set transient response type based on function in the network. | m; increasing rate:m<0; constant:m=0; decreasing:m>0 |
| 3 | For each neuron, set the duration of transient firing based on its function in the network. | τθ; Eq26 |
| 4 | For each neuron, calculate bias current. | Ibias; Eq21 |
| 5 | For each neuron, ensure fsp∈[0,Fmax]. | τmem; Eq23 |
| 6 | For each synapse, limit its non-linearity. | τs; Eq27 |
| 7 | For each synapse, set synaptic gain ksyn based on function in the network. | Gmax; Eq29 |

The equations, visually checked against the supplied table, are:

```
fsp,approx = (Fmax/R) Ubar_infinity
tau_theta = tau_mem_bar (1 - m/2)
I_bias = G_mem theta_0 / (2 - m)
tau_mem = (R/Fmax) (1 - m/2) / theta_0
tau_s <= -1 / (Fmax ln(delta))
G_max = k_syn R / [(E_s - k_syn R) tau_s Fmax]
```

Use the paper's voltage/conductance normalization and consistent time/rate units. `m` here is the GLIF adaptive-threshold coefficient, not the persistent-Na activation gate also named m in other models. `tau_theta` is not a persistent-Na inactivation constant or an adaptation-current time constant. The table is an analytical design procedure, distinct from Hunt's four-part training method and the mechanical curriculum.

## Status convention

**C**: exact table operation is implemented and verified for the scoped component. **P**: a related function/calibration exists, but the specified2020 analytical operation is incomplete or uses a different method. **M**: needed for a literal2020 spiking implementation, but no implementation/evidence was found. **—**: the variable is absent from this graded model, so the literal spiking step does not apply; this does not by itself invalidate the graded model.

No audited whole controller currently has evidence sufficient to claim **all seven2020 steps completed**. This does not mean none uses FSA: syn6 explicitly uses2017 FSA signal transmission, and newer spiking ports include analytical rate/filter mappings. It means those narrower results should be named accurately.

| Actual family / implementation | 1 scales | 2 transient m | 3 τθ | 4 bias | 5 rate range | 6 δ/τs | 7 functional Gmax | Concrete interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| s3k graded | P | — | — | P | — | — | P | Global voltage scale, chosen membrane/NaP timing, literature/group weights, zero/rest bias; no GLIF conversion. |
| w2lvar graded | P | — | — | P | — | — | P | Same graded base with W2L circuit-rule conductances and per-muscle feedback. Functional routing is explicit; exact function-derived whole-network conductances are not. |
| syn6 graded | P | — | — | P | — | — | P | R5mV and2017Eq18 are explicit for positive PF→MN weights; the other pathways use source/rule/tuned values. Multi-input dynamic fit exists as a separate analysis, not evidence that its output replaces every production connection. |
| biw2lvar-Hunt-PT graded | P | — | — | P | — | — | P | Intrinsic/sensory/motor/reflex fits and task-balance pathways exist; these are not2020GLIF-threshold design. |
| functional PF / four-cell RF-HAM graded | P | — | — | P | — | — | P | Explicit four-cell topology, fitted PF/MN gains, membrane-bank time constants, and tied-learning variants. A temporal-bank fit is not an Eq26τθ calculation. |
| W2L108 graded | P | — | — | P | — | — | P | Imported named types/fields and declared engineering changes, not a newly derived2020spiking circuit. |
| W2L704 actual spiking +80graded hybrid | P | P | P | P | P | M | P | Explicit scales, fixed threshold/reset/refractory, Ca gates plus adaptation current, spike-count decoding, inherited/scaled conductances. Fixed threshold is compatible with m=0, but the adaptation current is not an implementation of adaptiveθ. No δ-based Eq27 or Eq29 derivation for all paths. |
| rat134 graded source/engineering variants | P | — | — | P | — | — | P | Exact source transfer fields, own-receptor normalization, acceptedA intrinsic fit, and pairedB engineering weights. No threshold-rate conversion. The nominalB-only intermediate remains held-out rejected. |
| CAD262 original graded hybrid allocator | P | — | — | P | — | — | P | Rebuild of the stock graded network on54actuators; engineering motor-target allocator is not a Table2/FSA derivation or a calibrated BPA model. |
| Current Gait central spiking prototype | P | P | P | P | P | M | P | Copies sourceNaP/Cm/Gm/bias and uses fixed-threshold(m=0-compatible) population/reset/readout transformations. Source/hyperparameter fits are real, but not Eq21–29 GLIF design. |
| Older production spiking mirror (`AARL_NET=spiking`) | P | P | P | P | P | P | P | Has adaptive-threshold increments, rate/readout calibrations, and measured synaptic scaling. Its threshold increment/time constant policy differs from Table2m, bias is selected, and τs is not shown derived from δ. Do not conflate this implementation with the newer central population ports. |
| SNS-Simscape graded knee/beverage reflex demos | P | — | — | P | — | — | P | Explicit physical-unit gradedRC cells and muscle feedback, chosen operating points/gains, verified response; no spiking-rate derivation. |
| SNS-Simscape fixed-threshold spiking knee/beverage/RG-MN demos | P | P | — | P | P | P | P | Fixed threshold is the m=0 special case. Mean-conductance/rate matching exists, but no full2020network-wideFmax/range and δ/Eq29 calibration. Adaptive-threshold block is explicitly future work in the builder. |
| Reduced Deng14-cell ODE/Simulink circuit used in dissertation | P | — | — | P | — | — | P | Literature-derived gradedNaP RG/PF/IN/MN parameters and perturbations; no literal GLIF conversion. Fixed h-gate timing is not thresholdτθ. |

**Completion versus validation:** a formula being coded is not a complete end-to-end proof. In particular, spike-rate clipping or a configured maximum is not verification of an accurate0..Fmax encoding over the actual input domain. Likewise an average-current calibration is not automatic proof of appropriate transients under recurrent feedback.

## Source-code evidence and scoped missing work

All paths below are under `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\`.

- **s3k:** `spinal\params.py:26` defines E_HI5mV; lines30–43 define membrane/filter time constants; lines55–62 defineNaP parameters. `spinal\build_network.py:199` creates a gradedRC cell with Gm1,rest0,bias0; line206 defines graded saturation/reversal; lines419–429 create NaP RGs. These are source-selected graded parameters, notGLIFTable2 formulas.
- **w2lvar:** `spinal\build_network_w2lvar.py:1` documents copied gradedRG, hip and knee/anklePF pairs, and rule-based gains; it inherits the same cell construction. Literal Fmax/θ/m/τs variables are absent.
- **syn6:** `spinal\build_network_syn6.py:285` imports `analytical_conductance`; lines480–488 assign PF membrane time constants; lines642–653 wire positive NMF weights through Eq18. `spinal\fsa_backsolve.py:132–137` computes `k R Gm/(ΔE−kR)`; lines140 onward fit constant nonnegative conductances to the full membrane current equation. The dynamic fit and forward-requirements analysis are not automatically installed in syn6. `spinal\synergy_model.py:17–42` explicitly separates basis export, coefficient replay, and suspended open-loop demonstration.
- **Hunt/PT and RF-HAM:** `spinal\hunt_training\intrinsic.py:22–30` builds gradedNaP cells; lines95–100 fit period via τh. `controller.py:36–38` creates temporal-bank cells at20–700ms; lines155–178 load/consume fitted conductance arrays. `train_fourcell_hunt.py:19–52` checksA/B then fitsPF/BI/MN/reflex stages. `spinal\functional_gait2392\network.py:31–43` constructs explicit graded four-cell PF pairs. These establish actual function, fitting, and testing, not literal adaptive-threshold formulas.
- **W2L108/rat134:** `spinal\virtual_walker_20261031\native_kernel.py` implements source graded conductance/ion-gate dynamics. Rat134 actual editable/native fields and adapters are archived at `D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\Rat132_bounded_diagnostic\Curated134_ReceptorGeometry\`; acceptedA and full-source/engineering edge distinctions are documented in its records. No claiming a2020spiking derivation merely because the native module is named IntegrateFireSim.
- **W2L704:** `spinal\virtual_walker_20261031\spiking_central_runtime.py:12–22` uses fixed threshold/reset, adaptation-current decay/jump, and normalized spike filtering; lines31–37 select12central pools and copy source fields. Lines63–88 modify inherited PF→MN gains. `spiking_sensory.py:14–28` derives LIF drive from a requested firing rate, with reset/refractory and normalized count filtering. This is useful analytical encoding, but not Table2's GLIFm/τθ/Ibias/τmem/Gmax derivation.
- **Generic sensory-spiking port:** `spinal\virtual_walker_20261031\spiking_sns_inputs.py:20–23` genuinely derives a filter time constant from the source discrete decay: `tau=−dt/log(1−dt*Gm/Cm)`. It uses200Hz,5ms cell membrane, and1ms refractory policies. That filter match is distinct from Eq27's δ nonlinearity criterion.
- **Current Gait central:** `spinal\virtual_walker_20261031\gait_central_spiking.py:31–35` specifies threshold4mV(RG)/1mV(PF), reset2.5/0, and1ms refractory; line35 uses a25ms count filter. Lines43–50 copy sourceCm/Gm/Ibias/NaP into16central populations. Its later sharedNaP-hτ calibration is not Table2τθ. The inverse rate-to-voltage readout is a derived interface, not function-derived outgoingGmax.
- **Original CAD262:** `robot_torso\run_neural_robot.py:117` rebuilds the stock network from the actual54actuator list. Its motor-target allocator is separate from neural rhythm generation; it must not be counted as neurons supplying an analytically derived BPA pressure law.
- **Older spiking mirror:** `spinal\build_network_spiking.py:94–116` sets thresholds, bias, increments, and threshold timing; lines117–140 specify decay/rate/calibration factors. `calibrate_spiking.py` calibrates current/rate transmission. These empirical/mean-conductance calibrations provide partial functional evidence but not completion of every2020table equation.
- **SNS-Simscape graded/spiking:** actual source is `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\`. `sns_build_library.m:22–24,84–89` gives gradedRC equations/units; `sns_build_spiking.m:17–23` explicitly chooses fixed thresholds and leaves adaptive thresholds for future work. Its lines44–50 state interface units, and line109 implements the seconds-to-ms derivative conversion. `demos\sns_build_knee_spiking.m:18–33` describes first-order/hand-tuned mean-conductance matching,5ms sensory membrane,10ms synaptic decay, and chosen thresholds; lines67 onward install those values. `sns_units_test_spiking.m` and `Code\MuJoCo_SNS\spinal\virtual_walker_20261031\verify_existing_spiking_reflex.m` provide actual interface/behavior verification artifacts, distinct from Table2 derivation.
- **Reduced Deng14:** `spinal\deng_cpg_ode.py:23–35` defines5nF/1µS, absolute resting/reversal potentials, fixed350ms h-gate, and the14named cells; lines42–66 define inherited neural connections. `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos\sns_build_deng_cpg.m:45–73` defines the named conductances,5nF/1µS source units, and350ms h timing. The ODE and block equations are graded, despite “integrate-and-fire” terminology sometimes being used broadly for the library.

## Clock and capacitance-unit adapter: an important distinction

`params.py` states that its numerical time constants are **seconds**, and DT=.002s. With Gm numerically1, `build_network.py` passes those second-valued numbers directly into the SNS field named `membrane_capacitance`. For example, the MN field is numerically0.03 and the intended functionalτ is30ms. Under physical nA/mV/µS units and a **seconds-based** derivative, that numericC has unitsµF:0.03µF=30nF. Calling the number0.03 **nF** would imply0.03ms at1µS and introduce a1000-fold mismatch.

The SNS-Simscape block uses an explicit nF parameter and the derivative multiplier `1000/Cm` while Simulink integrates in seconds (`sns_build_library.m:84–89`; spiking equivalent in `sns_build_spiking.m:109`). Thus the equivalent30ms/Gm1 cell needs **Cm30nF**, not0.03nF, in that block. The spiking-knee example Cm2nF/Gm0.4µS correctly corresponds to5ms.

Native source-kernel imports use their source TimeConst in **milliseconds** and explicitly pass `dt_seconds*1000` to the membrane/gate kernel (`rat134_core.py` plus `replay_steps`; W2L `spiking_central_runtime.py:11`). This is a different but explicit clock adapter. A source TimeConst/relativeSize representation should not be relabeled as a measured physical capacitance without its normalization.

These observations distinguish an internally correct functional time constant from an unqualified biological-capacitance label. This audit makes **no numerical changes**. A transfer ledger must include clock units, physical versus normalizedC/G, absolute versus rest-relative reversal/threshold values, and source activation-range units; copying like-named fields alone is insufficient.

## Bandwidth: what the paper supports

The supplied2020paper does **not** prescribe one preferred locomotor control bandwidth or a universal ratio between spike rate and controller bandwidth. Fmax is an encoding limit, τmem is the membrane-response scale, τθ controls threshold-transient duration, and τs trades average-conductance linearity against fluctuation. The page11 example uses100Hz, R20mV, θ01mV and δ1%, givingτs2.171ms; those are illustrative design choices, not mandatory robot values. Section3.5 and Figures7–10 show spike timing/finite-population artifacts; §5.3 discusses when spiking versus graded models are useful.

For the **graded single-cell leak-only approximation**, `fc=1/(2πτ)` is a useful engineering calculation, not a claimed paper prescription: the current30msMN and20msafferent filters imply approximately5.3Hz and8.0Hz. PF80ms implies2.0Hz before phase-ranked multipliers. With synaptic conductance, `τeffective=C/(Gm+Σg)` changes; NaP gates, delays, muscle/tendon/BPA dynamics, and the closed mechanical loop introduce further poles/nonlinearities. A desired whole-model bandwidth must therefore be chosen from the task, measured closed-loop response, actuator limits, and robustness tests, not read offFmax or gait cadence alone.

## Additional syn6 versus Di Russo2023 comparison

[Di Russo2023](https://doi.org/10.1088/1741-2552/acfdcc), §2.2.2/Fig3, uses **five prescribed raised-cosine temporal primitives**, equally spaced at normalized phases0.1,0.3,0.5,0.7,0.9 with common half-width0.2. A contact-reset phase-oscillator pair supplies their phase. Weighted primitive inputs can excite or inhibit MNs. This is not an NMF of the currently simulated activation matrix.

Actualsyn6 instead obtains **six nonnegative spatial-muscle factors** and temporal coefficients from Gait2392 backsolved activations, then assigns their PF channels to RG-E/F families and staggered membrane filters. Positive PF→MN weights use2017Eq18; inhibition resides elsewhere in the circuit. Therefore “same synergy controller” or “Di Russo reproduction” would be inaccurate. Di Russo's three-stage imitation→dynamic/effort optimization→remove imitation is a useful training precedent, not a description of completed syn6training. Its published mechanical demonstration is sagittal;3D extension needs explicit frontal/transverse/trunk control and validation.

The next concrete implementation task, if authorized, is a **scoped2020design ledger for one population pathway**, with actual input domain, desired transient, selectedFmax/R/θ0, derivedbias/τmem/τs/Gmax, and forward transient/rate checks. Completing that on one pathway should not be relabeled whole-controller completion.
