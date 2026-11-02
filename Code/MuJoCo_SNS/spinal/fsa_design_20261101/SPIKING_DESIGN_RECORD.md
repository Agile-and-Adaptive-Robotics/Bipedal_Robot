# Analytical spiking design and actual Gait input ports

Machine: **easteregg2**. This is an electrical design record. The all-family registry remains authoritative: **no Hunt training, curriculum, optimizer, or physical trials are released**. Existing controller sources, saved standing winners, K6 synergy files, and accepted central-spiking checkpoints are preserved.

## Scope and equations

`table2_glif.py` implements the seven rules in Szczecinski et al. (2020), Table 2, with explicit seconds/Hz and mV/nA/µS/nF conversions. Membrane capacitance is `1000 G_m tau_mem` nF. Reset is voltage reset; adaptive threshold remains a separate state. Synaptic spikes **reset** conductance to its maximum, followed by exponential decay. Persistent-Na/Ca bursters are explicitly outside the pure GLIF derivation.

The existing Gait numerical capacitance convention is µF with seconds: C=0.05, G=1 means 50 ms and 50 nF. Treating it as 0.05 nF would introduce a 1000-fold clock error. Carrier clocks in the new ports use physical seconds, independently of any intrinsic rhythm scale.

`conductance_matched_design.py` is a separately labelled derivation, not a replacement attribution to Table 2. It solves exact constant-threshold Eq. 13 using total postsynaptic conductance and G×E to match firing rate; finite-δ Eq. 28 then supplies event conductance. It currently covers excitatory, constant-threshold receivers. Adaptive/inhibitory rate-function derivations remain unimplemented. Graded receivers instead receive the intended voltage-equilibrium conductance budget.

## Completed numerical checks

External results root: `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101`.

| Check | Result | Evidence subfolder |
|---|---|---|
| Literal input/rate mapping, R20 mV, Fmax100 Hz | Measured rates within 0.5 Hz of exact Eq. 13 | `gait_table2` |
| Adaptive threshold, requested 500 ms pole | Fitted 493.247 ms; threshold RMSE 0.000924 mV | `gait_table2` |
| Literal Table-2 unit-gain transmission | **Fails** registered 3% rate-gain tolerance: 1.12148 at N16 | `gait_table2` |
| Full-conductance rate counterpart | N16 gains 0.99749, 1.00375, 0.99374; two-input max relative error 0.438% | `full_conductance_glif` |
| Gait-like R5/θ1/E40 counterpart | Gains 1.01587, 0.99492, 1.00314 | `full_conductance_glif` |
| Scalar versus accelerated electrical state | Max 7.11e−15; actual events exact | `full_conductance_glif` |

The literal negative is preserved as an approximation/matching-scope finding, not a claim that the paper is wrong. Halving the clock retained the mismatch. A too-coarse clock was rejected before the Gait-like test; its rejection is retained.

## Actual compiled-family sensory subvariants

`event_conductance_port.py` copies the final intended compiled network, removes selected graded presynaptic axons on the copy, and delivers actual event-reset/decay conductance to existing graded receivers. Pure input cells use exact within-interval GLIF event times. No original graded input voltage carries the replacement output; poisoning that diagnostic voltage leaves all other outputs unchanged. The underlying RG, PF, MN, and other sensory pathways retain their declared graded/ionic implementations.

The three-afferent tests replace bilateral soleus Ia, II, and Ib inputs: six populations ×16 cells = **96 actual spiking cells**. The explicit carrier choice is 500 Hz, R20 mV, θ1 mV, m0, δ0.01. Each physical input is declared in **0–1 nA**; values outside that range raise an error instead of being clamped. This domain does not cover every possible physical-run afferent current, so these results do not establish full-family design eligibility.

| Actual family | Converted outgoing axons | Remaining graded cells | Net right soleus MN ΔIa / ΔII / minimum ΔIb (mV) |
|---|---:|---:|---|
| s3k | 14 | 404 | +0.162650 / +0.162650 / −0.225221 |
| w2lvar | 20 | 882 | +0.057870 / +0.056715 / −0.161065 |
| syn6 | 12 | 788 | +0.042542 / +0.041510 / −0.135840 |

Each 0.8 nA, 100 ms pulse produced **640 actual spikes**. Outgoing-off controls retained those spikes but removed all effects on nonencoded neurons exactly. All three families passed 128-step disabled parity, 128-step poisoned-shadow parity, and 32-step complete-state restoration. Their actual delivered-conductance phase lag at 7.1556 Hz was −0.91586°, with normalized RMSE 0.009959. Results are in `actual_gait_three_afferent_ports\{family}\RESULTS.json`; earlier Ia-only checks are separately retained under `actual_gait_ports`.

These are net responses through the actual competing network routes, not isolated sign proofs for a single branch. For example, w2lvar Ib first reaches RG_E, while syn6 Ib first reaches IBIN. A compiled MuJoCo model supplied the 92 actuator names, but **no body integration occurred**.

Opt-in `install_port_profile` consumes a complete profile through the actual `SpinalNetwork.compiled`/`step` API; unknown keys are rejected. Default existing builders and standing configurations are unchanged. A consumption fingerprint detects changed source parameters, event conductances, or policy. Reintroduced graded axons are rejected directly in the forward path.

The actual public class is `SpinalNetwork`; the original `actual_SpinalCord_step64_exact` result key denotes its `step` interface. Additional `API_CLOCK_GUARDS.json` checks confirm that installed interface matches direct backend execution exactly for64steps. Changed event peaks, reintroduced graded axons, and unused profile fields are rejected. Partitioning the same physical waveform into2ms versus1ms intervals gives event-time error below9e−16s and interval-mean conductance error below4.3e−14µS. This is carrier-only clock invariance, not whole-network integration parity. Full-input firing measured499.921875Hz; zero input produced no spikes.

The encoder is not lossless at arbitrarily small amplitudes. Its exact constant-threshold model is silent through θ/(2R)=2.5% of the declared input domain. For the guarded zero-origin, unsaturated source activation, elementary logarithm bounds give steady normalized conductance error between−2.5% and+δ/(1−δ)=1.0101% of full scale over the whole0–1 domain. A1001-point independent evaluation lies between−2.5% and+0.5791%. This is a DC bound; low-amplitude transient bandwidth is not established. Nonzero source activation thresholds or saturation inside the declared domain are rejected rather than silently encoded linearly.

## Legacy production spiking selector

The actual `AARL_NET=spiking` builder is audited separately in `production_spiking_mirror\PRODUCTION_MIRROR_APPLICABILITY.json`. A fresh build contains422 cells:310spiking plus112graded, including92gradedMN and12gradedreadout sinks; there are1040spike edges and158graded edges. The physical2ms host interval contains four500µs substeps. This changes numerical resolution, not biological cadence.

Its298 nonadapting LIF cells and12 post-spike-adapting RG/PF cells are not silently assigned Table2 compliance. Source bias remains16nA; conversion to a reset-centered voltage coordinate still leaves the effective bias different from the analytical Table2 bias. Membrane constants copy graded values; synapses add capped increments with5/20ms decay; nominal40/50Hz calibration points are not a declared maximum firing domain. The retained graded sink taps do not carry central presynaptic chemical output, which comes from actual spikes. Source files are unchanged, and full-family FSA verification remains false. No neural or body time integration was used for this audit.

## Accepted central-spiking anchor is a different design

The accepted standing anchor remains `D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\Gait_dynamic_readout_v4\TIMING_CANDIDATE.json`, with its current central/physical checkpoints, approved `BURST_METRIC_V2.json`, and 82-array fresh proof. It is **256 central spiking cells plus 778 graded cells**, using TimingZeroSFA and retained NaP/observer mechanisms. The 8 s protocol contains a 0.6 s pinned preparation and 7.4 s released standing under retained apparatus. No walking claim follows.

`accepted_v4\ACCEPTED_V4_APPLICABILITY.json` and `.md` audit that actual anchor. Its fixed thresholds, refractory time, heterogeneity, NaP gates, event-history observer, additive rate filters, and retained source weights are not silently relabelled as all seven Table-2 rules. **FSA design verification remains false** for this central variant; the new sensory subvariants do not inherit its physical pass or confer a new one.

Non-GLIF classification is not a functional failure. Companion `accepted_v4\ACCEPTED_V4_FUNCTIONAL_CONTRACT.json/.md` documents the appropriate ionic/observer contract, with existing scoped core evidence passing. This includes the independently derived passage-time inverse, actual event-marked conductance/channel histories, explicit support regularization, physical clock, source-normalized chemical conversion, approved waveform/phase proofs, and complete-state/fresh controls. Remaining requirements concern new operating domains/cadences, adoption and verification of local reflex budgets for this exact variant, and fresh confirmation of changed profiles after global release. The retained NaP architecture does not need fictitious GLIF threshold-adaptation rules to qualify its own declared functions.

## Reproduction

Use `D:\Anaconda\envs\myo\python.exe` on easteregg2. All scripts below are under `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101`.

- `verify_table2.py --run-design-probes --output <external-directory>`
- `verify_conductance_match.py --output <external-directory>`
- `verify_gait_event_ports.py --family s3k --all-afferents --output <external-directory>`; use fresh independent processes for w2lvar and syn6.
- `audit_accepted_v4.py` audits frozen central provenance without a physical run.
- `verify_event_contract.py --family syn6 --output <external-directory>` checks the installed public API, consumption guards, physical-clock partition, and steady DC bound; use fresh independent processes for the other families.
- `audit_production_mirror.py --output <external-directory>` audits the actual legacy spiking selector without time integration.

Remaining design work includes full declared physical afferent domains, adaptive/ionic compatibility, all required family pathway proofs, and independent central-spiking design. Passing this selected input-port subset does not release the global gate.
