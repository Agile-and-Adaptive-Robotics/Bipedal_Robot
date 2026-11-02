# Physical-domain sensory event ports

Machine: **easteregg2**. Electrical design only; no Hunt, fitting, curriculum, body integration, or standing/walking promotion. Earlier1nA proofs and the accepted centralv4 controller are unchanged.

The runner clips normalized velocity/length to±1 and normalized muscle force to0–1.5. At caller afferent scale1, gamma0, and vestibular proprioception0, the declared canonical input maxima are Ia2.88nA, II2.92nA, and Ib1.125nA. Ib includes the runner's factor0.5; zero force gives0nA. The requested normal Ib operating point0.225nA is distinct from zero force.

`fast_event_port.py` adds an opt-in preallocated exact-event/reset-conductance backend. Numba and scalar versions use the same equations, without fastmath. Default storage is counters and a bounded last-event block; full event history requires explicit capture. Pure-input source axons are replaced on a copy, while the other neurons remain graded/ionic. Source windows with nonzero lower thresholds or saturation inside the declared domain, including the new0.225nA Ib route window, are rejected. They require a separately derived nonlinear axon mechanism.

All raw evidence is under `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\physical_domain_v2`.

## Preserved comparisons

- **Homogeneous5nA domain:** R20mV, θ0.1mV, F500Hz, N16, δ0.01. State/reference parity reached2.22e−16 and event times were exact. Normal Ib onset after steady quiet was22.863ms, so a small DC dead zone did not establish fast low-current response. Evidence: `homogeneous\syn6`. Its source version is separately preserved in `homogeneous\source_snapshot`.
- **Midpoint-bias5nA domain:** deterministic `b_j=Gθ(j+0.5)/N`, with quiet equilibrium `U_j=b_j/G`. Every bias is below threshold; zero input is silent. Mean bias equals the Table2 bias, but this is an engineered population extension rather than literal per-cell Table2 bias. N16 normal Ib onset became1.391ms, and conductance phase passed10°. However, normal Ib passive-receiver voltage error was0.025155mV versus graded modulation RMS0.055245mV, or **45.53%**. This remains a failed dynamic-fidelity comparison. Evidence: `midpoint_population\syn6`, with its separate source snapshot.

The modulation-relative criterion was declared after observing the5nA/N16 shot-noise problem, before the typed/N128 extension: receiver voltage RMSE≤5% of the matched graded receiver's modulation RMS and extra phase≤10°. This does not rewrite the earlier fullscale-conductance gate. Receiver probes integrate actual source passive C/G/reversal parameters at a declared2.5mV operating bias, with other network currents excluded; they do not clamp a model or demonstrate whole-controller standing.

## Canonical typed/N128 extension

`canonical_typed_profile` uses the declared per-type maxima and N128, with the same θ0.1mV, R20mV, population-mean500Hz target, and δ0.01. Per-cell endpoint firing rates are separately derived; synaptic decay uses the actual maximum. The count-resolution choice gives a one-event fraction `1/(128×500×0.020)=0.00078125`; against a0.1nA modulation in the1.125nA Ib domain, this is0.879%. This count calculation is **not** a receiver voltage-error bound.

Changing encoding range does not strengthen the neural pathway. For each axon, the typed endpoint target conductance is the original source conductance multiplied by `typed_current_max/(source_Gm×source_window_max)`. Thus the source physical-current slope is retained. Reset peaks are normalized using the exact population endpoint rate distribution: `Gpeak=Gtarget/mean[tau_s f_j(1−exp(−1/(tau_s f_j)))]`. The carrier uses physical seconds independently of rhythm/cadence clocks.

The δ0.01 comparison did not pass all families: normal Ib receiver errors were5.746% for s3k,8.016% for w2lvar, and4.952% for syn6. Saved-data decomposition showed a dominant mean offset, consistent with the finite-tail gain approaching1/(1−δ) away from full firing rate. These comparisons remain under `canonical_typed_N128`.

## Fixed tighter-tail design

`canonical_tight_tail_profile` declares δ0.001, giving a0.1001% finite-tail gain budget instead of1.0101%. Decay time and reset peaks are rederived from the exact per-cell rate distribution. N128, bias distribution, membrane constants, physical input slopes, and original source conductances are unchanged. This is a fixed analytical design comparison, not a gain search.

All three fresh-process electrical cases passed the predeclared5%-of-modulation-RMS and10° extra-phase criteria at7.1556Hz:

| Family | Ia receiver relative RMS | II receiver relative RMS | Ib receiver relative RMS |
|---|---:|---:|---:|
| s3k | 1.541% | 1.546% | 2.099% |
| w2lvar | 1.541% | 1.677% | 2.001% |
| syn6 | 1.541% | 1.677% | 2.412% |

Ia and II use1±0.5nA; Ib uses0.225±0.1nA. The high-domain Ib wave also passes. Receiver absolute errors, modulation RMS denominators, and phase differences are retained in each result. The reference is the same actual passive membrane driven by continuous conductance proportional to held physical current; it excludes the presynaptic graded encoder pole and other network currents. This is not a full-network graded-versus-spiking equivalence claim.

Normal Ib earliest AP after exact quiet equilibrium is39.064µs, independently checked against the threshold-crossing equation; its first host output is at2ms. The previous1.997ms figure was the first *listed* cell-ordered event, not the earliest event; the verifier now takes the minimum timestamp. This is an earliest-population-event statistic, not every neuron's response delay. Actual conductance phase includes host timing, and the receiver comparison measures the transmitted waveform.

Scalar/JIT discrepancy is at most8.88e−16; event timestamps are exact. Disabled160-step, outgoing-off160-step, complete-state32-step, zero-input, canonical-bound rejection-before-mutation, and nonlinear-window rejection checks pass. Default event history remains empty; the last block is bounded. OriginalK6 and source files match their before/after hashes. The selected six input pools contain768 actual spiking cells; RG/PF/MN and other afferents retain their declared graded/ionic implementation.

Results and complete parameter profiles: `canonical_typed_N128_delta0p001\{family}`. The compact checkpoint is `D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\physical_domain_v2\canonical_typed_N128_delta0p001\CHECKPOINT.json`. Each completed case has artifact hashes in `COMPLETED.json`; its NPZ embeds protocol, profile, and result JSON. Failed repeats leave a separate pending protocol instead of replacing a completed one.

Matched scalar-port versus Numba-port **neural-only** benchmarks measured13.62–17.90× in the final runs, excluding warmup. Both use the already accelerated graded-network backend; this is not a comparison against an unmodified SNS-Toolbox backend or a full MuJoCo simulation. No body, standing, walking, broader-input, or robustness pass follows from these local electrical checks. The global training gate remains closed.

## Reproduction

Program: `D:\Anaconda\envs\myo\python.exe`.

Script: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_design_20261101\verify_physical_event_port.py`.

Use `--family syn6 --tight-tail --output D:\Users\Ben Bolen\NeuroCampaign_20261031\fsa_design_20261101\physical_domain_v2\canonical_typed_N128_delta0p001\syn6`, quoting the output path. Run each family in a fresh process. `--typed` selects the preservedδ0.01 comparison; `--population-bias` selects the distinct5nA/N16 midpoint comparison; no switch selects homogeneous. No NMF is run.
