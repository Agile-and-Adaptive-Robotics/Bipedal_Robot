# FSA activation-to-circuit backsolve

Source: `bsolve_out.npz['acts'] (converted-MuJoCo ridge/NNLS)`; 121 frames, 0.500-2.500 s.

Selected **6 NMF synergies per leg**. Per-side minimum ranks at centered VAF >= 90%: right 5, left 6.

These S1--S6 labels are NMF components, not PF neurons and not the existing ZCode E1/E2/F1/F2 early/late stance/swing channels. The direct S-to-MN FSA fit is an exploratory reduced mapping, not a claim of one PF neuron or one PF layer per synergy. The cycle plot maps heel strike to 0%, measured toe-off to 50%, and the next heel strike to 100%. Only one complete measured stride per side is available in this recording, so no between-cycle variance can yet be estimated.

The target mapping is `V_MN = 5 mV * activation`. Analytical conductances use Szczecinski et al. (2017) Eq. 18; the dynamic fit uses the implemented conductance-based LIF equation with nonnegative excitatory PF synapses and a nonnegative tonic bias.

## Side r

- NMF centered VAF: 0.9261
- Dynamic FSA centered VAF: 0.8626
- Dynamic FSA activation RMSE: 0.1251
- Analytical one-synapse gains outside Eexc/R limit: 0
- Mean frames needing net negative MN current: 0.161

## Side l

- NMF centered VAF: 0.9126
- Dynamic FSA centered VAF: 0.8350
- Dynamic FSA activation RMSE: 0.1332
- Analytical one-synapse gains outside Eexc/R limit: 0
- Mean frames needing net negative MN current: 0.168

## Interpretation

The upstream PF excitation/inhibition traces are exact inverse requirements for the PF leaky membranes. They are not evidence that the current two-half-center RG generates those waveforms. That forward closure must be tested separately.

References: Szczecinski, Hunt & Quinn (2017), DOI 10.3389/fnbot.2017.00037; Szczecinski, Quinn & Hunt (2020), DOI 10.3389/fnbot.2020.577804.