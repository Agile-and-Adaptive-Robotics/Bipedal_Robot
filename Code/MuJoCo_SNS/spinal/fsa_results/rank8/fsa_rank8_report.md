# FSA backsolve with 8 synergies (rank override)

What-if companion to `fsa_backsolve_report.md` (selected 6 by the smallest centered-VAF>=90% criterion: right 5, left 6). This run fixes the shared synergy count at **8 for both legs**. Source: `bsolve_out.npz['acts'] (converted-MuJoCo ridge/NNLS)`; 121 frames, 0.500-2.500 s, identical inputs to the rank-6 run.

The recorded rank-6 outputs in `fsa_results/` were not modified; every file of this variant lives in `fsa_results/rank8/`.

## Why 8 is not a selection of record

The centered-VAF criterion was already satisfied at rank 6 (right 0.9261, left 0.9126). At rank 8 the curve reaches 0.9528 (right) and 0.9462 (left), so components S7-S8 buy +2.68 / +3.35 percentage points of centered VAF. Only one complete measured stride per side exists, so the extra components are unconstrained between cycles and carry no functional labels. NMF components also carry scaling and permutation ambiguities; the S1-S8 ordering is not identified beyond this run.

## Side r (8 synergies)

- NMF centered VAF: 0.9528
- Dynamic FSA centered VAF: 0.8961
- Dynamic FSA activation RMSE: 0.1088
- Analytical one-synapse gains outside Eexc/R limit: 0
- Mean frames needing net negative MN current: 0.161

## Side l (8 synergies)

- NMF centered VAF: 0.9462
- Dynamic FSA centered VAF: 0.8747
- Dynamic FSA activation RMSE: 0.1161
- Analytical one-synapse gains outside Eexc/R limit: 0
- Mean frames needing net negative MN current: 0.168

## Rank 6 (recorded) vs rank 8 (this run)

| metric | r rank 6 | r rank 8 | l rank 6 | l rank 8 |
| --- | --- | --- | --- | --- |
| NMF centered VAF | 0.9261 | 0.9528 | 0.9126 | 0.9462 |
| Dynamic FSA centered VAF | 0.8626 | 0.8961 | 0.8350 | 0.8747 |
| Dynamic FSA RMSE | 0.1251 | 0.1088 | 0.1332 | 0.1161 |
| Invalid analytic synapses | 0 | 0 | 0 | 0 |
| Negative-current frames | 0.161 | 0.161 | 0.168 | 0.168 |

## Interpretation

The upstream PF excitation/inhibition traces are exact inverse requirements for the PF leaky membranes. They are not evidence that the current two-half-center RG generates those waveforms. That forward closure must be tested separately. S7-S8 have no provisional functional labels, unlike S1-S6. The direct S-to-MN FSA fit remains an exploratory reduced mapping.

References: Szczecinski, Hunt & Quinn (2017), DOI 10.3389/fnbot.2017.00037; Szczecinski, Quinn & Hunt (2020), DOI 10.3389/fnbot.2020.577804.