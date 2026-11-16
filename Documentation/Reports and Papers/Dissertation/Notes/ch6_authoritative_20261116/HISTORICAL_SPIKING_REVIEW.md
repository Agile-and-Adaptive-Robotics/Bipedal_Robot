# Historical spiking-comparison evidence review

Read-only review on easteregg2, 16 November 2026. No model was integrated, no neural network was evaluated, no simulation or fitting was run, and no figure generator or canonical manuscript was edited. `HISTORICAL_SPIKING_audit.py` reproduces the saved-array reanalysis with `D:\Anaconda\python.exe`; all inspected file SHA256 values and computed statistics are in `HISTORICAL_SPIKING_RESULTS.json`.

## Findings

The retained MuJoCo graded/spiking table and period comparison reproduce from the saved arrays, reference trajectories, and scoring source. No numerical replacement is required.

One narrow wording correction was communicated immediately: the AnimatLab bilateral variant should be described as **strongly asymmetric firing**, rather than strictly “unilateral firing.” The v1 saved record has right RG extensor/flexor rates 442.4/448.2 Hz and left rates 6.2/9.0 Hz. The rounded 442–450 Hz high-rate range is supported, but the low-rate side is not completely silent. The v2 conversion is a separate trial with high-rate flexor activity on both sides; its result must not be merged with v1 to imply a single unilateral trace.

## MuJoCo neural rhythm

Primary recording:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002\goal4_gate2_traces.npz`

SHA256 `6fdba0137ca402084ca43085f4dee3aff45dedcc7aff3c5efd6dbf3a9d093f9d`.

The 0–20 s saved neural-only comparison has constant 2.5 nA drive. Reapplying the stated analysis from `tools/check_selfsustain_spiking.py` gives:

| Statistic | Recomputed value |
|---|---:|
| Graded mean period | 0.3992244897959184 s |
| Spiking mean period | 0.6210666666666668 s |
| Relative period increase | 55.56827863544969% |

The graded readout is right RG-E minus RG-F voltage, with peak prominence 0.3 mV. The spiking period uses the difference of cumulative spike counts, converted to the original 100 ms count envelope; peaks use prominence 1 and minimum distance 150 ms. Means use interpeak intervals over the full 20 s, while the sustain criterion checks the last 10 s. This is the existing scoring filter, not newly smoothed simulation evidence. Rounded 0.399/0.621 s and +55.6% are correct.

This array file contains the recorded readout and cumulative counts, not a frozen compiled SNS object or complete launch-time source manifest. Current builder hashes cannot retrospectively establish historical compiled-network identity. The report and original analysis script bind the named constant-drive experiment and its conventions, while the reanalysis verifies the stored signal statistics.

## MuJoCo body comparison

Primary recordings:

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\spinal_run_spkbase_s3k.npz`, SHA256 `6127692c35e2dc5382a6df52cc49bb534f1a36e9e00841417440402fedacb75d`.
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\spinal_run_spkbase_spiking.npz`, SHA256 `0d12e7ab6b16ee0b166ba2da8ec1e1cfec033abb61faf32c12c88e95bfc388e3`.

Both store 8,000 samples at 2 ms, `t=0..15.998 s`. The comparison script uses the existing `kine_ref.py` with the exact current reference files `subject01_walk1_ik.mot` and `subject01_walk1_grf.mot`. This is offline scoring of stored arrays; no runner or network builder is invoked. The current scorer SHA256 is `dfda44a8bfb54fd2a5f95eb227291306aac4c04a4c428b5f56226912c825417a`. Its recomputed results agree exactly with the historical primary baseline/stretch logs, rather than merely approximately reproducing the rounded manuscript table.

| Metric | Graded s3k | Untuned spiking mirror | Definition |
|---|---:|---:|---|
| Kinematic score | -161.56754173676563 | -236.3496680668847 | Existing bilateral scorer and reference motions |
| Cycles R/L | 9/0 | 8/9 | Post-startup t≥2 s; own contact onsets, RG fallback only if contact onsets insufficient |
| Right mean-cycle knee minimum | -64.29205155222364° | -79.827318008436° | Minimum of the right leg's phase-resampled mean cycle, not whole-run joint minimum |
| Contact fraction R/L | 0.6367142857142857/0.9847142857142858 | 0.7988571428571428/0.8557142857142858 | Samples with foot force>20 N over 2≤t<16 s |
| Maximum absolute pelvis tilt | 22.0975375424625° | 19.22026766617569° | Full stored 16 s trial |

The primary logs are `reports_spiking_20261002/logs/baseline_s3k.log` and `logs/ground_eval_spiking.log`; `tools/baseline_eval.py` identifies the winner-loading flow, retained stage-4 settings, and untuned mapped spiking implementation. Both logs report finite states and no fall. The current table correctly treats the left graded leg's large contact fraction as compatible with zero detected cycles. It does not interpret sustained contact as successful stepping.

The numeric NPZs and their configuration payloads remain available. No historical compiled MJB, compiled-network snapshot, or complete launch-time source digest was identified for this October comparison. Recomputed saved-array scores and preserved run logs establish the reported bounded comparison; they do not prove that today's evolving runner would dynamically reproduce it. The manuscript's “historical” and “same selected parameters” scope is therefore appropriate.

## AnimatLab conversion table

Primary records are under:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_spiking_20261002\`

The audited inputs include `goal3_animatlab_spiking.md`, `tools/analyze.py`, `tools/make_spiking.py`, `tools/sameleg_xcorr.py`, seven `metrics_*.json` files, their named `.asim` copies, and available raw chart files. Exact source, metric, `.asim`, and raw chart hashes are retained in the audit JSON. The `.asim` copies are archived native-run inputs; no new native run was made.

| Manuscript row | Primary record | Result |
|---|---|---|
| Two-layer walker | `metrics_baseline_W2L_modern.json`; `metrics_W2L_modern_spiking.json` | Baseline RG 0.4432 s; v1 spiking 2 extensor/0 flexor spikes over the 0–5 s neural chart and no sustained rhythm. Left hip 6.146–6.939° over 0.5–5 s supports “pinned near +6.” |
| Bilateral variant | `metrics_baseline_BipRG_modern.json`; `metrics_BipRG_modern_spiking.json` | Graded channel periods 0.4614–0.4654 s and folded lag -0.424 cycle. V1 right RG rates 442.4/448.2 Hz, left 6.2/9.0 Hz: strongly asymmetric high-rate latch, not exclusively unilateral firing. |
| Reorganized standalone export | `metrics_Biped_standalone_stale_spiking.json` | 58 bursts per left half-center, median interval 0.0848 s (11.793 Hz), same-leg envelope r=1.000 at zero lag. This is explicitly the archived standalone export, not a new export of the later GUI project. |
| Thresholds-only ablation | `metrics_bisect_W2L_neuronsonly.json` | Raw physical chart ranges independently recompute: hip -15.1150..23.1531°, knee -0.8100..47.8966°, matching -15.1..23.2 and -0.8..47.9 rounding. Graded synapses remain; the result isolates the tested conversion mapping rather than demonstrating a fully calibrated spiking network. |

All available joint-angle ranges in the seven metric cases were independently recalculated from raw chart values, using the declared radians-to-degrees conversion and the original trailing-zero-fill removal. The original L Angles chart spans 0.5–5 s; L Hip Physical spans 0–5 s. Raw RG spikes/bursts were also remeasured with the original thresholds and refractory/gap definitions; results are retained in the audit JSON. The run stop time is 5.1 s, while scored charts end at 5 s. Those are distinct durations.

The native `.asim` files retain model/neuron/synapse definitions, but no installed-simulator binary digest or machine-compiled model identity was identified for the archived runs. This review verifies source-bound native verdicts and available raw recordings, not a fresh simulator validation. The table's negative conclusion is appropriately limited to the tested conversions; it does not establish that a recalibrated spiking realization could never reproduce the graded rhythm.
