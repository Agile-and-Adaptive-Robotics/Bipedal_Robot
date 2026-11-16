# PF timing: independent local evidence review

Reviewed on easteregg2, 16 November 2026. Canonical Chapter 6 and Appendix D were read before this review. This audit ran saved-array calculations only. It did not run a neural model, mechanical simulation, optimization, or training campaign, and did not modify manuscript sources or scientific recordings.

## Finding and required interpretation change

The original syn6 timing failure is real and reproducible. The manuscript's subsequent blanket statement that output fitting and body work remained held is obsolete. Later experiments established timing under progressively different replay conditions and then exposed a persistent live-body acceptance failure. These are useful subsequent results, not a replacement of the original failure by an unrelated success.

The supported-ground W2L108 and CAD examples remain different controller/body families. Neither promotes these syn6 experiments or the BiW2Lvar+ air candidate through a continuous physical curriculum.

## Verified sequence

All relative raw paths below are under `D:\Users\Ben Bolen\ControllerTraining_20261104\`.

| Configuration / evidence | Independent result | Appropriate claim |
|---|---|---|
| Original `waveform_edge_D_20261104\syn6\BEST_TRACE.npz` | Recomputed 0/12 channel passes; three centered-variance components per side explain 99%, versus six for intended primitives; bilateral S5 duties 1.000 versus 0.325/0.330 targets | A downstream MN fit concealed PF timing failure in this checkpoint |
| `s3k_relay_output_20261104\BEST_TRACE.npz` | Recomputed rank 2 per side; quantitative channel windows not declared | Rank describes the recorded PF outputs and does not certify an intended design |
| `syn6_PF_input_repair_mirrored_20261105\BEST_TRACE.npz` | Recomputed 0/12 passes | Gain-only repair remained unsuccessful |
| `syn6_PF_relays_20261105\combined\RESULT.json` and `FRESH_REPEAT.json` | Retained reports identify all 12 outputs passing in offline D replay, with RG trace and MN synapses preserved | Hybrid RG-event-triggered relay proof of timing concept; these reports do not contain a raw full-waveform array, so report preservation is weaker than later raw-array verification |
| Fixed native SNS relays, `syn6_native_20261105_v7\final_verified\` | Recomputed 48/48 channel-cycle checks at nominal 1.233 s; 16/48 at 1.048 s and 8/48 at 1.418 s; all saved replay arrays match fresh repeats exactly | Nominal neural timing works, but cadence generalization fails |
| Adaptive native clocks, `syn6_adaptive_20261105\replay\` and `replay_live_range\` | Recomputed 48/48 checks at each of six controlled replay cadences, approximately 0.863, 0.925, 0.986, 1.048, 1.233, and 1.418 s | Adaptive timing works on these replay afferents/cadences; this is not six body-stage passes |
| Adaptive live body, `syn6_adaptive_20261105\body_slow_handover\` and `live_full_native\TIMING.json` | Predictive knee guard at 9.814 s; 11/30 channel-cycle checks after full native handover; fresh captured arrays exactly equal | Transmitted native PF behavior does not retain its replay timing or pass the coupled-body stage |
| Adaptive PF plus 150 ms projected-knee premotor relay, `syn6_flexor_handover_20261105\predict150_stage2\` | Completes 16 s, but stored assessment rejected; 12 R/11 L knee cycles, period 0.796 s, CV 0.287, phase error 0.138 cycles; native timing recomputed 16/48; max native knee coordinates R/L -1.404/-1.452 deg; fresh repeat arrays exact | Knee-safe completion is useful progress, not accepted rhythmic air motion or PF timing |
| Same candidate `predict150_stage3\` and `predict150_stage4\` | Completes 8 s standing and 16 s supported-ground diagnostics, respectively; both stored assessments rejected; saved compiled identities verified | Completion is not standing or ground-stepping acceptance |
| Later `syn6_stage2_regularize_20261105\` | Retained records: S2 latch yields one S2 onset in all 11 counted intervals, but S3 still opens repeatedly and full score remains 16/48 | A bounded local repair did not solve overall PF timing |
| Latest `syn6_s3_single_20261105\{target,overlap,wide,late}\STOP.json` | Guards at 6.450/8.050/9.084/6.890 s for requested windows 0.26–0.41/0.20–0.65/0.16–0.88/0.40–0.88 cycle | No reviewed candidate meets single-burst timing and knee safety together |

The row for regularization is report-checked; its 16/48 raw rescoring was not repeated in this audit. The four last guard times were read directly from STOP files. The original, fixed-native, adaptive replay, adaptive live, and 150 ms live timing numbers above were independently recalculated from raw arrays.

## Scoring definitions and windows

### Original failure

`source_snapshot_pf_design_20261105\pf_timing_audit.py` scores transmission `clip(V / 5 mV, 0, 1)`, with active transmission **greater than** 0.2. It interpolates each channel at 200 equally spaced phase samples in each of cycles numbered 2–10 inclusive, using the original 1.233 s calibration clock. The sampled window is 2.466–13.556835 s and the nine cycle waveforms are averaged. There is no per-channel phase realignment. The composed syn6 target mirrors the right primitive onto the left with the measured contralateral loading phase offset; it does not use independently factorized left primitive labels.

The original combined gate requires absolute fundamental phase lag ≤30 deg, absolute duty difference ≤0.15, matching onset/offset counts, and transmission excursion ≥0.2. These are provisional engineering requirements. Original recomputed scalar features match the retained audit exactly (maximum difference zero). `PF_original_audit.py` reproduces this calculation and saves `PF_ORIGINAL_RESULTS.json`.

### Later relay/native timing

`syn6_s5_timing.py::measure` adds an explicit maximum onset and offset error of **5 percentage points** of cycle to the fundamental-lag, duty, and excursion requirements. This stricter timing definition must not be presented as an unchanged 0/12-to-48/48 denominator. `native_pf_check.py::assess` identifies ipsilateral RG-E onsets with a Schmitt trigger armed below 1.5 mV and triggered at ≥3.5 mV. It scores each of the last four complete RG intervals separately, at 200 phase samples per interval. Six channels × two sides × four intervals gives 48 channel-cycle checks. The reference origin is locked to the nominal reference RG phase, not independently adjusted for each channel. Each exact interval is preserved in `PF_AUDIT_RESULTS.json:offline[].windows`.

For the fixed-native nominal replay, four periods on each side alternate between 1.232 and 1.234 s. The shortened replay has 1.048 s periods; the lengthened replay has approximately 1.418 s periods. The cadence protocol scales legacy-cell integration and RG inactivation constants together with afferent replay cadence, while keeping native relay parameters and MN synapses fixed. Thus these are controlled neural cadence stress tests, not naturally emerging cadence changes in a freely moving model.

Live scoring first excludes samples before full native handover, recorded at neural time 6.241999999999535 s. In the longest adaptive capture, the remaining native interval spans 4.172 s and contains two complete left and three complete right RG intervals: 4/12 and 7/18 checks, respectively. In the 150 ms candidate, 9.758 s of native capture is available; the last four complete intervals per side produce 8/24 passes on each side. Neural time includes preparation that differs from reported experiment time, so the handover neural timestamp must not be substituted for the 9.814 s body guard time.

The saved body assessments for the 150 ms candidate are all rejected, and this audit independently verifies finite saved arrays, native knee maxima, raw repeat identity, and compiled-model hashes. It does not claim an independent reimplementation of all historical body acceptance gates, whose running evaluator has since changed. The retained source hashes in each recording identify that evaluator version.

## Controller, model, and source identities

The native repair adds graded SNS onset, shutoff, quieting, trigger, follower, and ramp cells. The adaptive implementation observes recurrence of the four bilateral E/F RG signals, adjusts clock integration rates, and hands over from legacy timing. It therefore combines SNS computation with an explicit timing scheduler. The 150 ms knee relay senses projected extension and recruits S3 while inhibiting S2; it does not change the frozen MN input synapses. These are new circuit/timing designs, not merely a gain-only retry of the original implementation.

All reviewed native replay and live TIMING records bind the frozen motor-input digest:

`41fbe663292846161c97852f2e2b0fdc6c30daaac231193f1648c439466c3a1e`.

For the 150 ms body tests, SOURCE.json binds:

- Reflex checkpoint: `waveform_edge_D_20261104\syn6`.
- Peak checkpoint: `waveform_edge_fit_20261104\syn6`.
- Timing checkpoint: `commissural_fit_v4\syn6\training\B`.
- Shared air expression 1.0, different from BiW2Lvar+ expression 0.1.
- No sustained direct RG current; F/E command interfaces remain; POST/BAL/FUNC channels are clamped off. Live Ia, II, and Ib return remains enabled.
- The repaired Gait2392 native biological-muscle body and its separate suspension remain. These trials are not CAD reference tensions or dynamically calibrated BPA actuation.

`PF_AUDIT_RESULTS.json:source_files` records hashes of all raw arrays, timing records, model binaries, configuration files, scorer sources, primitive arrays, and reference motion used in the recalculation. Its body rows verify the actual MJB bytes against each recording's declared model SHA256. No assumption that the present runner matches a historical runner was required.

Original hashes independently verified:

| File | SHA256 |
|---|---|
| Original s3k BEST_TRACE | `89f58384f9dc85d40866af8505518b16a9aa774414eac80cfd9216657fd73aea` |
| Original syn6 BEST_TRACE | `8c2bcf9a0f265a0ec9020765d4c3689d49e68b2a2377b4bc2f6ed578a3b19251` |
| Mirrored repair BEST_TRACE | `a83ec6920cc3a4f37a0c1b0cb5c9d027662456c7cc35072c7bc96a5f589f5efd` |
| Original retained scorer | `94e7686757f7b02283e6ed0e3a635456313ac9abe6cea61c66891981b338d81c` |
| FSA primitive arrays | `45519a6b166a6182c84eb262891ce21217c999fa298f74279c1220ab1517f470` |

## Revision ledger recommendations

| Existing statement | Recommended replacement/addition | Reason |
|---|---|---|
| Chapter 6: further downstream fitting held pending PF design; Appendix D: new fitting/body launches held | Original gain-only repair failed. Subsequent RG-triggered relay and native-SNS timing experiments passed specified offline checks, followed by unsuccessful live-body tests. | Historical hold is not current experimental status. |
| Original syn6 0/12, s3k ranks, saturated S5 | Retain unchanged and bind to original named checkpoints and original clock. | Raw-array audit reproduces them. |
| No later timing/body evidence | Add compact progression: fixed-native nominal 48/48; adaptive 48/48 at six controlled cadences; live failure; knee-safe 16 s candidate still 16/48 and rejected. | Distinguishes offline function, transmitted behavior, and physical acceptance. |
| Original failure appears to be final state | Add latest single-S3 tests could not satisfy timing and knee safety together. | Preserves unsuccessful follow-up and prevents favorable replay from becoming locomotion acceptance. |

Suggested main-chapter wording: “Subsequent RG-triggered relay designs reproduced the intended PF windows in offline replay. A native SNS implementation passed all 48 channel-cycle checks at the nominal cadence, and adaptive clocks extended that result to six replay cadences. The improvement did not transfer directly to the live mechanical loop: coupled trials either reached the predictive knee guard or completed without satisfying the retained rhythm and PF-timing requirements. The original failure and the later partial repairs therefore identify separate verification boundaries, from intermediate neural function to accepted coupled behavior.” Appendix D can provide the exact counts, criteria, replay scope, and negative live results above.

## Remaining limitations

- Hybrid-relay success is corroborated by retained numerical reports and fresh-repeat reports; no full raw hybrid-relay waveform was found in that directory. The stronger native replay claims do have raw arrays and independent rescoring.
- The post-5-November running body/evaluator sources have changes. This audit uses archived raw arrays and saved version bindings; it does not claim a fresh dynamic replay under current code.
- The reviewed later single-S3 attempts end at the knee guard. No later accepted syn6 body-stage result was found in the specified authoritative local training root.
- No figure was created or edited by this PF audit. The original failure figure remains scientifically valid when its checkpoint and offline scope are retained.

## Reproduction

Run from the repository with the verified environment:

```powershell
& 'D:\Anaconda\envs\myo\python.exe' 'Documentation\Reports and Papers\Dissertation\Notes\ch6_authoritative_20261116\PF_original_audit.py'
& 'D:\Anaconda\envs\myo\python.exe' 'Documentation\Reports and Papers\Dissertation\Notes\ch6_authoritative_20261116\PF_audit.py'
```

Both scripts read preserved recordings and write only their own JSON audit results in this Notes directory.
