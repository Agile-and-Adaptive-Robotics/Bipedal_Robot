# Historical implementation and literature-baseline review

Read-only source and saved-array review on easteregg2, 16 November 2026. No neural/mechanical simulation, benchmark, figure generator, training, or model build was run. No canonical source was edited by this reviewer. `HISTORICAL_implementation_check.py`, `HISTORICAL_native_check.py`, and `HISTORICAL_SPIKING_audit.py` reproduce the bounded calculations. Their JSON results record every inspected input path and hash.

## Material findings sent for correction

1. The W2L NumPy/Simulink **0.99996 minimum is best-lag correlation at a -2 ms comparison lag**, not unshifted same-clock correlation. The Simulink model is a **26-cell representative core**, compared with the four corresponding RG signals in the **95-cell full NumPy network**. It is not full-network equivalence. Retain the rounded period and burst counts, with this scope.
2. The Di Russo “native runtime unavailable” statement is obsolete. A subsequent source-native SCONE baseline was run, and its retained fresh-process repeat is independently verified here. That success does not validate the separate rejected 960-trial Gait2392 adaptation. Preserve both results with their chronology and model identity.
3. The AnimatLab bilateral spiking result is **strongly asymmetric firing**, rather than strictly unilateral: the low-rate side remains active. The high-rate side's 442–450 Hz range is supported.
4. The spiking MuJoCo contact fractions and cycle statistics use the post-startup **2<=t<16 s** window. Make this explicit where the table currently says only “post-startup.”

No other material numerical discrepancy was identified in the bounded implementation checks below. Source-level and raw-array evidence are distinguished explicitly.

### Canonical correction readback

The coordinator's final correction batch was read back after this review:

- `60-sns.tex:106` states 1 ms command-source and 5 ms mechanical/sensor rates.
- `60-sns.tex:160` identifies the 26-cell representative core, 95-cell NumPy network, post-startup 2–12 s burst window, and best-lag correlation with 2 ms relative alignment; it excludes complete-network and unshifted sample agreement. The appendix's corresponding passage has the same scope.
- `95-AppendixD.tex:352` replaces the obsolete Di Russo runtime blocker with the separate 10 s native baseline and exact saved fresh repeat, while retaining the rejected adaptation and excluding Gait2392/physiological validation.
- `95-AppendixD.tex:366` specifies the spiking ground table's `2<=t<16 s` contact window.
- `95-AppendixD.tex:410` now describes strongly asymmetric AnimatLab firing.

No material historical implementation finding remains open from this review. Missing historical compiled identities and unavailable paired raw Simulink arrays remain bounded evidence limitations, as recorded below; they are not silently reported as newly performed validation.

## W2L representative-core comparison

Primary files:

- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\campaigns\20260930\w2l_numpy_ref.npz`
- Same directory: `w2l_numpy_ref_20260930.py`, `w2l_numpy_ref.json`, `w2l_trace_compare_20260930.py`, and `w2l_trace_compare.json`.
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\SNS_W2L_CPG_run_20260930.mat`
- Same results directory: `SNS_W2L_CPG.slx`.
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev\build_sns_w2l_cpg_20260930.m`
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\logs\laptop_20260930\sns_w2l.log`

The source builder and saved run log identify four persistent-sodium half-centers plus 22 graded cells, with 42 synapses. The full NumPy reference has 95 cells and 208 synapses. Shared RG definitions and pathways are reproduced in a representative Simulink core; Renshaw cells, Ia/II/Ib chains, second PF pairs, toe contacts, and most MNs are not transferred into this core.

Both source generators prescribe tonic currents E=2/F=3 nA, a 10 nA/10 ms kickoff, and alternating 20 nA heel-contact trains at 1 Hz with active fraction 0.30. The reference JSON's literal protocol string “no contact” is stale: the executable archived generator explicitly supplies those trains. The manuscript's contact-driven description is supported by the executable source, not that incorrect metadata label. No historical source file was changed in this review.

The raw arrays independently reproduce all four stored correlations, RMSE values, lag values, and burst counts exactly:

| RG signal | Best-lag correlation | Lag (ms) | Unshifted correlation | Bursts Simulink/NumPy |
|---|---:|---:|---:|---:|
| Left extensor | 0.9999715821150922 | -2 | 0.9994734734775296 | 10/10 |
| Left flexor | 0.9999701745385268 | -2 | 0.9995913481609914 | 10/10 |
| Right extensor | 0.9999619662968479 | -2 | 0.9994315224019197 | 10/10 |
| Right flexor | 0.9999754082731023 | -2 | 0.9995541353808316 | 10/10 |

The burst gates use native t>=2 s, with Simulink samples through 12.000 s and NumPy through 11.998 s. The cross-correlation script truncates to 5,000 matched 2 ms samples over 2–11.998 s and searches lags within ±200 ms. Burst periods are 1.000000 or 0.9997777778 s, all rounding to 1.000 s. At the best lag, the correlation uses the overlapping samples; no scientific recording is phase shifted or overwritten.

The no-heel-train latch is documented in the archived source comments and earlier campaign records. A separate retained raw no-contact trace was not identified in these primary comparison files. The present check does not rerun that deletion. The saved SLX is retained, but no separately hashed historical compiled Simulink binary or fully frozen NumPy network object was found; source and saved signal agreement have that boundary.

## MuJoCo–Simulink bridge

Primary files:

- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge\logs\bridge_tests_abc.mat`
- Same directory: `gt_const.mat`, `gt_const.npz`, `test_abc.log`, and `test_abc_laptop.log`.
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge\matlab\run_bridge_tests.m`
- Same MATLAB directory: `bridge_groundtruth.py`, `e1_sns_reflex_loop.m`, and `e2_cpg_open_loop.m`.
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge\matlab\exp_out\e2_cpg.mat`
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\mujoco_bridge\logs\e2_cpg.log` and `e12_laptop.log`.

The five signals are right-knee position, right-knee velocity, and right-vastus-medialis actuator length, velocity, and force. The Simulink record contains 81 rows at 5 ms over 0–0.4 s. Excluding its initial pre-step row leaves **80 samples over 0.005–0.400 s**, matching the Python ground truth exactly for all five numeric arrays: every recomputed maximum absolute error is zero. The saved pass flags for tests a, b, and c are true.

The two-rate test uses a **1 ms command-source interval**, a rate-transition block, and **5 ms MuJoCo mechanics/sensor intervals**. Its 81-row grid and 5 ms mean interval independently recompute. These rates should replace the methods' vague “different stated rate combination” when that sentence is clarified.

The separate E2 source sends the full exported network's 92 MN command channels through the bridge during a 5 s fixed-step run at 2 ms. Its primary log reports all commands finite and 96% of channels recruited. Its saved `e2_cpg.mat` uses MATLAB opaque timeseries and is hash-preserved; this review does not reconstruct all those object arrays. Source and log establish the bounded signal-exchange demonstration. The same log contains an instability warning, so it must not be described as mechanically accepted walking or a clean physiological transfer. The current main text confines it to interface operation and separately reports the tonic/fixed-point limitation, which is appropriate.

No historical compiled MEX/MuJoCo binary digest was bound to these saved test arrays in the primary results. The source-built bridge's installation reports remain historical records; bit-exact saved readback does not establish the present machine's ability to rebuild or execute that bridge.

## Deng and the early 410-cell export

Primary Deng files:

- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\demos\sns_run_deng_demo.m`
- Same directory: `SNS_Deng_CPGDemo.slx` and its RG/PF builders.
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\results\deng_cpg_ref.mat`
- `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\deng_cpg_ode.py`
- `D:\GitHub\Bipedal_Robot\Code\Matlab\SNS_Simscape\logs\sns_run_deng2.log`
- `logs\desktop_capability_20261024_REPORT.md` records the subsequent rebuilt verification.

The raw NumPy clock is **milliseconds**, 0–20,000 ms. Converting that clock to seconds and applying the burst analysis yields ten bursts after 0.3 s and mean period **1.937866666666667 s**, supporting the rounded 1.94 s reference. The Simulink primary log reports ten bursts, **1.938 s** period, and extensor/flexor correlation **-0.857**, supporting -0.86. Those Simulink pointwise arrays are not retained in that log package, so the period/correlation are source-bound saved verdicts rather than freshly recalculated Simulink statistics.

The old `sns_run_deng_demo.m` pointwise reference interpolation compares a seconds clock directly with the exported millisecond clock. Its printed maximum voltage differences therefore do not represent a correctly aligned waveform comparison. The manuscript does not quote those differences and appropriately limits its conclusion to period/alternation and differing configurations. No pointwise-equivalence claim should be added from that output.

For the 410-cell export, the source `sns_verify_from_json.m` explicitly scores all cells at 0.3 s and nine traced cells throughout t<=0.3 s with a 2 ms Euler update. The primary `logs\desktop_capability_20261024_b.log` reports **4.201e-6 mV** maximum endpoint error and **4.198e-6 mV** traced error. The older `logs\sns_verify.log` reports failure over a 2 s comparison, with later trace errors of order millivolts. This supports the manuscript's early agreement/later divergence distinction. `results\verify_ref.mat` preserves the NumPy reference and the SLX export is retained, but the paired full Simulink verification arrays were not saved, so this audit cannot independently recalculate the 4.201e-6 verdict. No full-time or full-controller equivalence follows.

## Di Russo source-native baseline versus failed adaptation

The stale native-runtime blocker is explicitly cleared by:

`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261027_codex\native_baseline\NATIVE_BASELINE_20261027.md`

Primary native record directory contains `eval_run.log`, `eval_repeat.log`, `native_build_manifest.txt`, `summarize_native.py`, `native_baseline_summary.json`, and `eval_spinal_native.sto.sto`. The corresponding external native workspace is:

`D:\GitHub\_zcode_build_tools\dirusso_full`

This workspace preserves `eval_repeat.sto.sto`, `spinal_controller_main.par`, its source scenario, initial state, and the original reduced planar **9-DoF/18-muscle** model. It is a different body/runtime from the exploratory Gait2392 adaptation.

`HISTORICAL_native_check.py` independently reads the native and repeat arrays: **1,001 rows × 1,338 columns**, 0–10 s at 10 ms storage. The numerical arrays agree exactly, and every text line after the output-filename header is byte-identical. The stored source-model, optimized parameter file, and initial-state MD5 values match the build manifest exactly. The native log reports objective 17.2411 and 16 steps; raw trajectory progression is 11.5743 m, its 1–10 s fitted speed is 1.1740958918 m/s, and pelvis height is 0.855665–0.918060 m. These numbers need not be added to the manuscript to correct the chronology.

The native scenario still has nonzero gait and joint-limit penalties. It is evidence of a repeatable author-native source baseline, not physiological validation or validation of the Gait2392 port. The rejected 960-trial adaptation remains rejected. A concise source-grounded correction should state that the custom runtime was subsequently built and an independently repeated native baseline became available, then retain the failed adaptation as a separate outcome.

The build manifest binds author-repository revision `ccf46fad3e6651a8daaca11d3967b88ab556ae12` and the custom SCONE fork `58ea68762791e133b0b18f2f0db83e3eb2a18fe3`, plus component versions and executable MD5s. The surviving input files' hashes were verified; the historical runtime executable was not located or rerun in this review. A read-only attempt to ask git for the external workspace's current HEAD failed because changing to that external working directory was denied. This does not affect direct saved-array or input-file reads. The historical revision binding remains the manifest's declaration, rather than a fresh repository-state verification.

## Historical spiking comparisons

The separate `HISTORICAL_SPIKING_REVIEW.md` and `HISTORICAL_SPIKING_RESULTS.json` provide exact primary paths and hashes. The delegated reviewer independently recalculated all retained MuJoCo table values from `spinal_run_spkbase_s3k.npz`, `spinal_run_spkbase_spiking.npz`, the exact reference MOTs, and the scorer. Scores -161.5675417368/-236.3496680669, right/left cycles 9/0 versus 8/9, right mean-cycle knee minima -64.2921/-79.8273 degrees, contact fractions, and pelvis tilts reproduce. Contact/cycle scoring uses 2<=t<16 s, while maximum tilt uses the whole stored trial.

The saved neural-only record reproduces 0.39922449/0.62106667 s periods and **55.5683%** increase, with the original spiking 100 ms count-envelope analysis. The rounded manuscript values remain valid. The raw AnimatLab charts and archived native inputs support the displayed ranges and negative conversion outcomes; only “unilateral” requires the narrower “strongly asymmetric” wording. Missing historical compiled-model/network digests are documented by that reviewer, rather than replaced by present-source hashes.

## Verification boundary

The bounded review establishes saved-array statistics, saved repeat agreement, and source-bound historical verdicts. It does not claim new simulator validation, a current dynamic replay of evolving software, universal implementation equivalence, physiological validation, or unsupported walking by the dissertation's SNS controllers. Exact SHA256 records for this review are in `HISTORICAL_IMPLEMENTATION_RESULTS.json`, `HISTORICAL_NATIVE_RESULTS.json`, and `HISTORICAL_SPIKING_RESULTS.json`.
