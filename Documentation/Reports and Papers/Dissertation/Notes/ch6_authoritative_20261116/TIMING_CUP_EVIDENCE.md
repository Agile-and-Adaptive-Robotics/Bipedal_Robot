# Independent timing, cup, and knee evidence review

Reviewed on easteregg2, 16 November 2026. Both canonical Chapter 6 and Appendix D were read in full. This review reads saved recordings and configurations; it does not rerun a scientific simulation, training, or benchmark. No canonical manuscript or model source was edited by this reviewer.

`TIMING_audit.py` reproduces the numeric checks with `D:\Anaconda\python.exe`. `TIMING_RESULTS.json` records results, exact input hashes, source bindings, array shapes, equality checks, and timing identities. The audit took approximately 2 seconds for saved-data computation; process startup is additional. Scratch extraction `tmp/ch6_authoritative_20261116/TIMING_trace_arrays.npz` supplies unchanged numeric samples for figure generation and can be recreated by the audit.

## Material findings and revision ledger

| Existing statement | Required correction or clarification | Independent source and scoring scope |
|---|---|---|
| Gait2392 process timing 199.624/29.284 s, without controller/workload/repeat identity | Retain values and 6.82 ratio, but identify the 1,025-cell source-setting Hunt-family controller, native SNS forward versus sparse preallocated forward, one paired 16 s protocol. This is a distinct configuration from the current BiW2Lvar+ air candidate and later 1,033-cell compact pair. | `D:\Users\Ben Bolen\NeuroCampaign_20261102_9h\configurations\source_{native,baseline}_s1.{json,result.json}` and `timings\source_{native,baseline}_s1\RESULT.json`. External process times 199.6235725999577/29.284279999905266 s, ratio 6.81674852858. Protocol is 16 s; first 0.6 s is held-pose preparation. Captured model physics clock spans 0.002–15.400000000001812 s while logical sample clock spans 0.602–16.000 s. No repeated throughput estimate. |
| Neural-call timing 179.079/9.780 s, factor 18.31 | Retain, explicitly tied to that same 1,025-cell source-setting pair. | `compiled.neural_step_wall_s` in the above result JSONs: 179.0791120993672/9.779516393900849 s, ratio 18.3116531418. This sums time within instrumented neural calls, rather than source-process time. |
| Compact full diagnostic run called “complete-process” in Appendix D | Replace scope with **evaluator diagnostic wall time**. Retain 92.944/91.516 s and ratio 1.0156. The saved field is not an external source-process timer. | `D:\Users\Ben Bolen\NeuroAcceleration_20261104\body\{sparse,compact}_final.result.json`: 92.94409799994901/91.51630390016362 s. Same 1,033-cell stage-4 Hunt-PT+RC condition, 16 s protocol, one pair. Saved evaluator source hash `ecda3d37012e593b54ece1166cb9df440c9f961acbc172689db258535bf6ac33` exactly recovered from git `9d7b8d4d50d6e7434eb2bd9743bf160578231de3` with CRLF restoration. Timer begins immediately before `R.main`; it includes capture saving, mechanics saving, and assessment, but excludes interpreter import/pre-run setup and final result/sidecar serialization. |
| Compact kernel 0.895287/0.555062 s, factor 1.613, with only “1,033-cell circuit” identity | Retain, identify Hunt-family implementation and current sparse baseline. Median of three independently instantiated 3,000-update blocks after 20 warmup calls each; this is not whole-body timing. | `NeuroAcceleration_20261104\results\KERNEL_BENCHMARK_hunt.json` and `source\benchmark.py`: sparse median 0.8952872001100332, compact median 0.5550623000599444 s; ratio 1.6129490329524199. Independent state comparison uses 1,200 updates with nonzero, zero, and None inputs and mutable-edge zero/restore checks. Saved maximum voltage/gate/output errors are zero. This audit validates the stored repeats/median, not rerun kernel performance. |
| Historical W2L isolated step gain 1.425 without workload/repeat count | Retain only with bilateralrg2 W2L core and Ib-port configuration, dense SNS versus ExactFastForward, 10,000 timed updates per backend in a single measurement block; 5,000 electrical parity updates precede timing. | `Code/MuJoCo_SNS/spinal/reports_20261027_codex/benchmark_exact.{py,json}`: 145.75979000073858/102.29317999910563 microseconds per update, ratio 1.424921876531876. The separate 30 s body replay pair matches all 16 saved numeric arrays exactly; this is not the same timing workload. |
| Appendix cup: spiking sag 2.02 versus graded 2.00, compared with 7.85 without reflex, without a scoring-window distinction | Prefer a matched **10 s endpoint** comparison: spiking 1.78, graded 2.00, disabled 7.85 degrees, at 0.625 kg added mass. Preserve the original 2.02 as the spiking maximum over `t>2` through 10 s, not an endpoint. The spiking maximum occurs at 9.688429131176422 s, mass 0.6055268206985264 kg. | Actual MCOS numeric samples in `Code/MuJoCo_SNS/spinal/reports_spiking_20261002/goal2_beer_spiking.mat` and `Code/Matlab/SNS_Simscape/results/sns_beer_results.mat`. Spiking endpoint 1.7775121577393576 degrees; maximum after 2 s 2.0242912575529264 degrees. Graded endpoint 1.9965350504348054; maximum after 2 s 1.9973946758572858 at 9.821378199012235 s. Disabled endpoints 7.846412871381035 (spiking file) and 7.846412871422803 (graded file). |
| Cup source comments say 0.5 kg over 8 s | Retain chapter's 10 s/0.625 kg clarification, and apply it to the appendix spiking comparison and any caption discussing the endpoint. | Both raw files end at exactly 10 s and 0.625 kg. Saved mass arrays satisfy `mass=0.0625*t` exactly. `sns_build_beer_demo.m` sets StopTime 10, unclamped Ramp slope pourRate=0.0625 kg/s. Spiking builder copies the same elbow model and load schedule. |
| Knee fine-step maximum difference 5.21e-5 degrees | Correct rounding to 5.20e-5 degrees if retained. | Raw `P` and `F` arrays in `results/knee_transient_reflex_20261028/replay_summary.mat`: maximum angle difference 5.202734295076041e-5 degrees, matching saved summary. No model trajectory changed. |

## Independent saved-array agreement

The source native/sparse pair has **45 numeric arrays across six files**, all exactly equal including NaN masks. Its compiled `model.mjb` files have the same SHA256, `ccb41ca3ef88ef85ea90d68dd39fc7e34713bfbaa83f32d2eac3a73500c6e9db`.

The compact body pair has **62 numeric arrays across seven files**, all exactly equal including NaN masks. Its compiled-model SHA256 is `6cd0dca1adddb64be3f7711abda9a500b4ec5a254d081ae747997664c183476a`. Model identities differ between these two pairs; they must not be merged into a single workload. The later compact result retains a 1.0450 ratio within accumulated neural-call time (16.906603014096618/16.178931983420625 s), also distinct from its warm-kernel median.

The older W2L reference/ExactFastForward body files `trial458_d30.npz` and `trial458_exact_fast_d30.npz` match in all 16 numeric arrays. The nonnumeric `fall_t` object entry is excluded from this numeric audit without enabling pickle. The archived 1.425 timing record has only one timed block per backend and does not supply a repeat distribution.

Other family warm-kernel records remain separate: s3k 410 cells, w2lvar 888 cells, and syn6 794 cells. The saved sparse/compact ratios are 1.082, 1.434, and 1.476, respectively. These secondary measurements need not be added to the manuscript to establish the narrower reported result.

## Raw cup extraction and interpretation

The MAT files store MATLAB timeseries objects in an MCOS subsystem. Standard SciPy loads the objects as opaque references, so `TIMING_audit.py` reads the embedded MAT property table and resolves the stored timeseries IDs in their recorded order. Each extracted time/data pair has equal sample count, starts at zero, ends at 10 s, and retains every stored sample. No interpolation, resampling, smoothing, amplitude scaling, sign change, or decimation is applied. Both independently saved spiking `final` and `sagAfter2` scalars agree with recomputed statistics to 1e-12 degrees. The graded endpoint values match the laptop's separately MATLAB-exported `graphics_metrics.json`.

The spiking enabled recording has 27,088 samples, compared with 2,427 for the graded enabled recording; the differing variable-step clocks are retained. Matching endpoint statistics require no alignment interpolation. The disabled recordings have 19,644 and 657 samples, respectively. The arm's raw angle is positive sag. Plotting negative elevation is a display transformation and must remain explicit.

The knee and cup use different actuator laws. `sns_build_demo.m` connects the reduced knee to `SNS_Library/BPAForce`. `sns_build_library.m` describes the simplified law as `F=Fmax*A*max(0,epsMax-strain)`, with a zero-to-Fmax saturation. The cup builder connects `BPA_20mm` blocks, maps normalized activation to pressure with Pmax=620 kPa, and retains an activation time constant of 150 ms. That first-order command state is not a measured valve-pressure model.

## Knee transient recomputation

The existing `replay_summary.mat` contains 50,001 points from 0 to 50 s on the previously generated 1 ms interpolation grid. This audit does not interpolate again. The original `run_transient.m` locates a 2 nA pulse downstream of the spindle transducer, upstream of extensor Ia-neuron input, from 30.000 to 30.030 s. Two MN output ramp widths change from 1 to 10 mV, and four Ia synaptic ramp widths change from 0.5 to 10 mV. The physical/neural source configuration otherwise remains the stated engineering control.

Recomputed values:

| Statistic | Recomputed value | Window/definition |
|---|---:|---|
| Paper-sign baseline | -42.639228100789865 degrees | Mean control angle, 25≤t<30 s, with negative sign conversion |
| Peak extension versus matched control | 1.0907928305291834 degrees | Maximum after onset, t≥30 s |
| Peak latency | 0.11100000000000065 s | Relative to 30 s onset |
| Sustained recovery time | 6.355000000000004 s | Last threshold violation after pulse offset followed by displacement≤0.01 degree and speed≤1e-4 rad/s through 50 s |
| Late mean residual | 9.30109875646625e-7 degrees | Stimulated angle relative to baseline, 45≤t≤50 s, paper sign |
| Maximum fine-step angle discrepancy | 5.202734295076041e-5 degrees | Entire 0–50 s, existing P versus F arrays |

The altered ramp widths, pulse duration, and simplified force law are numerical functional-test conditions, not physiological calibration. The raw late residual depends on the declared baseline window; it must not be replaced with a single endpoint difference from the control.

## Verification attempts and boundaries

An initial read-only MATLAB extraction command was attempted:

`& 'D:\Program Files\MATLAB\R2025a\bin\matlab.exe' -batch "addpath('Documentation/Reports and Papers/Dissertation/Notes/ch6_authoritative_20261116'); TIMING_export_mat"`

It failed before script execution because the default preferences folder was not writable. A retry with `MATLAB_PREFDIR=D:\GitHub\Bipedal_Robot\tmp\ch6_authoritative_20261116\matlab_prefs_timing` passed that boundary but failed with **MathWorks Licensing Error 63 / -63.2**, “You are not a license administrator.” No model was loaded or simulated, and no numerical result was generated by MATLAB in this session. The successful raw-data verification instead uses SciPy under the base Anaconda Python, as documented above.

The saved current source `physical_eval.py` differs from the compact experiment's locked source. That discrepancy was resolved for timer-scope verification by recovering the exact historical bytes from git, restoring CRLF, and checking the saved SHA256 before reading timer placement. No historical source was installed over current code. The full-process source benchmark has its own external wrapper timing records, independent of inner evaluator `wall_s`.

No new scientific simulation, optimizer, retraining, or timing benchmark was performed. No model points, transforms, stiffnesses, bounds, flags, source recordings, or adopted controller weights were changed. Root editor owns manuscript revision, figure placement, full build, and rendered-page review.
