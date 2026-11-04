# Acceleration experiment record — Session 1

Machine: easteregg2. 4 November 2026. CAD simulation remains on hold pending Ben's route revision. No training, tendon/damping changes, manuscript edits, commit, or push in this pass.

## Implemented

Optional CompactPairwiseForward in D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_cpg\compact_pairwise_eval.py. The existing current Hunt path already used SparsePairwiseForward, preallocated I/O, and Numba edge calculation. This increment packs nonempty synaptic reduction groups, reuses state/current/output workspaces, caches constant gate functions, and evaluates ionic force on nonzero ionic channels. The original reduction tree, fixed-tau equations, inactive gate evolution, timestep, and mutable edge weights are preserved. No fastmath or physics changes. This is Numba plus NumPy, not a new C++ or GPU implementation.

It removes unused dense g/ge wrapper buffers. Parent construction still transiently allocates these before removal; compiled SNS objects retain other dense matrices. Some indexed expression temporaries remain. Do not call this allocation-free or fully linear-memory execution. Spiking and electrical channels remain unsupported by this backend; use their existing implementations.

## Fresh-process electrical comparisons

Each family ran in its own fresh process. 1,200 comparison steps include nonzero/zero/None inputs and active edges temporarily set to zero and restored. All voltage, b/c gate, and output states match exactly. Warmed timing is the median of three 3,000-step runs, versus the already accelerated SparsePairwiseForward, not unmodified SNS Toolbox.

| Controller | Cells | Current sparse (s) | Compact (s) | Incremental speedup |
|---|---:|---:|---:|---:|
| hunt | 1033 | 0.895287 | 0.555062 | 1.613x |
| s3k | 410 | 0.353856 | 0.327003 | 1.082x |
| w2lvar | 888 | 0.781493 | 0.545040 | 1.434x |
| syn6 | 794 | 0.664793 | 0.450304 | 1.476x |

Hunt wrapper buffers removed: 17,073,424 bytes (16.28 MiB), excluding additional packing savings. Electrical benchmark files: results\KERNEL_BENCHMARK_{hunt,s3k,w2lvar,syn6}.json. The earlier KERNEL_BENCHMARK.json is exploratory/superseded and must not supply the headline numbers.

## Fresh-process Gait2392 body pair

Final configs: body\sparse_final.json and body\compact_final.json; each evaluates the same saved 16 s stage-4 Hunt-PT+RC condition, without refitting. Both accepted. All 62 numeric arrays across ordinary traces, assistance, balance, mechanics, regional contact, contact events, and aligned neural/mechanical capture are exactly equal. The model.mjb hashes are identical. Current accounting is available for all 7,700 captured samples in each run, with independent voltage-update residual exactly 0 mV. This establishes parity on the tested condition, not all possible networks/states or a new curriculum pass.

Full diagnostic run wall time: sparse 92.944098 s, compact 91.516304 s (1.0156x). Recorded neural-step wall time: 16.906603 s versus 16.178932 s (1.0450x). These are one paired timing observation with capture/diagnostics, not a repeated end-to-end throughput estimate. Warm kernel speedups must not be presented as whole-simulator speedups. An earlier reviewed-profile pair also matched, with 81.706790/81.248882 s wall time; timing varies. Neither pair is a speed claim against bare MuJoCo, AnimatLab, or Simscape.

## Selection and provenance

Set "backend": "compact" in a Hunt physical_eval configuration. It defaults to the previous sparse implementation. physical_eval selects the implementation before FSA install creates mutable edge caches; the original constructor is restored immediately afterward. Do not replace the wrapper after install: zero-but-later-enabled pathway slots would otherwise be lost. Current capture recognizes the actual compact forward binding and checks numerical update closure. Call refresh_parameters() after changing intrinsic channel parameters; edge gmax values remain mutable. Returned buffers are reused, so copy histories.

D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\hunt_training\controller.py is unchanged, SHA256 00112d00d7e65c700c1c89e651f57354588d88335386209e5803b67f7154b57c. Historical integrated profiles/checkpoints were not edited. body\FINAL_BACKEND_PROFILE.json preserves their locks and adds hashes for the new backend/evaluator/capture files. An early profile rejection occurred while an optional selector was placed in controller.py; that approach was removed, the exact original file restored, and the final paired trials repeated successfully. No provenance gate was disabled.

## Files and restart

Root: D:\Users\Ben Bolen\NeuroAcceleration_20261104
- results\VERIFIED_SUMMARY.json: compact final comparison, timing, closure, model hashes.
- source\benchmark.py: fresh-family benchmark; argument hunt, s3k, w2lvar, or syn6.
- source\summarize_verification.py: read-only comparisons; writes compact summary/reviews.
- body\*_final.result.json and matching .npz files: final body outputs. Result JSON includes large compiled matrices; inspect selected fields, do not print it in full.
- body\{sparse_final,compact_final}\capture\: aligned traces, metadata, actual model.mjb, and independent SPARSE_CURRENT_REVIEW.json.
- body\FINAL_BACKEND_PROFILE.json: additional source locks and original fitted derivation.

Run from D:\GitHub\Bipedal_Robot with D:\Anaconda\envs\myo\python.exe. Example: python Code\MuJoCo_SNS\spinal\hunt_training\physical_eval.py <full config path>. Redirect stdout because this runner prints full compiled arrays. The exact final config paths are above; avoid overwriting them when making a new experiment.

Remaining performance work: profile biomechanics, target allocation/linear algebra, capture, serialization, and startup separately. Optimize those measured hotspots next; this pass does not claim completed C++ encoder/Jacobian/solver kernels or faster spiking execution. GPU/Orin hardware performance remains unmeasured.
