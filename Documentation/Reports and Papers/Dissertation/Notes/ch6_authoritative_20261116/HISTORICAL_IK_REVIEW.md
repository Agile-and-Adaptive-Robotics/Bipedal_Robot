# Historical human-model IK comparison: independent audit

Read-only reanalysis on easteregg2. No simulation, fitting, plot generation, controller change, or canonical manuscript edit was performed. Reproduce with `D:\Anaconda\python.exe HISTORICAL_IK_audit.py`. `HISTORICAL_IK_RESULTS.json` records all results, exact source hashes, loading onsets, retained periods, and comparisons. Calculation took approximately 0.4 s.

## Verdict and recommended prose clarification

The retained numerical claims are correct. All six absolute-angle RMSE values, both cycle counts, both mean periods, both model stance fractions, both reference periods, and both reference stance fractions reproduce **exactly** from saved raw recordings and reference files. The complete model-mean, model-SD, and reference columns in the source CSV also reproduce with maximum difference **0.0**.

The Appendix calls −145.07 a “cycle-shape score.” This is too narrow: `kine_ref.compare` includes waveform shape, range, mean posture, phase, knee minimum, period, period variability, contact, double support, and bilateral terms. Use **“composite kinematic score”** or **“recorded kinematic score.”** The separate six-cell table is correctly described as absolute-angle RMSE between unshifted mean cycles; do not replace it with the scorer's mean-removed shape RMSE.

The source analysis uses **2≤t<16 s** (saved samples through 15.998 s), which should be stated where the historical comparison is introduced. Counts and periods come from qualifying loading-onset pairs within that window. The reported model stance fractions are time fractions over the complete post-2 s window, not means of stance fractions computed only from retained cycles.

## Primary evidence and identities

| Exact path | SHA256 |
|---|---|
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261028_publication_flow\human_ik\human_fresh.npz` | `f5e3ad24ae780af51f947af9388f39a53ebbecf99e649f4cdbb38df6fc617014` |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\finalist_eval.npz` | `223e1b8699639247ec9ad8d08af66107f1ef4f47ef8ff0f23fcba2e033844d03` |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\kine_ref.py` | `dfda44a8bfb54fd2a5f95eb227291306aac4c04a4c428b5f56226912c825417a` |
| `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\subject01_walk1_ik.mot` | `7380251a0e2434bea1ad62aa8ece037e3dc8ff014f8b747395b63c4c8a4d55e2` |
| `D:\GitHub\Bipedal_Robot\Solid_Models\OpenSim\Gait2392_Robotbody\subject01_walk1_grf.mot` | `69a5e8633f8666bc192abd8069e6812fe7b4c8bb8968dc3254d4824abf263662` |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261028_publication_flow\human_ik\gait2392_ik_overlay.csv` | `64925c1d6544b13f5687baf4cd0ae61f9bdb378f88fea7f3f0a6b345fe995d7f` |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20261016\wean_human\round8_best3.json` | `2b9cab27c4ce7e5de496f73f20ff8c461c74771c3f26d4ce85454bdc41acf415` |

Additional original generator: `Documentation\Reports and Papers\Dissertation\Figures\60-sns\Components\gait2392_ik_overlay\publication_ik_overlay_20261028.py`, SHA256 `d135c5cd7ea3d206d246f4c8f655487617c0939d46832d90ec14f5feb6764f3c`. The later active-figure generator `Figures\95-AppendixD\Components\generate_gait2392_ik_overlay.py` reads that CSV directly and changes presentation only; its retained provenance binds the same CSV hash. Both generators were inspected but not executed here.

The raw publication replay is an independent saved recording of the historical round8_best3 configuration at supplementary drive scale **0.3807379786697019**. `human_metrics.json` preserves its selected drive knobs and feedback gains. Its eight numeric arrays (`t`, `act`, `q`, `qfull`, `com`, `neuro`, `contact`, and `footz`) match the original `finalist_eval.npz` **byte for byte**, with equal shapes and dtypes. This independently supports the saved repeat claim without running the model again.

## Recomputed statistics and exact windows

The recording contains samples **0.000–15.998 s** at 2 ms. Analysis begins at **2.000 s**. Loading is defined by normal force **>50 N**, separately for each side. Adjacent loading onsets define candidate cycles; retained durations are **0.35–2.5 s**, with at least eight samples. Each retained cycle is interpolated to the same 100-point 0–100% grid; no fitted phase shift, mean subtraction, or amplitude scaling is applied to the absolute-angle table or overlay.

| Quantity | Right | Left |
|---|---:|---:|
| Retained model cycles | 10 | 10 |
| Mean model cycle period (s) | 1.2848 | 1.2538 |
| Model time in stance, 2≤t<16 s | 0.7101428571 | 0.6485714286 |
| Reference cycle period (s) | 1.2334 | 1.2133 |
| Reference 100-phase-sample stance fraction | 0.61 | 0.61 |
| Absolute hip RMSE (deg) | 28.03326159 | 24.65338202 |
| Absolute knee RMSE (deg) | 27.43993204 | 28.44518307 |
| Absolute ankle RMSE (deg) | 8.71393173 | 7.37790473 |

Reference loading cycles are **0.6283–1.8617 s right** and **1.2550–2.4683 s left**, selected as the first complete loading-onset interval inside the IK recording. The reference stance fraction is the fraction of the 100 normalized phase points above 50 N. The figure's bottom model strip instead shows the fraction of retained cycles in stance at each phase. These are related but distinct summaries.

Recomputing `kine_ref.compare` from `human_fresh.npz` and the named references gives **−145.07280220616852**, identical to `human_metrics.json`; difference **0.0**. This verifies the stated composite score and does not change its experimental scope. Mean periods conceal variable loading onsets: retained period CVs are approximately **0.3222 right** and **0.1707 left** in the same saved calculation. These values are diagnostics, not proposed new headline results.

## Limits

- No compiled-model binary or immutable compiled-model hash was identified beside this historical raw replay. The independently matching raw arrays, recorded controller settings, score, reference files, and generating scripts establish this comparison's numerical identity; they do not establish complete historical mechanical-source recovery.
- Current `kine_ref.py` reproduces every adopted score/overlay statistic exactly, but no separate launch-time scorer-hash manifest was found in this evidence folder. Exact numerical reproduction supports the retained calculation without claiming that all launch-time source files were recovered.
- Support, absence of forward progression, and retained shared descending drive remain separate model conditions; this comparison is not harness-free walking or physiological validation.
- No numerical replacement is recommended. Only the score name and explicit post-2 s window require clarification.

## Closure after canonical correction readback

The affected current Appendix D paragraph now states **composite kinematic score −145.07**, explicitly identifies the **2–16 s post-startup analysis window**, and explains that the score combines several kinematic terms rather than only cycle shape. The retained periods, stance fractions, 50 N loading criterion, ten cycles per side, and absolute-angle RMSE table are unchanged. These amendments resolve the two prose findings of this review without changing the independently reproduced numerical results or raw arrays. This closure is a readback only; no further source search, generator, or simulation was run. Full-build and rendered-manuscript verification remain separate coordinator checks.
