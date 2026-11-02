# Syn6 K6/K8 read-only design review

Role: **Gait neural pathway composition reviewer**. Machine: easteregg2.
No refit, basis replacement, body run, physics change, or training was performed.

The reviewed report is
`D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\fsa_results\rank8\fsa_rank8_report.md`.
Its K6/K8 inverse-reconstruction comparison is matched at the dataset and
factor level, but is not a matched forward production-Syn6 experiment.

Array-exact checks establish that production `spinal\synergy_basis.npz`
contains the recorded FSA K6 factors: `W_r/l == r/l_gain.T`,
`H_r/l == r/l_source.T`, and identical muscle-name order, targets, and time.
Rank8 retains exactly the same targets, muscle-name order, and time as K6.
Its gain arrays are8x43 instead of6x43. Both use121 frames at0.500–2.500 s,
from the converted-MuJoCo ridge/NNLS activation backsolve. These are computed
activation targets, not measured EMG. Each side has one complete cycle.

| Metric | Right K6 | Right K8 | Left K6 | Left K8 |
|---|---:|---:|---:|---:|
| NMF centered VAF | 0.9261 | 0.9528 | 0.9126 | 0.9462 |
| Dynamic FSA centered VAF | 0.8626 | 0.8961 | 0.8350 | 0.8747 |
| Dynamic FSA activation RMSE | 0.1251 | 0.1088 | 0.1332 | 0.1161 |
| Fraction requiring negative MN current | 0.161 | 0.161 | 0.168 | 0.168 |

The added components improve reconstruction on this same cycle. They do not
resolve the negative-current requirement, which follows the desired MN
transient rather than NMF rank alone. S7/S8 have no identified functional
labels, and component scaling/permutation are not uniquely identified.

Production `build_network_syn6.py` loads this K6 spatial basis and maps W
through its Eq18 conductance rule. It also adds actual RG/PF/laminated/contact
and local-afferent circuitry. The rank8 report instead evaluates an inverse
PF/MN reconstruction with separately derived/fitted requirements. Its exact
upstream excitation/inhibition traces do not show that the production RG
generates those traces. Forward RG-to-PF closure and closed-loop/body
performance remain untested by that comparison.

The optional canonical-right K6 design is another explicitly different
condition: it copies actual right-side compiled functions to the left without
refactoring or overwriting K6. Thus neither that design nor the running
production controller inherits a K8 performance claim. K6 remains unchanged.

Hashes of read-only sources:

| Source relative to `Code\MuJoCo_SNS\spinal` | SHA256 |
|---|---|
| `synergy_basis.npz` | `63f79c8803b023c9e5ebd5c4c62d571f88d5334df0637058ddc73bfe76c07f6b` |
| `fsa_results/fsa_backsolve.npz` | `45519a6b166a6182c84eb262891ce21217c999fa298f74279c1220ab1517f470` |
| `fsa_results/rank8/fsa_backsolve.npz` | `b952a81293ffe1d721e5f927a6b335fda28e7b82e9f6765f479107e86e3c3878` |
| `fsa_results/rank8/fsa_rank8_report.md` | `1a9a129d026e5b5b4aac2be167caed3272f3f9ae32e262adf87e5ee1bae2f4b0` |

Training requires Ben's explicit approval. Free3D walking is not a prerequisite
for scoped neural design checks; unsupported/control-assisted directions must
instead remain explicitly labeled in each condition's coverage ledger.
