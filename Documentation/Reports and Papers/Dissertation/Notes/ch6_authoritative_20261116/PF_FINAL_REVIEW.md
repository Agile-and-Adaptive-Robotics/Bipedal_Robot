# Independent PF evidence, prose, and figure review

Review scope: final PF-related revisions in canonical Chapter 6 and Appendix D, plus the replacement raw-array PF figure. Both manuscripts were read in full at the start of the audit, and their current changes were read back for this review. This reviewer did not edit the canonical TeX. Scientific support was independently recalculated from the authoritative local arrays; no neural/body simulation or fitting was run.

## Verdict

No material PF evidence or prose findings remain in the reviewed text. The original failed configuration, the later replay successes, and unsuccessful coupled-body transfer remain separately identifiable. Later replay success is not presented as successful standing, supported stepping, physiological validation, or harness-free walking.

The following completeness findings were raised and resolved by the coordinating author before this readback:

1. The original nine-cycle averaging clock and 2.466–13.556835 s sampled scoring window are explicit.
2. The new figure caption identifies a raw two-cycle display, distinct from the full averaged audit behind the 0/12 verdict.
3. The script map now includes the actual native timing and checking entrypoints.

## Evidence-to-text check

| Current source | Review result |
|---|---|
| Chapter 6, approximately lines 208–210 | Preserves original syn6 0/12 and s3k variance findings; identifies gain-only failure; separately reports 48 channel-cycle replay checks and negative live transfer. Matches independently recalculated evidence. |
| Chapter 6 Discussion, approximately line 259 | Interprets the replay-to-live gap as a verification boundary. It does not infer biological validity from a fitted output or nominal rhythm. |
| Appendix D, approximately lines 241–254 | Correct outgoing transmission, original unshifted clock, nine-cycle window, active threshold, original criteria, later Schmitt-trigger detection, four-interval scoring, additional edge tolerance, and changed denominator. |
| Appendix D progression table, approximately lines 264–268 | Fixed-native 48/48,16/48,8/48 at stated cadences; adaptive 48/48 at six controlled cadences; adaptive body guard9.814 s and11/30 after handover;150 ms candidate16 s but16/48 and rejected;8 s standing/16 s ground diagnostics rejected. All retain their scopes. |
| Appendix D, approximately line 273 | Three right/two left intervals explain30 checks. MN synapses remain frozen. Later single-S3 failures are retained without unnecessary catalogue detail. |
| Appendix D script map, approximately line 463 | Native timing, measurement, and hybrid relay sources map to existing repository files. |

No claim that the original twelve channel-mean checks and later48 channel-cycle checks have the same statistical denominator is made. Cadence-controlled replay remains distinguished from spontaneous freely moving cadence adaptation. The adaptive timing scheduler is disclosed rather than presented as exclusively neuronal phase integration.

## Replacement figure

`PF_timing_authoritative.pdf` was regenerated from `waveform_edge_D_20261104\syn6\BEST_TRACE.npz` and the raw `fsa_backsolve.npz` primitive arrays. It retains all1,234 saved samples per displayed channel over2.466–4.932 s. The voltage-to-transmission transform is the actual outgoing `clip(V/5mV,0,1)` law. The target uses the original right-primitive time mapping; there is no individual phase adjustment, smoothing, averaging, or decimation.

Visual inspection of the PNG found readable common labels, distinct solid-blue/dashed-black curves, aligned time axes, separated panels, no clipped labels, and a visible0.2 threshold. Independent PDF text extraction with the machine's base Python/PyMuPDF measured:

- Page396.8504×331.6535 pt, corresponding to14.0×11.7 cm.
- Minimum text10.0 pt.
- Panel identities rightS1/rightS2/rightS5, matching caption and raw columns.
- Figure SHA256 `820e98ea7913c2ce632dce6b7229a8a8ad4fa1356e8cbaa7522e49892ab19d32`.

The new caption correctly separates the display window from the nine-cycle scoring audit and identifies the older-readout air candidate as a different condition.

At the initial final-review readback, the TeX include had been updated but the new active PDF had not yet been copied from Notes into `Figures\95-AppendixD\`. That integration action is owned by the coordinator. This review does not claim a complete manuscript build or rendered-page check; those require the final active assets and full-build result.

## Verified native source paths

Absolute paths checked directly with `Test-Path`, `Get-Item`, and SHA256 on easteregg2:

| Source | SHA256 |
|---|---|
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\native_pf.py` | `40dcbdbbff04b89e34505993454b8b1813b9da7b0613d5987e636051cc593839` |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\native_pf_check.py` | `7676dc1bb38aa35210c4986d9e42db6cff0b6980fe1fade44dbb5a956c92ff73` |
| `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\controller_profiles_20261104\syn6_s5_timing.py` | `30753d05fbae0b46e5fb44b8d85eac04ec73ee69a483f2fd74c06a96e36999a4` |

These three files are present in the current repository. No copies were found recursively under `D:\Users\Ben Bolen\ControllerTraining_20261104`. They must not be redirected to invented external snapshot paths. The original PF scorer does have a retained external snapshot, documented separately in `PF_EVIDENCE.md`.

## Remaining verification boundary

This independent review validates the scoped PF text against saved-array recalculation and the replacement figure's scientific content and physical typography. It does not replace the coordinator's full canonical LaTeX build, final page inspection, active-asset move manifest, unrelated evidence audits, or final source-hash conservation check.

Reviewed source hashes before any later coordinator edits outside this PF scope:

- `60-sns.tex`: `394b89d2dbc77783efa013bfbe8c1a6acd4827fa228bb2354258fc38ff4ece6b`.
- `95-AppendixD.tex`: `fa5eb35c130d8d9fce7e5a80902612a2473fc1069ca7854d4519c804b0ea2dc4`.

Reproducible numerical evidence and exact raw hashes are in `PF_ORIGINAL_RESULTS.json`, `PF_AUDIT_RESULTS.json`, and `PF_figure_provenance.json`.

## Active-asset integration closure

After the coordinator installed the figures, the exact active path
`D:\GitHub\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Figures\95-AppendixD\PF_timing_authoritative.pdf`
was rechecked. Its SHA256 equals the reviewed Notes PDF byte-for-byte:
`820e98ea7913c2ce632dce6b7229a8a8ad4fa1356e8cbaa7522e49892ab19d32`.
The active PNG also matches the reviewed Notes preview
(`1fb5f6cc71cd489947155c44ad1954935bad2025d913d2a08f166adebc22d83f`), and its alt-text companion is present.

The current canonical include retains `width=14cm`; the caption still identifies raw right S1/S2/S5, the2.466–4.932 s display, solid/dashed curves, the0.2 threshold, all stored samples, and the distinct full nine-cycle0/12 audit. The figure integration finding is closed. No generator rerun occurred. The full build is in progress at this closure; compiled-page verification is not claimed here.
