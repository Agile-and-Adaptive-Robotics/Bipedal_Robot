# Canonical newer W2L: historical LH style and reference-parameter checkpoint

**Session 1, easteregg2, 31 October 2026.** This checkpoint precedes the newly requested expansion from four to twelve type-II channels. The MuJoCo agent owns that next functional edit; do not regenerate this checkpoint over its later work.

## Editable model and standalone

- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.aproj`
- `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer_Standalone.asim`
- Matching runtime: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\W2L_VirtualWalker_Newer.asim`

W2L here means the newer bilateral-RG/contact/type-II model. No original 80-neuron model, V3-zero candidate, 1 nA candidate, or earlier trained result was substituted. The physical model, current 96 neurons, all 204 chemical endpoint pairs, commissural paths, contacts, and existing type-II channels are retained at this checkpoint.

## What was actually copied

Reference: historical uniarticular LH of `Biped_2xCPG_wSubs.aproj`, git `23b0e785`, preserved at `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_20261031\build\subnetwork_rebuild_20261031\23b0e785.aproj`.

The correction transplants **148 actual node geometries** and **156 actual link drawings**, including their saved control points and ArrowDst/ArrowOrg properties. Matching right-side local motifs use the same LH geometry. RG/PF pages retain the reference spacing instead of the previous enlarged grid. The K&A overview again presents subsystem icons; PF-to-MN links use destination MN-page OffPage connectors as in the reference. The corrected PNG/SVG previews follow stored drawing head styles, including filled inhibitory circles. These are reconstructed previews of saved drawing data, not screenshots of a verified GUI opening; physical-body icon images are not reproduced.

**92 matched local neural edges use 19 exact historical named synapse definitions**, including equilibrium potential, maximum conductance, presynaptic threshold, and saturation. This is a numerical model change expressly requested after the first geometry-only pass. The full per-edge/type before-and-after ledger is `REFERENCE_PARAMETER_CHANGE_LEDGER.json` in this folder. For example, the canonical PF-to-MN default type was 5 uS and is now the historical 2.749 uS. Other reference edges receive their historical pathway-specific types.

The remaining model-specific newer connections are not falsely described as historical copies. Their endpoints, link payloads, and effective type numbers are verified unchanged. Their drawings use explicit control points and historical OffPage/class styling, with additional pathways on overview/contact/commissural pages. No historical clone folder was created.

## Exact local type values at this checkpoint

| Historical named type | Equilibrium (mV) | Maximum (uS) | Threshold (mV) | Saturation (mV) |
|---|---:|---:|---:|---:|
| Hip_Ib_ext_Excite | -10 | 0.59 | -60 | -40 |
| PF to MN Excite | -10 | 2.749 | -60 | -50 |
| RG Inhibit | -70 | 2.749 | -60 | -25 |
| Knee KEbE-5 | -10 | 0.589 | -60 | -40 |
| RG_HIP_EX | -40 | 0.176 | -60 | -40 |
| PF Inhibit | -70 | 2.749 | -60 | -25 |
| RG_K&A_EX | -50 | 0.333 | -60 | -50 |
| AnkleZ AEbE-5 | -10 | 0.5 | -60 | -40 |
| Hip_Ib_flx_Excite | -10 | 1 | -60 | -40 |
| RE to MN Inhibit | -70 | 0.5 | -60 | -40 |
| Knee KFbE-6 | -10 | 0.5 | -60 | -40 |
| AnkleZ AFbE-6 | -10 | 0.5 | -60 | -40 |
| RG Excite | -40 | 2.749 | -60 | -25 |
| Knee Flx MN Excite | -10 | 1.516 | -60 | -50 |
| Knee Ext MN Excite | -10 | 4.93 | -60 | -50 |
| RE Inhibit | -70 | 0.5 | -60 | -40 |
| PF Excite | -40 | 2.749 | -60 | -25 |
| RG to RG Excite | -40 | 0.1 | -65 | -40 |
| MN to RE Excite | -40 | 0.5 | -100 | -10 |

The full GUI type chooser inventory, including added newer types, is `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\subnetwork_review\NAMED_SYNAPSE_TYPES.csv`.

## Validation and actual limits

All 18 embedded drawing pages parse and have verified object counts and interleaved AddFlow endpoint handles. The 204 GUI chemical edges exactly match the synchronized standalone IDs/endpoints/types. All exported intrinsic and used type numerical fields match, polynomial adapter coefficients match, and native object IDs are unique. Earlier standalone-only edges were removed and GUI-only edges exported; the complete difference ledger is `CANONICAL_EXPORT_SYNC.json`.

A fresh native **one-second suspended/deafferented probe** completes with exit 0 and 496 finite samples for all 96 neurons, spanning **−98.541 to −44.9764 mV**. No RG-flexor threshold crossings after 0.1 s were observed in this short probe. It is a load/finite-response test, **not a rhythm, air-stage, ground-stepping, or curriculum pass**. The first probe also completed but its chart did not flush because chart and simulation ended together; the recording endpoint was corrected in the probe only and the test repeated.

Canonical contact polynomial Actual C=1 A/contact, D=0 and the inherited left/right RG-extensor tonic 2/0 nA remain as supplied. They were not silently converted to a candidate 1 nA/20 nA contact model or symmetrized. Therefore the successful earlier assisted simulations cannot be assigned to this newly copied-reference-parameter configuration. GUI adapter metadata versus drawn adapter targets also has a separately flagged discrepancy; no unnoticed sensory rewiring was made in this style pass.

Actual AnimatLab GUI opening remains unverified: earlier CUA access failed, and the user's existing unsaved window was not touched. No owned native process remains running.

## Review and provenance paths

- Comparison: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\subnetwork_review\LH_reference_comparison.png` (and `.svg`).
- Every saved page: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\subnetwork_review\exact_*.png` / `.svg`.
- Verification: `D:\GitHub\Bipedal_Robot\Neuromechanical_Models\W2L_VirtualWalker_Newer_20261031\REFERENCE_STYLE_PARAMETER_VERIFICATION.json`.
- Before-edit bytes and source hashes: `D:\GitHub\Bipedal_Robot\tmp\w2l_lh_style_20261031\before\`.
- Before numerical reference-type changes: `D:\GitHub\Bipedal_Robot\tmp\w2l_lh_style_20261031\before_reference_parameters\`.
- Native probe and raw 96-neuron data: `D:\GitHub\Bipedal_Robot\tmp\w2l_lh_style_20261031\native_reference_smoke\`.
- Exact style/route correspondence: `D:\GitHub\Bipedal_Robot\tmp\w2l_lh_style_20261031\STYLE_REPAIR_AUDIT.json`.

Checkpoint hashes (subsequent twelve-II augmentation will change these):

- Project: `9761a1fa848c500c6b5597143c2dbdaf547b508ad61d5aebec6b34b9dd80baa4`.
- Standalone: `85318b4d3e32c64256b7966f2615750e9195ccd244335dd0f577bb1223644ef5`.

No commits or pushes. Older SUBNETWORK_REBUILD/SUBNETWORK_VERIFICATION files describe the superseded first layout pass; this record supersedes those layout claims. Do not rerun the old blanket layout builder or generated-variant builders over the canonical project.
