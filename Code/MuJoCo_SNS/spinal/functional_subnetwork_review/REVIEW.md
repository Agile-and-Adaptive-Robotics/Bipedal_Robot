# Functional-subnetwork review: source intake and static validation

Latest drawing revision: `drawings/RG_PF_Feedback_revised.pdf` (10pages). Per Ben's instructions, the updated set explicitly shows RG/PF interneurons, actual Ia/II/Ib/heel/toe central projections, multi-muscle channel counts and member lists, and the triangle base against the receiving object. W2L108's absent direct toe route and the saved syn6 condition's zero contact-input gains remain labelled rather than silently changed. The represented bundle and relay fields were checked against279saved chemical-edge records; all10SVGs parse and all10PDFpages contain selectable text. The earlier PDF remained open/locked, so it was preserved and the revision was saved under the explicit revised filename. No circuit execution or controller changes.

Status: **drawing work in progress; no circuit execution**. Ben clarified that the deliverable is readable connection drawings, especially sensory feedback into RG/PF, using his supplied subnetwork example and established Micklem/Deng motor-layer conventions. The first central-feedback pages are in `drawings/W2L108_L_central_feedback.svg` and the matching R page, with PDF/PNG versions and separate named-synapse keys. `drawings/RG_PF_Feedback.pdf` combines those four pages. Work on EB475WS4, with read-only retrieval from easteregg2 authorized by Ben. No training, simulation, calibration, optimization, model import/compile, manuscript edit, GUI action, or browser-storage operation was performed. Python files under `Gait_TimingZeroSFA/source_snapshot` are unchanged reference copies from the saved ZIP, not edited controllers or executed circuits.

## Completed checks

`STATIC_VALIDATION.json` records the exact saved condition, source hash, census, endpoint checks, and checkpoint hash for s3k, w2lvar, syn6, and W2L704/80. All saved chemical-edge endpoints resolve by both index and name. All 15 checked production-family builder/runner/parameter hashes match their saved values. These are the specific saved **standing candidates** and supported W2L hybrid condition in their ledgers, not a claim that these are the latest winners or that standing completes a walking curriculum.

Each family folder contains lossless node and chemical-edge CSV tables, input tables, static input reachability, a byte-exact editable JSON checkpoint, and runtime context. Nested source parameters are retained as JSON cells. Reachability is graph-theoretic only: it does not establish signal polarity, physiological identity, enabled phase, effective transmission, or causality. The extracted tables are an intake layer, not the completed pathway-by-publication audit requested in the handoff.

The source-manifest check has three failing records, all referring to the same shared `Gait2392_standing_families/PROVENANCE_CORRECTION.md` hash. The graph and figure hashes checked against these manifests pass. The repeated note mismatch is not three independently changed graphs. The manifests were not rewritten to hide the discrepancy.

The W2L108 checkpoint was retrieved from easteregg2 and preserved byte-for-byte locally at:

`D:\Users\Ben Bolen\NeuroCampaign_20261031\subnetworks\W2L108_verified`

Its editable project is `editable/W2L_VirtualWalker_Newer.aproj`, with both native exports, dependencies, six saved-page PDF/SVG/PNG sets, and original verification ledgers. All 58 manifest entries pass. Fresh static checks find 108 native neurons, 216 chemical connections, 56 adapters, no duplicate project object-definition IDs, and no issues across all 24 independently parsed CDATA drawings. InLinks/OutLinks IDs are references and are excluded from duplicate-definition counts. Native tables preserve every serialized field and resolve chemical endpoints and named type payloads. Full fresh semantic GUI/native parameter parity remains pending; the earlier paired-run evidence is retained as historical evidence only. No native load or new run was attempted.

Historical LH style source recovered directly from git: revision `23b0e785`, `Neuromechanical_Models/Biped_2xCPG_wSubs/Biped_2xCPG_wSubs.aproj`. The recovered bytes and hash are in W2L108. No source original was replaced. Existing historical type-comparison evidence is in the retrieved checkpoint; it has not yet been independently repeated for every current family.

## Findings to retain during design

- Keep the W2L hybrid's four ankle PF target swaps separate from its editable source. `runtime_context.json` preserves `ankle_runtime_target_edits`, recruitment changes, central-output changes, physical membranes, clocks, and readout contract. A target-swap record's original conductance is not necessarily the final scaled edge conductance; use the effective edge row for the executed value.
- The three graded production candidates retain defined zero-conductance pairs separately from nonzero chemical edges. Do not delete these or build additional zero placeholders during this audit. Whether a source gain prevents construction versus zeroes an existing compiled pair requires the actual branch/configuration review.
- Native W2L uses named type `SynAmp` as effective conductance; per-connection `G` is preserved but must not be treated as the effective tuning parameter. Native adapters retain their actual polynomial and channel fields. Physiological labels alone do not establish what a receptor measures.
- Source numeric capacitance on the seconds clock is not directly nF: `Cm_nF=1000*C_numeric` when conductance is uS and current is nA. Thus C=.03 corresponds to 30 nF. Do not silently rescale model parameters while correcting labels.
- The retrieved `CONTROL_DESIGN_LITERATURE_REVIEW.md` supersedes the earlier unresolved Micklem title/spelling evidence in the older literature review. The actual supplied Micklem PDF hash is `5c7ef65fc06440e0b8b87aad8506fc52fb924b8f1c8caa88836ea5c2cc3d5a87`. PDF page 3 was extracted and visually inspected: Figure 2 places graded RG/PF/motor circuits before a separate spiking output layer; Figure 3 uses physiological neuron, inhibitory-circle, and feedback conventions. This is a drawing reference, not an equivalent graph for our central-spiking Gait implementation.
- Syn6's six NMF factors remain distinct from Di Russo's five raised-cosine primitives and Stark's Ward.D2 clusters. Existing spatial synergy-weight maps are not muscle-by-stride activation heatmaps. No NMF basis was regenerated.
- Future learning requires tied absolute bilateral weights for verified homologues; unmatched/source-asymmetric paths remain frozen unless separately approved. This review has not changed parameterization or weights.

## Required condition selection and remaining work

Ben supplied authoritative anchors after the initial intake: Rat134 Curated134_ReceptorGeometry (refinement12 trial3 remains a diagnostic comparator); Gait_dynamic_readout_v4 TimingZeroSFA; CAD262 support085 mixed pressure pair via its reproduction manifest; and the fixed-h direct ODE Deng14 reference. These replace the earlier condition-selection questions. Named source evidence was retrieved without executing it. Rat134's supplied standalone hash,134/274census, and21saved drawing pages pass static checks. The Gait frozen configuration agrees byte-for-byte with its snapshot and all five listed controller dependencies match. All32entries in its file manifest now match retrieved evidence.35non-gait source-evidence entries match the retrieved sources; several EB475WS4 live files differ, so the local live source is not substituted for the saved evidence. Details are in NAMED_CONDITION_VALIDATION.json and NAMED_SOURCE_PARITY.csv. These checks support accurate drawings; they are not numerical circuit tests.

| Requested family | Intake status / remaining evidence |
| --- | --- |
| s3k, w2lvar, syn6 | Saved standing-condition ledgers extracted and source hashes checked; individual-pathway and scaling review pending. |
| biw2lvar-Hunt-PT | Need designated consumed configuration/runtime graph; no September template substituted. |
| Functional four-cell RF/HAM | Local Session1 manifest and editor export located. README explicitly says conductances combine two expression views and exclude disabled edges, so this is not a complete current runtime reconstruction. Need consumed source condition. |
| W2L108 graded | Verified native checkpoint retrieved; static checks completed, detailed physical-path and literature review pending. |
| W2L704/80 hybrid | Saved effective graph and runtime overlays extracted; seven-step spiking audit and expanded physiological drawings pending. |
| rat134 | Curated receptor-geometry anchor selected;134/274native node/edge/adapter tables extracted;21drawing pages pass. |
| CAD262 | Support085 mixed pair selected; six tied recruitment edges and physical pressure-conversion path drawn separately from the52hybrid routes. |
| Current Gait central spiking | TimingZeroSFA selected; frozen snapshot retained unchanged;256soma identities plus778graded identities tabulated; central interface drawing distinguishes proxy-voltage chemical connections from spike events. |
| SNS-Simscape graded/spiking reflex | Local builders/snapshots located; remote `spiking/reflex_simulink` exists. Need exact demonstration selection and static block/parameter audit. |
| Reduced Deng14 | Fixed-h ODE selected;14nodes/18edges extracted by static syntax inspection and drawn; no model import or integration. |

For every selected family, remaining deliverables are full physical sensor/conversion/receptor-state pathways, source-author versus engineering-rule comparisons with exact endpoints and publication/table references, operating-voltage sign analysis, seven-step GLIF formula/field checks or the correct graded analogues, overview and expanded editable physiological subnetworks, native/export parameter parity, and print-size inspection of all vector figures. Only one retrieved II preview and the Micklem source page were visually inspected in this intake; no claim of complete figure QA is made.

No new circuit family was invented. Editable JSON copies preserve existing exports; the W2L native checkpoint preserves existing saved pages. New connection drawings are independent SVG/DOT files, with their source edges in companion tables. They do not modify native projects, controllers, or browser drawings. Twenty initial overview/detail pages and four more focused W2L108 RG/PF feedback/key pages are saved. These are a partial drawing set, not completion of every requested family. Native loading and simulation checks remain unperformed; training remains gated.

## Reproduction

Run `audit_saved_graphs.py` with EB475WS4's myo Python for the four saved JSON ledgers. Run `audit_native_checkpoint.py` after retrieving the W2L108 checkpoint to the documented intake location under root `tmp/functional_subnetwork_remote`. Both scripts perform file parsing and static graph checks only. Neither imports a controller nor starts a simulator. Generated hashes and tables accompany this report. Remote source originals, unsaved user applications, and browser localStorage remain untouched.
