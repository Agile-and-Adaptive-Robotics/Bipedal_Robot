> Superseded by [the current committee revision](../committee_revision_20261027/REVIEW_REPORT.md). Working source is now ProofFinal; do not follow the old source-directory or outstanding-figure instructions below.

# Dissertation committee review copy

The review candidate is **Bolen_Dissertation_Committee_Review_20261027.pdf**, with matching source in **Bolen_Dissertation_Committee_Source_20261027.zip**. It is 299 pages and builds without errors, unresolved citations, unresolved references, duplicate labels, or rerun warnings. One 0.67-point overfull text line remains in Appendix C; it is visually negligible. The deadline used here is **October 28, 2026, 9 a.m. Pacific, for advisor and committee review**, as Ben clarified. This is not an eTD submission package or a claim that every scientific result has been independently reproduced.

**One mechanics audit remains before calling the work finished:** the current flexor evaluator applies a point transformation to force vectors. Its effect on the adopted fit and predictions must be measured. The exact evidence and a bounded easteregg2 task are in `EASTEREGG2_HANDOFF.md`. No fitted coefficients, stored results, or source simulation code were changed in this pass.

## Which source was used

The downloaded `Dissertation/Bolen_Dissertation.zip` already contains the advisor-facing introduction, SNS, Discussion, conclusion, and Appendices D and E, as well as recent torque material. Its SNS wording was checked against the open Overleaf editor. This is newer than the split-tree account in the project digest. The digest's statement that every available PDF must lack one set of revisions is therefore not a reliable description of this download.

The candidate was built from that complete ZIP. It was not assembled by replacing newer chapters with October 26 snapshots. Original ZIP SHA-256: `12ae70f667dfe6fdb524b796447480151e1ea91adfe3da678c04dd43821982af`. The original archive, original downloaded PDF, local ProofFinal, and live Overleaf were left intact. Any Overleaf edits after the download still need reconciliation before uploading candidate files to the live project. No comments were marked resolved in Overleaf.

## Changes and directive checks

Line references below refer to `source/chapters/` in this review folder. The complete before/after record is `changes_from_download.diff`; it includes every displaced passage rather than silently discarding it. `review_verification.json` records file hashes and build checks.

| Request or finding | Disposition and evidence |
|---|---|
| Leave Table 4.3 alone | Honored. `tab:FlxPinPreGof` is exactly unchanged. All existing Chapter 4 table environments compare exactly with the download. Hunt comparison corrections P06/P07 were excluded. |
| Preserve newer advisor work | Used the complete latest download, retaining its topic-chapter structure and Appendices D/E. Introduction and Background are unchanged. |
| Explain which torque fit is adopted | `30-torque.tex:458` now leads with the adopted two-bracket fit, its training/validation split, and its existing adopted values. The displaced initial single-bracket values are preserved in Appendix C's supporting fit section. No existing table values changed. |
| Physical torque figures and claims | Applied the previously prepared P08/P09 figures. The truncated flexor caption is complete. `30-torque.tex:698` onward distinguishes nine extensor and fifteen flexor measurements, selected single-muscle targets, prediction at measured pressure, and nominal margins. The eight-above/one-below prediction result is stated correctly. Abstract and conclusion are aligned with the sampled evidence. |
| Do not confuse prediction with demonstrated performance | The physical-test section explicitly limits conclusions to sampled angles and pressures and states that measurement-uncertainty bounds are absent. The Discussion labels the full-envelope claim as a design prediction (`70-overall.tex:12`). |
| Resolve Discussion contradiction | `70-overall.tex:49` now describes the adopted separate staggered flexor paths, replacing the stale statement that both actuators share one path. The whole-model geometry claim is also narrowed to the actual checks. |
| Appendix A mechanics and biomechanics | Corrected pennation in fiber–tendon equilibrium (`92-AppendixA.tex:9`), the tensile moment-arm sign convention, MuJoCo fiber-length/velocity normalization and pulling-force sign (`:52`), and missing origin-bracket compliance in the series expression (`:173`). |
| Appendix A force-comparison interpretation | `92-AppendixA.tex:62` explains that the script uses static-optimization reference forces, matched activations, sampled positions, and no assigned trajectory velocities. The saved statistics are retained. Claims that this comparison demonstrates compliant-tendon dynamics or physical delay were withdrawn. Conventional R-squared is distinguished from the MATLAB sample-scaled FVU (`:119`). |
| Appendix B equations and appearance | The synaptic-current sign is consistent with inward-positive membrane current (`93-AppendixB.tex:51`). Wide parameter/calibration/conversion tables wrap within the page. Network diagrams are enlarged on rotated-content pages. Broad assertions about every toolbox implementation and documentation being unable to drift are narrowed. |
| Appendix C coordinate equation | `94-AppendixC.tex:10` includes the translation needed to express a point from the femur frame in the tibia frame. This matches the point transformation used to construct the route direction. This correction is distinct from the unresolved misuse of that helper for force vectors. |
| Database readability and references | The repeated SADb explorer and supporting tables now appear once, in existing Appendix E, with Chapter 5 references. This removes three multiply-defined labels. Novelty and automatic physiology-to-controller claims are narrowed; `60-sns.tex:80` explains what database evidence supplies and what still needs parameter selection and testing. |
| Bibliography and navigation | Existing `thesis.bib` bytes were preserved exactly, with two documentation entries appended. Long URLs and appendix paths wrap. Bibliography and Appendices bookmarks now point to the correct pages. Title-page month is October for this review copy. |
| Chapter 3 publication content | Only the malformed S'arosi author-name accent was corrected. Scientific content and numbers were preserved. |

## The later column

| Item | Review-copy decision |
|---|---|
| Capability v3 and vastus-medialis actuator-count study | Already present in the downloaded Chapter 4, including the unsuccessful five-angle envelope result. Retained; no duplicate section added. |
| Independently replayed W2L result | Included in Appendix D, Section D.5.1 (`95-AppendixD.tex:233`). The independent replays, perturbation checks, analysis window, and residual harness force are explicit. It remains supported stepping. |
| Exact accelerated neural backend | Included in the same subsection, limited to the isolated neural-step speed improvement and trajectory parity. No whole-simulation speed claim. |
| CAD-derived AH001 robot | Included in Appendix D, Section D.5.2 (`95-AppendixD.tex:239`), as a supported stance/integration test. The 94.68% body-weight support and reference-tension rather than calibrated-BPA interface prevent it from being mistaken for independent balance or completed actuator validation. |
| Deng original/translated comparison | Included at `95-AppendixD.tex:54`. Different periods/parameterizations are distinguished from qualitative agreement in perturbation response. The subsequent correction of the false suppression interpretation is incorporated. |
| Direct Li MuJoCo port | Not promoted to a walking result. The report has unilateral contact and no left-foot engagement. It is an unfinished replication task. The existing original-simulator Li result remains separately scoped. |
| Di Russo adaptation | Not added as a positive locomotion result: the duration/speed search selected a somersault, which the report rejects. Original-controller execution remains unavailable locally. It is not needed to complete the dissertation's tools contribution. |
| Further optimization, bilateral-RG restructuring, and support removal | Future work. They need new experimental results, not just prose inclusion. No last-minute tuning was launched. |
| Defense deck, Hyfydy, and unaccepted optional schematics | Deferred. They do not repair the committee copy's scientific argument. Existing figure candidates are not equivalent to approved new dissertation figures. |

The digest's unused-47.9-cm main-text sentence and “40.2 N m at 804 N” clause were not present in the latest downloaded locations described by the digest. No speculative replacement was made. The historical unused-test bookkeeping remains in the appendix.

## Evidence for added or corrected claims

- Physical extensor/flexor records: `Notes/advisor_revision_20261023/REVIEW_INBOX/easteregg2/20261026_1700_torque-physical-audit/REVIEW.md`, its `csv/Ext20mm_series_audit_20261026.csv`, and P01/P08/P09 patches and figure assets. The extensor target is vastus medialis, consistent with `Code/Matlab/Mesh_Optimization/buildKneeExtContext20mm.m`, target selection around lines 119–140.
- W2L replays and exact backend: `Code/MuJoCo_SNS/spinal/reports_20261027_codex/REPORT.md`, `verification.json`, and companion support/perturbation records. The 2.246 m displacement belongs to the 4.5–30 s analysis interval, not an unqualified whole-run average-speed calculation.
- CAD robot: `Code/MuJoCo_SNS/spinal/reports_20261027_robot_control/SESSION_1_UPDATE_20261027.md` and `Code/MuJoCo_SNS/robot_torso/controller_assessment.json`.
- Deng: `reports_20261026_dengoriginal/DENG_ORIGINAL_20261026_REPORT.md`, together with the superseding `reports_20261026_dengmujoco_20261026/TWIN_VERDICT_CORRECTION.md` under the spinal folder.
- Force-comparison scope: `Code/MuJoCo_SNS/spinal/_muscle_force_compare.py`, especially its second matched-activation pass and saved CSV. The script comments overstate dynamic interpretation; its executed data assignments are the evidence used here.
- Primary technical references: [MuJoCo 2.3.7 muscle documentation](https://mujoco.readthedocs.io/en/2.3.7/modeling.html#muscle-actuators), [MuJoCo 2.3.7 implementation](https://raw.githubusercontent.com/google-deepmind/mujoco/2.3.7/src/engine/engine_util_misc.c), [OpenSim static optimization documentation](https://opensimconfluence.atlassian.net/wiki/spaces/OpenSim/pages/53085189/Working+with+Static+Optimization), and the already-cited Thelen and Millard papers. The OpenSim documentation citation has no invented publication date.

## What remains

1. **Mechanics audit:** establish the numerical effect of the force-vector transformation and resolve the adopted-run bracket-frame lineage. The current code executes pitch-only frames at both brackets; its header and standing context describe a two-rotation origin frame. The appendix follows executable code, but that alone does not prove which version generated the adopted fit. Do not select a convention by guessing or change Table 4.3.
2. **Measurement uncertainty:** the nominal minimum margins are small. Calibration, angle, load-cell placement, and repeatability information are needed to quantify uncertainty. The committee copy now states the limitation; it does not invent error bars. The worksheet's alternative torque rows should not be substituted for the audited record without resolving their different reductions.
3. **Figure-label cleanup:** the supplied extensor graphic still uses the shorthand “vasti” and “Back-calculated” in its embedded legend; the caption states the precise vastus-medialis and pressure-specific-prediction meanings. Appendix A's inherited plot still says “compliant” in its legend, and its corrected caption explains the static-optimization source. Regenerate those legends from the original data without changing curves or statistics. Dense connectome and route figures remain vector-zoomable but retain small source labels; their full redesign is separate from the completed page-layout pass.
4. **Live-source acceptance:** review the candidate and reconcile any Overleaf changes made after the downloaded ZIP. The 137-entry advisor matrix and live comment acceptance state have not been independently closed in this pass. The document can be reviewed as a committee draft; calling every advisor concern resolved would overstate this audit.

## Build and visual checks

Built on EB475WS4 with the full MiKTeX executable paths in a fresh `source/build_final` directory, absolute BIBINPUTS/BSTINPUTS, and five successful passes including BibTeX. Final log has no unresolved items or label-rerun warning. All appendix pages were inspected in rendered contact sheets; changed equations, wide tables, torque figures, and enlarged diagrams were also sampled at full-page resolution. This is a layout and targeted scientific-consistency review, not an independent rerun of every experiment or a line-by-line certification of all 299 pages.

No commits, pushes, emails, Overleaf writes, or remote jobs were performed. Preserve the source ZIP, PDF, this report, handoff, diff, and verification JSON. Temporary `build_*` and QA render folders need not be committed; Ben controls any commit through GitHub Desktop.
