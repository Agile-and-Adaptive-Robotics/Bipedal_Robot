# Revision E — SNS chapter review

Local revision of Chapter 6 against Revision D at fetched commit `a988508a638d43181c70be9deca31c2b57079609`. Branch `KneeTestSetup_BenBo_stw` and its upstream matched after fetch. Recent work reviewed includes `36f60475` (SNS validation), `af290071` (Simscape figures/demos), `07ca92f7` (dissertation figures), and `a988508a` (SNS update). No merge, commit, push, or Overleaf edit was made.

## Main finding

The chapter needed both a clearer argument and a stricter account of what was actually evaluated. New modeling is substantive, but a populated result file or a successful build does not establish that all intended parameters or pathways were used. The repeated clarity complaints remain open for Ben and the advisor to judge.

**Critical:** at the reviewed source version, the W2L evaluator discards **55 of 61** winner fields. Only `ch_ank`, `ch_contra`, `comm`, `cref`, `ctn_amp`, and `ibnA` survive the accepted-key filter. The runner also does not pass the individual fields to the builder. See `review/parameter_plumbing_audit.json`. This invalidates attributing these outputs to a fully applied per-target conductance search. It does not erase every recorded mass-change result; it requires identifying and replaying the effective controller. A newer uncommitted easteregg2 source could differ.

The deletion figure also duplicates left-heel counts into the right-side column. Its `kneeL_deg` field selects hip column 0 and converts already-degree values again. Large stored values are a reporting bug, not evidence of physical knee rotations. The displayed knee traces use the separate correct knee column, so this finding does not itself invalidate those traces. Several deletion switches also fail to reach the builder. Unchanged runs cannot demonstrate biological redundancy when the intended deletion was not applied.

## Changes and retained evidence

- Replaced the five-task list/table with connected explanations of the question, why it matters, the tool, its test, and its limitation. The practical gaps are scoped to this workflow, without claiming that no comparable library exists.
- Incorporated Ben's latest 26 October reply: one shared four-neuron RG, separate hip and knee/ankle PF groups for each leg, and added afferent feedback at the MN level. Preserved the collaborative “We.” Latest replies and tracked edits are recorded in `review/overleaf_latest_sns_replies.json`; original conversion replies remain in the earlier review snapshot.
- Split the robot-first conversion and concrete repairs into readable paragraphs. Retained six inert actuator entries within the 92-entry human-model conversion. Removed broad claims that every source muscle or routing had been validated.
- Added the Deng-derived brief-pulse result, five-neuron and 410-neuron checks, adapted Li circuit, and BPA air-step test. Preserved the distinctions between final-state agreement, waveform agreement, imposed inputs, and true coupling. The Deng log's approximately 63 mV waveform discrepancy is explicit; period matching alone is insufficient.
- Preserved the ten-condition mass-change record, with its criterion and displacement range, as recorded supported progression in W2L. The support setting is nominal and the harness includes feedback and planar constraints. The support ladder permits re-tuning. The per-synapse records give failed basin results of **9/12 at w25 and 7/12 at w10**; the separate earlier class-search **10/12** is not substituted for these.
- Added the separate Di Russo adaptation's completed 120-trial stability failure. No successful original-paper walking replication is asserted.
- Replaced the stale blanket tonic-drive discrepancy with current short-horizon evidence and its limits. The 0.3 s result does not establish long-term oscillator or coupled-model equivalence.

## Explicit removals and number changes

The capability table (`tab:sns_capability`) and the defective perturbation figure (`fig:w2l_perturbation`) were removed from the compiled chapter. Their complete prior source remains in `review/replaced_capability_table.tex` and `review/excluded_perturbation_figure.tex`; original image assets remain copied in `figs`. Neither label has an active incoming reference in Revision E. This is a documented exclusion, not deletion of the historical evidence.

The unverified 620 kPa table claim, strong 40 N push/recovery interpretation, and “load-bearing pathway” conclusion were not carried forward as accepted results. Their original text and numbers remain in the preserved source/diff. New numerical claims are tied to the sources below and the evidence manifest. `review/conservation.json` records all numeric-token changes, including layout dimensions and duplicated table/caption values. All original citation keys were retained; one existing Di Russo citation was added. Bibliography bytes and all other chapter sources are unchanged from D.

## Self-scoring and compliance

Scores are deliberately provisional: **2 = revised local draft with evidence; acceptance 0 = not accepted by Ben/advisor**. Earlier repeated-failure scores are preserved in their historical matrix. A build does not promote a clarity score. The precise section, source line, changed text, rationale, limitation, and score for each concern are in `review/self_scoring_matrix.json`.

- **E01 — Repeated unclear gaps and contribution**, § Knowledge Gaps and the New Tools, `chapters/60-sns.tex:27`: Replaced the five-task inventory and confusing table with question, reason, tool, result, and limit in connected prose. **Open:** Prior repeated clarity failures remain unaccepted. Ben/advisor comprehension review required.
- **E02 — Ben reply: per-leg PF organization and MN feedback**, § The AnimatLab Bipedal Walker, `chapters/60-sns.tex:8`: Included four-neuron shared RG; separate hip and knee/ankle PF groups per leg; We added MN feedback. **Open:** Live reply captured; clarify who performed later tuning if needed.
- **E03 — Actual conversion failures and robot-first attempt**, § Model Conversion and Actuator Coupling, `chapters/60-sns.tex:43`: Separated human reference from robot routing; split repairs into readable paragraphs; retained 92 entries / six inert. **Open:** Single-actuator checks do not prove full source-model equivalence.
- **E04 — New graphical implementation evidence**, § Matched Neural Implementations, `chapters/60-sns.tex:105`: Added pulse, tutorial, short 410-neuron, and adapted Li tests with their actual comparison scope. **Open:** 63 mV Deng waveform discrepancy; Li neuron substitution; 0.3 s network horizon.
- **E05 — Coupled mechanics and air stepping**, § Closed-Loop Model Coupling, `chapters/60-sns.tex:115`: Added eight-BPA air-step result and 0.31 s mechanical / 1.938 s RG mismatch. **Open:** No entrainment or imported-body control established.
- **E06 — W2L mass changes, support, and robustness**, § Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:119`: Preserved recorded ten-condition result; defined criterion; separated nominal support, actual force, and re-tuning. Included failed basin outcomes 9/12 and 7/12. **Open:** No harness-free walking; source/result provenance requires reconciliation.
- **E07 — Per-target optimization not applied**, § Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:125`: Static audit found 55/61 winner keys discarded; removed acceptance of per-target optimization claims. **Open:** Priority ZCode repair and effective-weight replay required.
- **E08 — Deletion and perturbation figure overclaims**, § Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:127`: Excluded defective figure; qualified recorded left counts; identified unapplied switches and missing recovery measure. **Open:** Need corrected bilateral logs and actual edge-removal checks.
- **E09 — New Di Russo modeling**, § Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:129`: Reported the completed 120-trial negative stability search as an adaptation. **Open:** Not a successful paper replication; no new simulation run here.
- **E10 — Comparison requested by advisor**, § Capabilities and Remaining Limits, `chapters/60-sns.tex:137`: Avoided claiming a universally unmet library need. Asked Ben for Hunt/Fletcher library identity. **Open:** Comparison remains open; do not score as satisfied.

## Evidence pointers

- `Code/MuJoCo_SNS/spinal/reports_20261016/w2l_wean/_w2l_wean_eval.py:128`
- `Code/MuJoCo_SNS/spinal/w2l_mujoco/test_w2l_ground.py:291`
- `Code/MuJoCo_SNS/spinal/reports_20261025_persyn/_poster_panels.py:93`
- `Code/MuJoCo_SNS/spinal/reports_20261025_persyn/_poster_panels.py:95`
- `Code/MuJoCo_SNS/spinal/w2l_mujoco/test_w2l_ground.py:718`
- `Code/MuJoCo_SNS/spinal/w2l_mujoco/test_w2l_ground.py:613`
- `Code/Matlab/SNS_Simscape/sns_verify_from_json.m:51`
- `Code/Matlab/SNS_Simscape/logs/sns_run_deng2.log:5`

Full source paths and SHA-256 hashes are in `review/evidence_manifest.json`. The Di Russo run log records no passing stage-2 trial. The older connectome audit is used as historical analysis, not as a current certification after subsequent changes.

## Build and visual review

Fresh MiKTeX pdflatex → bibtex → pdflatex → pdflatex build: **262 pages**, zero errors, zero unresolved references/citations, and zero missing figures. There are **26 overfull boxes** across the full dissertation; these remain a document-wide layout limitation. The rendered Chapter 6 pages (printed 119 onward) were checked for clipping and figure placement. See `review/contactsheet.png` and `review/build_result.json`. No new model evaluation or independent reproduction of simulation outcomes was performed.

Output: `Revision_E_local.pdf`. Reviewable source: `chapters/60-sns.tex`. Build directory: `build_20261026_030543`. The chapter preserves existing figures except the explicitly excluded perturbation figure. Other chapters remain the Revision D baseline, including any inherited issues outside this SNS pass.

New uncommitted v4 demo figures and their generator appeared during this session. They were inspected selectively and left untouched, outside the compiled E chapter. The knee-reflex v4 image reports a maximum 3.6° difference from its NumPy twin; it should not be called exact agreement. Approval/interpretation of these evolving figures remains pending.

## Evidence questions sent to Ben

1. Are there newer original-paper protocol comparisons, especially Deng perturbations and Li closed-loop contact feedback?
2. Are there matched feedback-on/off or BPA-model comparisons beyond the committed reports?
3. Is there a corrected perturbation figure or raw bilateral log on easteregg2?
4. What are the exact Hunt/Fletcher libraries that the advisor requested for comparison?
5. Does easteregg2 have uncommitted evaluator/runner changes that apply the individual conductances? Recover exact source and output hashes together.

These questions remain pending. E scopes the claims to available evidence and does not invent missing results. `ZCODE_HANDOFF.md` specifies repairs and the smallest useful checks before additional optimization.
