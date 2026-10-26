# Application rule — targeted patches, not chapter replacement

The live Overleaf document is the baseline. Do not upload/paste Revision E chapter files wholesale: they lack current tracked deletions, changes, and comment anchors. Use `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/targeted_reconciliation/OVERLEAF_TARGETED_PATCHES.md`. It contains exact-text replacements, comment/date evidence, and the current table/deletion decisions. Ben’s deleted sentence must stay deleted. Supply modeling evidence separately from prose patches; do not launch parallel chapter rewrites. No Overleaf changes or inter-agent messages were made by this reconciliation.

---

# Latest correction — original Figure 6.1 and conversion-repair appendix

Current output: `Revision_E_corrected.pdf` (270 pages). The earlier Revision_E_local.pdf was open and locked and remains the earlier version. Figure 6.1 now has (A)/(B) labels and uses Ben’s exact supplied original Connor–Ben controller as panel B. The main text distinguishes the original shared RG from Ben’s later per-side RG extension with RG/PF and contact feedback. The heading is “Block Programming Neural Pathways.” Main-text conversion details now link to a consolidated Appendix D.1 repair account, complementing the existing Appendix A force comparisons. See `review/figure61_correction/CORRECTIONS.md` for source lines and build evidence. Earlier pass reports below are preserved as history.

---

# Latest pass — submission clarity

The current PDF and chapter sources include a submission-focused prose pass across Introduction, Background, Database, SNS, Discussion, and Conclusion. Read `review/clarity_pass/CLARITY_PASS.md` and its line-based comment ledger for the newest editorial changes. The scientific evidence history below remains preserved. No Overleaf edits or new model runs were made.

---

# Revision E — updated evidence review

This version incorporates Ben's follow-up answers and a read-only inspection of easteregg2's newer source and outputs. It supersedes the first E draft's statement that the corrected replay was unavailable. The first E text, report, handoff, and matrix remain preserved in `review/60-sns_before_evidence_followup.tex`, `review/before_followup_REPORT.md`, `review/before_followup_ZCODE_HANDOFF.md`, and `review/self_scoring_matrix_before_followup.json`.

## What the chapter now establishes

The argument remains: we needed to separate neural behavior, translation fidelity, actuator/mechanical coupling, and supported locomotion before attributing observed motion to a mechanism. The confusing five-task table stays replaced by prose explaining the question, reason, tool, result, and limit.

- **Rybak:** 13/13 tests of selected published neural behaviors passed. PF perturbation preserves timing (0.001-period residual); RG perturbation resets it (0.416). Source inspection shows this is an activity-based reduction with calibrated global constants, changed pulse scaling, and integration settings. It reproduces selected qualitative neural results, not an exact spiking model or a walking experiment.
- **Deng:** our Simulink-versus-our-NumPy comparison verifies a translation chain and sustained bursting. The NumPy reference is also our translation. The original AnimatLab project was mined but not run for that comparison; the published perturbation triad remains untested. Period matching is not voltage-waveform identity.
- **Li:** matched imposed-contact and contact-removal responses check the contact-driven design principle. The 58/98 Simulink transcription substitutes graded dynamics. Neither result is a figure-by-figure quantitative reproduction of the original closed-loop experiment.
- **Di Russo:** implemented partial attempt; stage-2's 120 trials did not achieve stable gait. The absence of a stable result is not explained away as insufficient budget, and no larger search was launched here.
- **BPA:** the saved W2L gate demonstrates 12 Festo-model actuators operating with the controller: air-stepping and retained-support standing, peak 372 kPa under a 620 kPa ceiling. No matched Hill/BPA A/B exists. The bridge report tests the converted **human** model, not a BPA force law; the bridge's 549.9 N/87.8° muscle response must not be relabeled BPA validation.
- **Human feedback wean:** the added hand-drive scale reaches zero with fixed re-tuned feedback, but descending DRIVE remains active and the model stays supported and in place. It is not feedback intact versus fully disabled, and not zero total neural drive.

## Corrected W2L evidence and its limits

The initial search really used only six of 61 stored fields. The 55 inert fields were not optimized conductances. The old 10/10 mass-change and nominal-support ladder records belong to that old effective controller. Replaying the stored candidate with all fields applied changes the operating point; it does not retroactively optimize those fields.

I copied the newer uncommitted remote runner and raw logs into `review/easteregg2_followup/`, without changing remote files, and independently recomputed heel events using the panel's 20 N threshold and 0.45 s event-merging rule. Full 16 s records give:

| Condition | L/R contact cycles | Interpretation |
|---|---:|---|
| Control | 17/16 | Baseline of corrected replay |
| Crossed heel edge removed | 16/16 | Cycling retained |
| Ia chain removed | 13/13 | Cycling retained; not unchanged gait |
| Renshaw removed | 16/14 | Cycling retained |
| Heel cascade removed | 17/15 | Cycling retained |
| Commissural removed | 13/12 | Earlier necessity claim withdrawn |
| Ib feedback removed | 0/0 | No heel-contact cycles in this configuration |
| Mild / strong RG stimulus | 18/17; 19/17 | Strong stimulus reaches 62.8° tilt |
| Forward / backward push | 18/18; 18/16 | Tilt reaches 49.6° / 47.7° |

The stronger stimulus and both push records cross the runner's **45° fall threshold**, despite continued contact cycling. In panel F, “fell” applies to the deletion conditions only. It does not certify the push conditions. The corrected figure has been restored with a caption stating that limit. Do not summarize this as robust gait recovery.

The saved `changed_body_TRUEWINNER_w25.json` records **1/10 survival**, with **9/10 exceeding 45° tilt**. Contact-cycle counts range 32–38 and minimum pelvis height is **0.8981 m**, not “never below 0.90 m.” The currently retained `changed_body_test_w25.json` is the earlier 10/10 record. These files represent different effective configurations and must not be interchanged. The chapter calls the corrected result persistent cycling with inadequate posture robustness.

The remote runner now forwards gains-named fields already present in KV to the builder. The evaluator still filters out fields absent from KV, however. Direct panel runs populate KV themselves; every optimization and fresh-verification entry point still needs explicit registration and an effective-weight signature check. The mass replay's exact initialization/command must be preserved before claiming it can be reproduced from the saved driver alone. No renewed per-target optimization has been accepted.

## Advisor comparison and evidence questions

Ben confirms the Hunt/Fletcher libraries were never supplied and does not want to request them. That question is closed as unavailable evidence. The chapter avoids asserting a universally unmet graphical-library need. No message to the advisor or ZCode was sent. No new modeling or optimization was queued.

Shevtsova's reduced-model checks remain a recorded DESIGN.md result, without a surviving standalone log in the supplied inventory; they are not promoted to new independently reproduced evidence in this chapter. The initial broad “no reproductions exist” answer is superseded by the later Rybak inventory, with the source's reduction caveats retained.

## Source and provenance

Initial E baseline: `a988508a`; desktop after this follow-up fetch: `653523e3`. Remote easteregg2 HEAD: `f1514e9b`, with uncommitted runner/results changes. Snapshot paths and hashes: `review/easteregg2_followup/snapshot_manifest.json`. Independent raw-log audit: `review/easteregg2_followup/independent_audit.json`. User clarification precedence: `review/followup_answers_and_precedence.json`.

Additional local evidence: `spinal/replication/verify_rybak2006.py`, its `.log`, and `rybak2006_standalone.py`; `spinal/campaigns/20260930/bpa_run/bpa_gate.log`; `spinal/reports_20261016/wean_human/WEAN_CAMPAIGN_20261016.md`; `Code/Matlab/SNS_Simscape/mujoco_bridge/BRIDGE_REPORT.md`. The stage-2 Di Russo log remains the record of its failed stability search.

## Per-directive compliance

Section, exact line, changed text, rationale, remaining limit, and provisional scores are in `review/self_scoring_matrix.json`. Twelve entries now cover the earlier review plus Rybak and human-wean interpretation. All remain **draft score 2 / acceptance 0**, not a self-awarded clearance of repeated clarity complaints.

- **E01**, Knowledge Gaps and the New Tools, `chapters/60-sns.tex:27`: Replaced the five-task inventory and confusing table with question, reason, tool, result, and limit in connected prose. Open: Prior repeated clarity failures remain unaccepted. Ben/advisor comprehension review required.
- **E02**, The AnimatLab Bipedal Walker, `chapters/60-sns.tex:8`: Included four-neuron shared RG; separate hip and knee/ankle PF groups per leg; We added MN feedback. Open: Live reply captured; clarify who performed later tuning if needed.
- **E03**, Model Conversion and Actuator Coupling, `chapters/60-sns.tex:43`: Separated human reference from robot routing; split repairs into readable paragraphs; retained 92 entries / six inert. Open: Single-actuator checks do not prove full source-model equivalence.
- **E04**, Matched Neural Implementations, `chapters/60-sns.tex:111`: Added pulse, tutorial, short 410-neuron, and adapted Li tests with their actual comparison scope. Open: 63 mV Deng waveform discrepancy; Li neuron substitution; 0.3 s network horizon.
- **E05**, Closed-Loop Model Coupling, `chapters/60-sns.tex:123`: Added eight-BPA air-step result and 0.31 s mechanical / 1.938 s RG mismatch. Open: No entrainment or imported-body control established.
- **E06**, Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:133`: Distinguished old six-field optimization from fully applied candidate; 9/10 mass cases fail tilt gate. Open: Persistent contact cycling is not posture robustness; candidate has not been optimized after mapping correction.
- **E07**, Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:127`: Preserved original 55/61 loading failure as history; inspected newer remote runner and saved outputs. Open: Evaluator still filters non-KV keys; verify fresh-process campaign/verification drivers before redo.
- **E08**, Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:129`: Independently recounted remote NPZ logs; restored corrected figure with 17/16 control, 13/12 commissural, and 0/0 Ib. Open: Claims restricted to this operating point; push and strong-stimulus posture failures are explicit.
- **E09**, Perturbation Robustness and Pathway Deletions, `chapters/60-sns.tex:144`: Reported the completed 120-trial negative stability search as an adaptation. Open: Not a successful paper replication; no new simulation run here.
- **E10**, Capabilities and Remaining Limits, `chapters/60-sns.tex:156`: Libraries were never supplied; Ben does not want to request them. Avoided unsupported novelty claim. Open: Unavailable comparison is a documented limitation, not an action awaiting Ben.
- **E11**, Verification of the Implemented Interfaces, `chapters/60-sns.tex:94`: Selected Rybak published neural behaviors reproduced with documented reductions Open: Scoped evidence; human/advisor acceptance pending.
- **E12**, Verification of the Implemented Interfaces, `chapters/60-sns.tex:148`: Human wean removes added hand drive, not all descending drive or all feedback Open: Scoped evidence; human/advisor acceptance pending.

## Build, figures, and conservation

Fresh MiKTeX pdflatex → bibtex → pdflatex → pdflatex build: **265 pages**, zero errors, zero unresolved references/citations, zero missing figures. **26 document-wide overfull boxes** remain. Rendered SNS pages were checked for clipping and caption readability. `Revision_E_local.pdf` is the output; `build_20261026_040330` is the build directory.

The old defective perturbation figure remains preserved but unused. Its replacement is `figs/60-sns/w2l_perturbation_fully_applied.pdf`; the original figure label is restored. The confusing capability table stays excluded, with prior source in `review/replaced_capability_table.tex`. No active incoming reference to its removed label remains. Number changes and exact source diffs are in `review/conservation.json`, `review/changes.diff`, and `review/followup_changes.diff`. Other E chapter files and bibliography were not edited in this pass. New v4 demo files outside E remain untouched.

Overleaf was not edited. No commit, push, remote modification, or new simulation run was made. Local records and chapter were revised; raw outputs were inspected and recounted rather than newly simulated. `ZCODE_HANDOFF.md` now separates completed remote repairs from remaining evidence and entry-point checks.
