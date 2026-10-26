from pathlib import Path
import json,hashlib,re,difflib,shutil
from collections import Counter
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review';F=V/'easteregg2_followup'
s=(E/'chapters/60-sns.tex').read_text(encoding='utf-8');b=(V/'60-sns_revision_D.tex').read_text(encoding='utf-8')
def li(t):return next(i for i,l in enumerate(s.splitlines(),1) if t in l)
def nums(t):return Counter(re.findall(r'\d+(?:\.\d+)?',t))
labels=lambda t:set(re.findall(r'\\label\{([^}]*)\}',t))
(V/'changes.diff').write_text(''.join(difflib.unified_diff(b.splitlines(True),s.splitlines(True),fromfile='Revision_D/chapters/60-sns.tex',tofile='Revision_E/chapters/60-sns.tex')),encoding='utf-8')
(V/'followup_changes.diff').write_text(''.join(difflib.unified_diff((V/'60-sns_before_evidence_followup.tex').read_text(encoding='utf-8').splitlines(True),s.splitlines(True),fromfile='E_before_followup',tofile='E_updated')),encoding='utf-8')
audit={'baseline':'review/60-sns_revision_D.tex, copied at a988508a','removed_labels':sorted(labels(b)-labels(s)),'added_labels':sorted(labels(s)-labels(b)),'removed_numeric_tokens':list((nums(b)-nums(s)).elements()),'added_numeric_tokens':list((nums(s)-nums(b)).elements()),'bibliography_modified_by_this_pass':False,'other_chapters_modified_by_this_pass':False}
(V/'conservation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
assert 'fig:w2l_perturbation' in labels(s) and audit['removed_labels']==['tab:sns_capability']
result=json.loads((V/'build_result.json').read_text());B=Path(result['build']);assert result['errors']==0 and result['undefined_references']==0 and result['missing_figures']==0
shutil.copy2(B/'main.pdf',E/'Revision_E_local.pdf')
snapshot={'machine':'easteregg2','remote_head':'f1514e9b8cef18ab95feb906a91f5d87b201d716','remote_dirty':True,'desktop_head_after_fetch':'653523e3','acquisition':'Read-only SSH/SCP after hostname verification; no remote file edits or model runs','files':[{'path':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(F.iterdir()) if p.is_file() and p.name not in ('snapshot_manifest.json',)]}
(F/'snapshot_manifest.json').write_text(json.dumps(snapshot,indent=2),encoding='utf-8')
confirm={'source':'Ben supplied follow-up answers in this chat, 26 October 2026','precedence':['Initial duplicated-count results','Count-corrected but old effective-controller results (13/8 control)','Runner-forwarding corrected replay (17/16 control), governs current chapter'],'confirmed':['No direct Deng perturbation reproduction; original AnimatLab project mined, not run for comparison','Di Russo implemented and staged, but no stable-gait reproduction','Rybak has 13/13 selected neural-behavior checks; source explicitly uses activity-based reduction and calibrated stand-ins','Hunt/Fletcher libraries were not supplied; Ben does not want to request them','Original 61-field search affected only six fields; 55 were inert noise','New replay applies stored fields, but those values were not optimized under the corrected mapping'],'independently_checked':'Local Rybak script/log; bridge report; human wean report; BPA gate; copied remote source, JSONs, and bilateral raw logs','limitations':['No new simulation run here','Remote source is uncommitted; preserve hashes and actual command before redo','Mass TRUEWINNER JSON lacks preserved raw body-change logs/driver initialization in this snapshot; treat as saved results requiring replay provenance','Evaluator still filters keys absent from T.KV; runner fix alone does not make every campaign entry point safe']}
(V/'followup_answers_and_precedence.json').write_text(json.dumps(confirm,indent=2),encoding='utf-8')
report=f'''# Revision E — updated evidence review

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

'''
matrix=json.loads((V/'self_scoring_matrix.json').read_text())
for row in matrix:
 if row['id']=='E08':row['text']=s.splitlines()[li('After correcting the runner')-1];row['line']=li('After correcting the runner')
 report+=f"- **{row['id']}**, {row['section']}, `chapters/60-sns.tex:{row['line']}`: {row['change']} Open: {row['limit']}\n"
(V/'self_scoring_matrix.json').write_text(json.dumps(matrix,indent=2),encoding='utf-8')
report+=f'''
## Build, figures, and conservation

Fresh MiKTeX pdflatex → bibtex → pdflatex → pdflatex build: **{result['pages'][0]} pages**, zero errors, zero unresolved references/citations, zero missing figures. **{result['overfull_boxes']} document-wide overfull boxes** remain. Rendered SNS pages were checked for clipping and caption readability. `Revision_E_local.pdf` is the output; `{B.name}` is the build directory.

The old defective perturbation figure remains preserved but unused. Its replacement is `figs/60-sns/w2l_perturbation_fully_applied.pdf`; the original figure label is restored. The confusing capability table stays excluded, with prior source in `review/replaced_capability_table.tex`. No active incoming reference to its removed label remains. Number changes and exact source diffs are in `review/conservation.json`, `review/changes.diff`, and `review/followup_changes.diff`. Other E chapter files and bibliography were not edited in this pass. New v4 demo files outside E remain untouched.

Overleaf was not edited. No commit, push, remote modification, or new simulation run was made. Local records and chapter were revised; raw outputs were inspected and recounted rather than newly simulated. `ZCODE_HANDOFF.md` now separates completed remote repairs from remaining evidence and entry-point checks.
'''
(E/'REPORT.md').write_text(report,encoding='utf-8')
handoff='''# ZCode handoff — Revision E evidence follow-up

Read Revision_E/REPORT.md, review/followup_answers_and_precedence.json, the remote snapshot manifest, and review/self_scoring_matrix.json first. The newest evidence supersedes the earlier 13/8 control and commissural-collapse account. Keep those records as history, not current claims.

## Current verified scope

- Original persyn search used six effective fields; 55 were inert. It was not a 61-field optimization. The stored values applied by the repaired runner form a replayed candidate, not an optimized winner under corrected dynamics.
- Independently recounted remote raw logs: control 17/16; crossed heel 16/16; Ia 13/13; Renshaw 16/14; heel cascade 17/15; commissural 13/12; Ib 0/0. Withdraw commissural necessity and all earlier no-op/zero-sensitivity redundancy claims. Ib dependence is specific to this operating point.
- Strong stimulus tilt 62.8°; pushes 49.6° and 47.7°. These cross the runner's 45° fall threshold, regardless of continued contacts. Panel F's “fell” entries cover deletions only.
- changed_body_TRUEWINNER_w25.json: 1/10 survival, 9/10 tilt failures, minimum pelvis 0.8981 m, cycles 32–38. Do not interchange with the preserved earlier 10/10 file. Store the exact full-field initialization and command used to generate this replay.

## Remaining technical work before any new search

1. Preserve the uncommitted easteregg2 runner and outputs with hashes. Do not overwrite them from the desktop baseline. The copied snapshot is in E/review/easteregg2_followup.
2. The runner forwards KV entries naming builder gains. `_w2l_wean_eval.evaluate` still filters keys absent from T.KV. Make input registration explicit and consistent in campaign, verification, sensitivity, body-test, and panel entry points. Unknown search keys should fail loudly. A direct panel run that seeds all KV keys does not prove a fresh campaign subprocess does so.
3. For every searched or deleted field, prove the intended effective edge weights change, unrelated weights remain fixed, and the complete compiled configuration reloads in a fresh process. Save signatures and effective weights. Preserve default baseline gates and the s3k anchor.
4. Keep old six-field and corrected full-field results separate. Rerun optimization only after these entry-point checks, valid posture/contact/phase gates, and an agreed objective. Do not weaken tilt gates to turn persistent cycling into walking.
5. Use matched fixed-controller feedback-intact/disabled tests only after correct application is proven. They are absent on the walking BPA model; human added-drive weaning is not that ablation. Explicitly set push/stimulus defaults and record actual harness forces.

## Interpretation and prose

Rybak's 13/13 result reproduces selected qualitative neural behaviors in an activity-based reduction with calibrated constants and modified pulses/integration. Do not call it a full exact reproduction of the spiking paper. Deng/NumPy period matching verifies our own translation chain; original-model perturbation comparison remains open. Li's deafferentation latch is a design-claim check, not full published-experiment replication. Di Russo is an implemented partial attempt with failed stage-2 stability, not reproduced stable gait. Shevtsova has a DESIGN entry but no surviving standalone log in the supplied inventory.

The bridge test uses converted human-model muscles; do not call its 549.9 N response BPA force validation. Separate W2L BPA gate: 12 Festo actuators, air stepping plus retained-support stand, peak 372 kPa within 620 kPa; no identical-mechanics Hill/BPA A/B. Human hand-drive scale zero retains descending DRIVE and support, with in-place activity.

Ben was never given Hunt/Fletcher's libraries and does not want to request them. Treat that comparison as unavailable and avoid an unsupported novelty claim. Do not assign Ben that follow-up again. No advisor or inter-agent message is authorized by this handoff.

## Workflow

Check hostname first and use only that machine's AGENTS.md paths. Read dissertation, animatlab, mujoco, matlab, neuroscience, sns-toolbox, and opensim skills from local ZCode_Skills and use relevant configured plugins. Ben prefers GLM-5.3/Flash for code/model/testing work and ChatGPT Astra for final prose; heavy computation belongs on easteregg2. No new runs or optimization were launched or queued by this review. Ben controls initiation. No Overleaf edits, commits, or pushes.

For the next evidence packet, provide exact source/model hashes, effective parameters, protocol, raw bilateral logs, posture/phase/contact/support criteria, outputs, limitations, and failed cases. Every local dissertation edit ends with a clean build and section/line/text compliance matrix. Clarity remains unaccepted until Ben/advisor judge it; a successful build is not a clarity score.
'''
(E/'ZCODE_HANDOFF.md').write_text(handoff,encoding='utf-8')
header='''# Revision E follow-up — corrected remote evidence, 26 October 2026

The latest E supersedes the earlier “corrected evidence unavailable” status below. Read-only easteregg2 retrieval and independent raw-log recount confirm corrected control 17/16, commissural deletion 13/12, and Ib deletion 0/0. The runner forwarding repair changes the operating point; it does not retroactively make the old six-field search a per-synapse optimization. Nine of ten corrected mass-change cases fail tilt, and both push conditions exceed the 45° fall threshold. E now states those limits and restores the corrected figure.

Rybak's selected neural-behavior tests (13/13) are included with activity-based-reduction caveats. Deng and Di Russo experiment-level reproduction remains open. The human wean removes added hand drive, not all descending drive or all feedback. Ben does not want to request the never-supplied Hunt/Fletcher libraries; that is no longer an outstanding request to him.

Updated source, PDF, report, handoff, snapshot hashes, independent count audit, and self-scoring matrix are under `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/`. Earlier versions are preserved in its review folder. No Overleaf edits, remote edits, new model runs, commits, pushes, or external messages.

---

'''
p=R/'CHATGPT_REPORT.md';old=p.read_text(encoding='utf-8');
if not old.startswith('# Revision E follow-up'):p.write_text(header+old,encoding='utf-8')
p=R/'ZCODE_HANDOFF_20261023.md';old=p.read_text(encoding='utf-8')
if not old.startswith('# Latest follow-up — Revision E'):
 p.write_text('# Latest follow-up — Revision E\n\nRead the updated `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/ZCODE_HANDOFF.md` first. Remote forwarding repair and new logs are preserved; commissural necessity is withdrawn, Ib removal gives 0/0, and 9/10 corrected mass cases fail tilt. Campaign/fresh-process input registration still needs checking before any new search. Hunt/Fletcher libraries were not supplied; Ben does not want to request them. Earlier handoffs below are historical.\n\n---\n\n'+old,encoding='utf-8')
print(json.dumps({'pdf':str(E/'Revision_E_local.pdf'),'build':result,'source_rows':len(matrix),'snapshot_files':len(snapshot['files'])},indent=2))
