from pathlib import Path
import re,json,difflib,hashlib,shutil,subprocess
from collections import Counter
R=Path(r'D:\Github\Bipedal_Robot'); P=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023'
D=P/'Revision_D'; E=P/'Revision_E'; V=E/'review'; S=R/'Code/MuJoCo_SNS/spinal'; M=R/'Code/Matlab/SNS_Simscape'
def line(p,text):
 return next(i for i,l in enumerate(p.read_text(encoding='utf-8',errors='replace').splitlines(),1) if text in l)
def ref(p,text): return f'`{p.relative_to(R).as_posix()}:{line(p,text)}`'
a=(E/'chapters/60-sns.tex').read_text(encoding='utf-8'); b=(D/'chapters/60-sns.tex').read_text(encoding='utf-8')
(V/'changes.diff').write_text(''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='Revision_D/chapters/60-sns.tex',tofile='Revision_E/chapters/60-sns.tex')),encoding='utf-8')
num=lambda s:Counter(re.findall(r'\d+(?:\.\d+)?',s))
label=lambda s:set(re.findall(r'\\label\{([^}]*)\}',s))
audit={'removed_labels':sorted(label(b)-label(a)),'added_labels':sorted(label(a)-label(b)), 'removed_numeric_tokens':list((num(b)-num(a)).elements()),'added_numeric_tokens':list((num(a)-num(b)).elements()),'bib_bytes_identical':(D/'thesis.bib').read_bytes()==(E/'thesis.bib').read_bytes(),'other_chapters_identical':all((p.read_bytes()==(E/'chapters'/p.name).read_bytes()) for p in (D/'chapters').glob('*.tex') if p.name!='60-sns.tex'),'original_citation_keys_retained':set(re.findall(r'\\cite\w*\{([^}]*)\}',b))<=set(re.findall(r'\\cite\w*\{([^}]*)\}',a))}
(V/'conservation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
assert audit['bib_bytes_identical'] and audit['other_chapters_identical'] and audit['original_citation_keys_retained']
for label_ in audit['removed_labels']:
 assert not any(re.search(r'\\(?:ref|autoref|cref)\{'+re.escape(label_)+r'\}',p.read_text(encoding='utf-8',errors='replace')) for p in (E/'chapters').glob('*.tex'))
B=Path((V/'latest_build.txt').read_text().strip()); result=json.loads((V/'build_result.json').read_text())
assert str(B)==result['build'] and result['errors']==0 and result['undefined_references']==0 and result['missing_figures']==0
shutil.copy2(B/'main.pdf',E/'Revision_E_local.pdf')
sources=[S/'reports_20261025_persyn/_poster_panels.py', S/'reports_20261016/w2l_wean/_w2l_wean_eval.py',S/'w2l_mujoco/test_w2l_ground.py',S/'reports_20261025_persyn/changed_body_test_w25.json',S/'reports_20261025_persyn/bodytest_trials.csv',S/'reports_20261025_persyn/basin_w10.json',S/'reports_20261025_persyn/basin_w25.json',S/'dirusso_gait2392/stage2_run.log',M/'sns_verify_from_json.m',M/'reports_20261024/DISSERTATION_EVIDENCE_INVENTORY_20261025.md',M/'logs/sns_run_deng2.log']
manifest=[{'path':str(p.relative_to(R)).replace('\\','/'),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
(V/'evidence_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
live={'project':'https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb','read_date':'2026-10-26','scope':'SNS review entries, including historical threads displayed at EOF; no source edits or review actions','new_reply':{'advisor':'ajh26','advisor_time':'23 October, 10:49 am','comment':'for each leg?','author':'bbolen83','time':'26 October, 2:11 am','reply':'yes. one hip PF for each leg and one PF for Knee and Ankle, per leg, all controlled originally with a single four neuron RG layer. ','position':1364,'anchor_status':'live character offset; NOT a line number'},'new_tracked_edits':[{'author':'bbolen83','time':'26 October, 2:13 am','position':1371,'from':'I','to':'We'},{'author':'bbolen83','time':'26 October, 2:13 am','position':1398,'to':'at the MN level.','source_in_context':'We added afferent feedback at the MN level.'}],'historical_replies':'Verbatim conversion replies remain in Notes/advisor_review_20261023/ben_conversion_replies_20261023.json. Live 23 October 2:03 and 2:05 replies were re-opened; 2:07 reply preview and historical snapshot were read.','uncertain_anchors':'Several 21 October threads share live position 14917 (EOF). Do not treat these as located source passages. Historical matrix anchors remain authoritative only for their frozen source snapshot.'}
(V/'overleaf_latest_sns_replies.json').write_text(json.dumps(live,indent=2),encoding='utf-8')
rows=[
('E01','Repeated unclear gaps and contribution','Knowledge Gaps and the New Tools','The suspended walker raised','Replaced the five-task inventory and confusing table with question, reason, tool, result, and limit in connected prose.','Prior repeated clarity failures remain unaccepted. Ben/advisor comprehension review required.'),
('E02','Ben reply: per-leg PF organization and MN feedback','The AnimatLab Bipedal Walker','For a course project','Included four-neuron shared RG; separate hip and knee/ankle PF groups per leg; We added MN feedback.','Live reply captured; clarify who performed later tuning if needed.'),
('E03','Actual conversion failures and robot-first attempt','Model Conversion and Actuator Coupling','My first conversion','Separated human reference from robot routing; split repairs into readable paragraphs; retained 92 entries / six inert.','Single-actuator checks do not prove full source-model equivalence.'),
('E04','New graphical implementation evidence','Matched Neural Implementations','A separate Simulink','Added pulse, tutorial, short 410-neuron, and adapted Li tests with their actual comparison scope.','63 mV Deng waveform discrepancy; Li neuron substitution; 0.3 s network horizon.'),
('E05','Coupled mechanics and air stepping','Closed-Loop Model Coupling','The Simulink air-stepping','Added eight-BPA air-step result and 0.31 s mechanical / 1.938 s RG mismatch.','No entrainment or imported-body control established.'),
('E06','W2L mass changes, support, and robustness','Perturbation Robustness and Pathway Deletions','The coupled W2L','Preserved recorded ten-condition result; defined criterion; separated nominal support, actual force, and re-tuning. Included failed basin outcomes 9/12 and 7/12.','No harness-free walking; source/result provenance requires reconciliation.'),
('E07','Per-target optimization not applied','Perturbation Robustness and Pathway Deletions','The parameter-loading','Static audit found 55/61 winner keys discarded; removed acceptance of per-target optimization claims.','Priority ZCode repair and effective-weight replay required.'),
('E08','Deletion and perturbation figure overclaims','Perturbation Robustness and Pathway Deletions','Exploratory network-deletion','Excluded defective figure; qualified recorded left counts; identified unapplied switches and missing recovery measure.','Need corrected bilateral logs and actual edge-removal checks.'),
('E09','New Di Russo modeling','Perturbation Robustness and Pathway Deletions','A separate controller derived','Reported the completed 120-trial negative stability search as an adaptation.','Not a successful paper replication; no new simulation run here.'),
('E10','Comparison requested by advisor','Capabilities and Remaining Limits','The resulting contribution','Avoided claiming a universally unmet library need. Asked Ben for Hunt/Fletcher library identity.','Comparison remains open; do not score as satisfied.')]
matrix=[]
for id_,concern,section,anchor,change,limit in rows:
 matrix.append(dict(id=id_,concern=concern,section=section,file='chapters/60-sns.tex',line=line(E/'chapters/60-sns.tex',anchor),text=a.splitlines()[line(E/'chapters/60-sns.tex',anchor)-1],change=change,limit=limit,draft_score=2,acceptance_score=0,meaning='2 = local revision with evidence; 0 acceptance = unaccepted, not certified clear or scientifically validated'))
(V/'self_scoring_matrix.json').write_text(json.dumps(matrix,indent=2),encoding='utf-8')
report='''# Revision E — SNS chapter review

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

'''
for row in matrix:
 report+=f"- **{row['id']} — {row['concern']}**, § {row['section']}, `chapters/60-sns.tex:{row['line']}`: {row['change']} **Open:** {row['limit']}\n"
report+='\n## Evidence pointers\n\n'
for p,t in [(S/'reports_20261016/w2l_wean/_w2l_wean_eval.py','if k in T.KV'),(S/'w2l_mujoco/test_w2l_ground.py','gains = dict'),(S/'reports_20261025_persyn/_poster_panels.py','cycR=cycles'),(S/'reports_20261025_persyn/_poster_panels.py','kneeL_deg='),(S/'w2l_mujoco/test_w2l_ground.py','log["q"][i]'),(S/'w2l_mujoco/test_w2l_ground.py','Fz = max'),(M/'sns_verify_from_json.m','all 410 neurons at'),(M/'logs/sns_run_deng2.log','max|dV_ext|')]: report+='- '+ref(p,t)+'\n'
report+='\nFull source paths and SHA-256 hashes are in `review/evidence_manifest.json`. The Di Russo run log records no passing stage-2 trial. The older connectome audit is used as historical analysis, not as a current certification after subsequent changes.\n'
report+=f'''\n## Build and visual review

Fresh MiKTeX pdflatex → bibtex → pdflatex → pdflatex build: **{result['pages'][0]} pages**, zero errors, zero unresolved references/citations, and zero missing figures. There are **{result['overfull_boxes']} overfull boxes** across the full dissertation; these remain a document-wide layout limitation. The rendered Chapter 6 pages (printed 119 onward) were checked for clipping and figure placement. See `review/contactsheet.png` and `review/build_result.json`. No new model evaluation or independent reproduction of simulation outcomes was performed.

Output: `Revision_E_local.pdf`. Reviewable source: `chapters/60-sns.tex`. Build directory: `{B.name}`. The chapter preserves existing figures except the explicitly excluded perturbation figure. Other chapters remain the Revision D baseline, including any inherited issues outside this SNS pass.

New uncommitted v4 demo figures and their generator appeared during this session. They were inspected selectively and left untouched, outside the compiled E chapter. The knee-reflex v4 image reports a maximum 3.6° difference from its NumPy twin; it should not be called exact agreement. Approval/interpretation of these evolving figures remains pending.

## Evidence questions sent to Ben

1. Are there newer original-paper protocol comparisons, especially Deng perturbations and Li closed-loop contact feedback?
2. Are there matched feedback-on/off or BPA-model comparisons beyond the committed reports?
3. Is there a corrected perturbation figure or raw bilateral log on easteregg2?
4. What are the exact Hunt/Fletcher libraries that the advisor requested for comparison?
5. Does easteregg2 have uncommitted evaluator/runner changes that apply the individual conductances? Recover exact source and output hashes together.

These questions remain pending. E scopes the claims to available evidence and does not invent missing results. `ZCODE_HANDOFF.md` specifies repairs and the smallest useful checks before additional optimization.
'''
(E/'REPORT.md').write_text(report,encoding='utf-8')
handoff='''# ZCode handoff — Revision E

Read `Revision_E/REPORT.md`, `review/parameter_plumbing_audit.json`, `review/self_scoring_matrix.json`, and Ben's replies before changing prose. The chapter's point is to establish what the tools let us test about neural control, feedback, and mechanical differences. Completing a diagram, import, search, or build is not acceptance of that scientific claim.

## First: recover the exact running source

Run `hostname` first and use only that machine's AGENTS.md entry. Paths below are repository-relative. On easteregg2 the checkout is normally `D:\\GitHub\\Bipedal_Robot`; do not copy desktop environment paths. Preserve local/uncommitted work and record source hashes, backend versions, model hash, actual command, and output hashes. The fetched desktop baseline is a988508a. A newer easteregg2 runner may explain the discrepancy; do not overwrite it with this baseline.

Read the `dissertation`, `animatlab`, `mujoco`, `matlab`, `neuroscience`, `sns-toolbox`, and `opensim` skills from the local ZCode_Skills repo/junctions. Use relevant configured plugins where they provide direct source or model access. Read `spinal/DESIGN.md` before runner/builder edits. Preserve the s3k anchor and its existing contracts; separate W2L and Di Russo paths. Do not equate an adaptation with original-paper reproduction.

## Priority 1: parameter plumbing before optimization

At the reviewed commit, `reports_20261016/w2l_wean/_w2l_wean_eval.py:evaluate` keeps only keys already present in `T.KV`. The selected w25 file has 61 fields; **55 are dropped**, leaving six legacy knobs. `w2l_mujoco/test_w2l_ground.py:run_phase` constructs family-level gains and does not forward the individual per-target fields, despite their existence in `w2l_cpg/build_w2l_net.py`. The sensitivity file's many exact-zero effects therefore cannot be treated as evidence that those pathways are unimportant.

1. Confirm the exact source used on easteregg2. If its loading differs, save that version with its original outputs first.
2. Make each intended key explicit in accepted input, builder forwarding, saved winner, and reload. Reject unknown search keys loudly; never silently discard them.
3. For every searched key, perturb that key and prove the expected compiled edge conductance changes, while unrelated edges retain their values. Record edge IDs, signs, units, and source-specific conductance interpretation. Distinguish encoding gains from synaptic conductances. Shared L/R parameters must be labeled reductions.
4. Save the effective compiled configuration, not only the intended JSON. Reload in a fresh process and compare signatures. Check default settings against the retained baseline.
5. Replay the existing w25 and w10 candidates with the exact old effective configuration and corrected intended configuration. Do not silently replace earlier records or claim earlier scores apply to corrected dynamics. Only then decide whether additional search is needed.

## Priority 2: deletion and perturbation evidence

`reports_20261025_persyn/_poster_panels.py:cycles` always uses `heelF[:,0]`; the right count calls the same function on the same matrix. Give each leg its own trace. Count after the settled window with the same threshold/merge conventions as the evaluator. Its `kneeL_deg` uses hip column 0 and converts degrees again; log q is already degrees in JNTS order, knee_L index 1. Correct names and units rather than interpreting the large numbers physically.

Before every deletion, verify the target edges exist and that their effective compiled weights become zero. Restore the baseline between cases. Several current deletion keys do not reach the builder. Do not infer redundancy from an unapplied deletion. Ib's `ibnA` encoding removal is distinct from deleting specified synaptic edges.

Recover saved raw logs first. If absent, rerun only the necessary baseline/deletion/stimulus/push conditions after the plumbing checks. Record bilateral contacts, joints, support force, forward displacement/speed, phase, duty, fall/finite status, and post-perturbation windows. Define recovery tolerance and horizon before interpreting the push. State the 0.15 s pulse duration, force direction, stimulus amplitude, and 1 s duration. Show the actual rather than nominal harness contribution. Do not treat nominal 10% wsup as measured 10% total support.

The mass-change evaluator independently counts both legs and its ten records remain useful, but match them to the effective controller/source. Its driver inherits other runner defaults; audit `pushN` and every perturbation switch, including whether a default 40 N pulse was present. Record unperturbed controls explicitly. Preserve body/inertia changes and required physiological joint/contact gates. This is W2L, not the 92-actuator BPA robot model.

## Priority 3: explain what the additional work establishes

- Deng-derived Simulink: ten bursts / 1.938 s period after one pulse is an intrinsic-rhythm test under fixed tau_h. The log also shows about 63 mV max waveform differences. Diagnose alignment, solver, and equations before calling it a dynamic replica. Do not replace source tau_h(V) in other literature models without labeling that change.
- Full 410-neuron comparison: all final voltages at 0.3 s and nine short traces pass; extend to cycle statistics and closed-loop behavior separately. Short agreement does not prove long-term transfer.
- Li: 58/98 topology with graded adaptation replaces original spiking stance cells. Scripted antiphase contacts and correlation -0.271 are an adaptation check; verify the original closed-loop experiment separately.
- Air-step: fixed pelvis, no contact, eight BPA elements; mechanics about 0.31 s versus RG 1.938 s. Establish neural-to-mechanical entrainment and feedback effects before calling a curriculum stage or gait validated.
- Di Russo: committed stage-2 search has 120 failed trials and best objective -100. Preserve that negative result. Inspect validity gates, muscle mappings, balance and feedback equations, receptor delays, stance/contact, and stage-1 excitation fit before spending more optimization. Do not weaken gates to manufacture a pass.
- Hunt/Fletcher library comparison remains pending identity/source from Ben. Compare editable neuron/synapse coverage, sensory/actuator coupling, mechanical interfaces, and actual verification. Avoid an unsupported novelty claim.
- Newly produced local v4 demo figures are not in E. Preserve their files and Ben's figure instructions. Knee-reflex v4 reports max 3.6° discrepancy; label approximation and show timing/solver evidence. Do not call those twins exact.

## Work allocation and acceptance

Ben prefers GLM-5.3 / GLM-5.3-Flash for code, modeling, testing, and optimization, and ChatGPT Astra for final prose. Use the smallest evidence-recovery and plumbing checks first; heavy runs belong on easteregg2. This handoff authorizes no new external message, Overleaf edit, commit, or push. Ben will initiate ZCode on the other machine.

Return one evidence packet mapping each claim to source version, model, neural variant, effective parameter file, protocol, raw outputs, criterion, result, and limitation. Include failed runs and substitutions. Then update the local E draft with precise results. Keep repeated clarity concerns unaccepted until Ben/advisor can explain the goal and result back without build-order jargon. Every dissertation revision ends in a clean build plus line-based compliance evidence.
'''
(E/'ZCODE_HANDOFF.md').write_text(handoff,encoding='utf-8')
header='''# Revision E — SNS review, 26 October 2026

Built the local Revision E SNS chapter against fetched commit `a988508a`. **Critical finding:** 55 of 61 per-target winner fields are discarded by the committed evaluator, and the runner does not forward them. The perturbation figure also duplicates left-heel counts into its right column; several nominal deletions are not applied. These are implementation/evidence issues, not just wording issues. A newer easteregg2 source may differ and needs recovery with its outputs.

Revision E replaces the confusing five-task capability table with gaps, reasons, tests, results, and limits in prose. It incorporates Ben's 26 October per-leg PF / shared four-neuron RG clarification and MN-feedback wording, recent Simulink tests, recorded mass-change outcomes, failed conductance robustness probes, and the negative Di Russo stability search. The defective perturbation figure is explicitly excluded and its source retained. Historical review scores and Ben's replies are preserved. Clarity is still pending Ben/advisor acceptance.

Files: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/Revision_E_local.pdf`, `REPORT.md`, `ZCODE_HANDOFF.md`, and `review/self_scoring_matrix.json`. The report includes exact section/line/text evidence, source hashes, numerical changes, and build results. Overleaf was read only; no commit, push, new modeling run, or external message was made. Other chapters are the Revision D baseline.

---

'''
p=R/'CHATGPT_REPORT.md';old=p.read_text(encoding='utf-8');
if not old.startswith('# Revision E — SNS review'): p.write_text(header+old,encoding='utf-8')
p=R/'ZCODE_HANDOFF_20261023.md';old=p.read_text(encoding='utf-8')
if not old.startswith('# Priority update — Revision E'):
 p.write_text('# Priority update — Revision E\n\nRead `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/ZCODE_HANDOFF.md` first. It supersedes earlier SNS evidence assumptions: 55/61 winner keys are dropped in the fetched code; right-side deletion counts duplicate the left. Recover any newer easteregg2 implementation before repairing or replaying. Do not start a larger search before effective-parameter and actual-deletion checks.\n\n---\n\n'+old,encoding='utf-8')
print(json.dumps({'pdf':str(E/'Revision_E_local.pdf'),'build':result,'compliance_rows':len(matrix),'conservation':audit},indent=2))
