from pathlib import Path
import re,json,difflib,shutil,hashlib
from collections import Counter
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review/clarity_pass'
def find(file,anchor):
 lines=(E/'chapters'/file).read_text(encoding='utf-8').splitlines();i=next(i for i,l in enumerate(lines,1) if anchor in l);return i,lines[i-1]
audit={};diff=''
def cites(s):return set(k.strip() for m in re.findall(r'\\cite\w*(?:\[[^\]]*\])*\{([^}]*)\}',s) for k in m.split(','))
for old in sorted((V/'before').glob('*.tex')):
 new=E/'chapters'/old.name;b=old.read_text(encoding='utf-8');a=new.read_text(encoding='utf-8')
 diff+=''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='before/'+old.name,tofile='after/'+old.name))
 labels=lambda t:re.findall(r'\\label\{([^}]+)\}',t)
 images=lambda t:re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',t)
 nums=lambda t:Counter(re.findall(r'\d+(?:\.\d+)?',t))
 audit[old.name]={'labels_preserved':set(labels(b))<=set(labels(a)),'added_labels':list(set(labels(a))-set(labels(b))),'figures_preserved':images(a)==images(b),'citation_keys_preserved':cites(b)<=cites(a),'removed_numeric_tokens':list((nums(b)-nums(a)).elements()),'added_numeric_tokens':list((nums(a)-nums(b)).elements()),'sha256':hashlib.sha256(new.read_bytes()).hexdigest()}
 assert all(audit[old.name][k] for k in ('labels_preserved','figures_preserved','citation_keys_preserved')),(old.name,audit[old.name])
(V/'changes.diff').write_text(diff,encoding='utf-8');(V/'conservation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
catalog=json.loads((R/'tmp/clarity_catalog.json').read_text(encoding='utf-8'))
maps={
'M-015':('15-background.tex','\\caption[Rybak two-level','Prior Ben reply reports removal of the confusing graphical abstract. Current retained schematics have explicit captions; no further figure was deleted. Exact original figure anchor remains uncertain.'),
'M-017':('15-background.tex','Moment arms relate','Verified the earlier replacement: signed fixed-axis moment arm and tendon-excursion definition are stated directly; ambiguous “stressed” wording is absent.'),
'M-018':('15-background.tex','Knee modeling requires','Verified the explicit distinction between prescribed Gait2392 knee kinematics and the four-bar robot knee; ambiguous “written against” phrasing is absent.'),
'M-026':('15-background.tex','OpenSim \\citep','Verified that the Simulation Tools paragraph describes OpenSim geometry and torque analysis directly, rather than a repository digression.'),
'M-048':('15-background.tex','Figure~\\ref{fig:decerebrate_cat} documents','Replaced the accounting metaphor with the literal purpose of the biological-preparation figure and scoped model relevance to validation.'),
'M-050':('15-background.tex','The experiments this research program','Retained the earlier explicit explanation of repeated controlled model interventions and its warning that repeated models are not independent biological subjects.'),
'M-074':('50-database.tex','\\caption[SADb citation','Rewrote the caption to state exactly which papers cite Grillner and which papers Grillner cites; retained counts 57 and 13.'),
'M-082':('50-database.tex','Records were retained','Explained inclusion, archive decisions, and consistent classification tags without vocabulary shorthand.'),
'M-083':('50-database.tex','Records were retained','Defined a sensory pathway as a sensory signal influencing neural activity or motor output, with the stance/swing example.'),
'M-092':('50-database.tex','When I began this work','Verified that the resource gap and the source-linked collection built to address it are stated in the chapter opening.'),
'M-101':('60-sns.tex','For a course project','Retained Ben’s per-leg PF organization, shared four-neuron RG, and MN-feedback clarification.'),
'M-102':('60-sns.tex','For a course project','Retained “We added” for the collaborative feedback work.'),
'M-103':('60-sns.tex','We adjusted','Used collaborative attribution for hand tuning.'),
'M-104':('60-sns.tex','We adjusted','Kept the concrete suspended-motion outcome and its support limit; removed build-order shorthand from the main rationale.'),
'M-105':('60-sns.tex','I subsequently expanded','Linked bilateral circuit construction directly to the subsection reporting its coupled tests.'),
'M-106':('60-sns.tex','The suspended motion left','Stated the unresolved loading/contact question and why adjustable support and recorded feedback were needed.'),
'M-114':('60-sns.tex','The project lacked','Replaced the five-system inventory with the concrete experimental capability gap.'),
'M-120':('60-sns.tex','I addressed this gap','Used “framework,” model conversion, force calculation, and named experimental controls in place of vague systems terminology.'),
'M-121':('60-sns.tex','I addressed this gap','Reviewed the same formerly awkward inventory as connected prose.'),
'M-122':('60-sns.tex','Existing tools supply','Explained the distinct roles of OpenSim, MyoConverter, MuJoCo, SNS-Toolbox, Simulink, Simscape, and AnimatLab.'),
'M-124':('60-sns.tex','My first conversion','Explained the actual robot-first conversion failure and the human-reference/robot-routing distinction.'),
'M-125':('60-sns.tex','The human-model conversion','Retained concrete small-muscle, frame, patella-path, and contact repairs, rather than an unexplained audit claim.'),
'M-126':('60-sns.tex','The source represents patella','Retained actual fixes supplied by Ben and recovered from source; advisor’s hypothetical examples are not treated as events.'),
'M-127':('60-sns.tex','For the BPA-routed robot','Specified inputs (length, pressure), Chapter 3 force calculation, applied force, sensory return, and omitted pressure dynamics.'),
'M-129':('60-sns.tex','The actuator and torque models','Stated the chapter’s purpose and implemented contribution in the opening paragraph.'),
'M-130':('60-sns.tex','The framework is evaluated','Explained circuit, interface, and coupled-model tests in terms of the questions they answer.'),
'M-131':('60-sns.tex','Existing tools supply','Explained why environments are complementary; graphical imports are distinguished from controlled coupled models.'),
'M-132':('60-sns.tex','The project lacked','Placed the gap and implementation before the methods; kept implications in Capabilities and Remaining Limits.'),
'M-133':('60-sns.tex','The framework addresses','Connected the delivered interfaces to the specific experiments they enable.'),
'M-134':('60-sns.tex','\\caption[Simulation-to-robot','Retained the previously corrected pipeline asset; checked rendered caption/placement, with no replacement figure introduced.'),
'M-135':('60-sns.tex','The MuJoCo simulations allow','Defined support, curriculum, stimulation, mechanical perturbation, recorded outputs, and search role directly.'),
'M-136':('60-sns.tex','Existing tools supply','Expanded operational purpose in the support/search subsection while retaining the roles and limits of existing tools.'),
'M-137':('60-sns.tex','Existing tools supply','Libraries were never supplied; Ben does not want to request them. No comparison or universal novelty claim is invented.')}
ledger=[]
for row in catalog:
 if row['id'] not in maps:raise RuntimeError('Unmapped comment '+row['id'])
 f,anchor,treatment=maps[row['id']];ln,text=find(f,anchor)
 ledger.append(dict(id=row['id'],file=f,line=ln,text=text,advisor_comments=[m for m in row['messages'] if m['author']=='ajh26'],ben_replies=[m for m in row['messages'] if m['author']!='ajh26'],original_highlight=row['highlight'],original_anchor_status=row['anchor_status'],treatment=treatment,status='Unavailable comparison; no request to Ben' if row['id']=='M-137' else 'Editorial treatment supplied; human clarity review pending',draft_score=2,acceptance_score=0))
(V/'comment_compliance.json').write_text(json.dumps(ledger,indent=2),encoding='utf-8')
result=json.loads((E/'review/build_result.json').read_text());assert result['errors']==0 and result['undefined_references']==0 and result['missing_figures']==0
B=Path(result['build']);shutil.copy2(B/'main.pdf',E/'Revision_E_local.pdf')
lines=['# Submission clarity pass — Revision E','',f'Edited 17 principal passages across six chapters, plus collaborative attribution and a direct results cross-reference. Reviewed {len(ledger)} clarity-related comment threads from the preserved advisor matrix. This is a prose and consistency pass; no simulation was run, optimization started, or Overleaf edit made.','',
'The chapter now leads with its concrete contribution: a framework for connecting neural controllers to the human reference and separate BPA models, varying support and training stages, applying neural or mechanical interventions, and measuring their outcomes. It distinguishes capabilities that exist from a universal biological result or unsupported walking. Optuna is described as the implemented search tool, without claiming a global optimum or superiority.','',
'The numerical findings and their scientific limits remain. The corrected W2L candidate retains contact cycling under several deletions; Ib removal gives 0/0; nine of ten mass-change cases fail posture. These results are not relabeled robust walking. The framework contribution does not depend on resolving those control problems before submission.','',
'The Discussion and Conclusion now summarize this contribution consistently; the stale tonic-mismatch summary was replaced. The Introduction states the same scoped framework gap. Background prose explains neural/mechanical relationships directly. Database prose defines sensory pathway and citation direction.','',
'## Comment-by-comment treatment','',
'Historical character offsets and some frozen-source contexts are not exact current anchors. The new lines below identify the present treatment, not a claim that an old highlight was precisely reconstructed. M-015’s exact historical figure remains uncertain; no additional figure was deleted. M-137 is an unavailable comparison, consistent with Ben’s instruction not to request those libraries. Acceptance stays pending, rather than being awarded because a checklist or build passed.','']
for row in ledger:
 comments=' / '.join(m['text'] for m in row['advisor_comments']);lines.append(f"- **{row['id']}**, `{row['file']}:{row['line']}` — {comments}\n  {row['treatment']}")
lines+=['','## Conservation and build','',
'Every existing figure asset reference, label, and citation key in the six edited chapters was retained. One subsection label was added for experiment controls. No bibliography file was edited. The prior chapter sources are preserved in `before/`; exact prose edits are in `changes.json` and `changes.diff`, and token-level numerical changes in `conservation.json`. The removed numerical occurrence is a repeated trace-correlation threshold in the Discussion; the value remains in the Conclusion and detailed SNS results. No measured result or adopted parameter value was changed.',
'',f"Fresh build: **{result['pages'][0]} pages**, no errors, no unresolved citations/references, and no missing figures; **{result['overfull_boxes']} document-wide overflow warnings** remain. A prior attempt encountered a PDF file lock; a fresh build directory completed all four passes. Build: `{B.name}`. The PDF output is `Revision_E_local.pdf`.",
'','Rendered review covers the revised SNS chapter and edited Background/Database passages. The compact multi-panel perturbation figure retains small labels; its underlying full-size PDF remains available. No claim of publication-ready figure typography or resolution of every dissertation layout warning is made.',
'','The scientific evidence report and remote snapshot history remain in `REPORT.md` and `review/easteregg2_followup/`. ZCode’s ongoing modeling need not block review of this frozen local prose draft. New results require evidence review before being incorporated; they are not pulled into this pass automatically.']
(V/'CLARITY_PASS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
header='''# Latest pass — submission clarity

The current PDF and chapter sources include a submission-focused prose pass across Introduction, Background, Database, SNS, Discussion, and Conclusion. Read `review/clarity_pass/CLARITY_PASS.md` and its line-based comment ledger for the newest editorial changes. The scientific evidence history below remains preserved. No Overleaf edits or new model runs were made.

---

'''
p=E/'REPORT.md';old=p.read_text(encoding='utf-8');
if not old.startswith('# Latest pass — submission clarity'):p.write_text(header+old,encoding='utf-8')
root_header='''# Revision E submission clarity pass — 26 October 2026

The local PDF now leads with the concrete SNS framework: Gait2392 and separate BPA models, neural-controller alternatives, adjustable simulated harnesses, staged tuning, neural stimuli, external perturbations, recorded responses, and bounded Optuna searches. It makes no solver-superiority or global-optimum claim. Existing results, including posture failures and the scoped Ib finding, remain qualified.

Seventeen principal prose changes across six chapters align the Introduction, Background, Database, SNS, Discussion, and Conclusion. The recorded clarity-related comments have line-based treatments in `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/clarity_pass/comment_compliance.json`; `CLARITY_PASS.md` is the readable review. Earlier texts are preserved. Acceptance remains for Ben/advisor to judge. No Overleaf changes, model runs, commits, or pushes.

---

'''
p=R/'CHATGPT_REPORT.md';old=p.read_text(encoding='utf-8')
if not old.startswith('# Revision E submission clarity pass'):p.write_text(root_header+old,encoding='utf-8')
# Refresh the earlier evidence matrix's current line pointers without changing historical scores.
p=E/'review/self_scoring_matrix.json';matrix=json.loads(p.read_text(encoding='utf-8'))
anchors={'E01':'The project lacked','E02':'For a course project','E03':'My first conversion','E04':'A separate Simulink','E05':'The Simulink air-stepping','E06':'The corrected candidate was also','E07':'The W2L tests examine','E08':'With the stored individual','E09':'A separate controller derived','E10':'The present contribution','E11':'An activity-based reduction','E12':'The human-model feedback-tuning'}
for row in matrix:
 row['line'],row['text']=find('60-sns.tex',anchors[row['id']])
p.write_text(json.dumps(matrix,indent=2),encoding='utf-8')
print(json.dumps({'comments_reviewed':len(ledger),'build':result,'conservation':audit},indent=2))
