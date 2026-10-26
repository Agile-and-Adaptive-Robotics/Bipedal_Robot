from pathlib import Path
import json,shutil,difflib
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review/targeted_reconciliation';V.mkdir(exist_ok=False)
p=E/'chapters/60-sns.tex';b=p.read_text(encoding='utf-8');shutil.copy2(p,V/'before_60-sns.tex')
deleted=' This suspended result concerns the original shared-RG controller in Figure~\\ref{fig:animatlab_model_architecture}(B).'
assert deleted in b;a=b.replace(deleted,'');p.write_text(a,encoding='utf-8')
(V/'local_deletion.diff').write_text(''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='before/60-sns.tex',tofile='after/60-sns.tex')),encoding='utf-8')
patches=[
dict(id='P00',action='KEEP DELETED; no new text',comment='ajh26, 23 October 10:51 am: Huh? Just delete all this?; Ben 26 October 5:20 am: deleted',old=' This suspended stepping result preceded the later feedback and bilateral-RG extensions; it does not demonstrate walking with those additions.',new='',status='Deletion verified live. Do not restore this sentence or a paraphrase. Retain the rest of the outcome paragraph unless Ben directs otherwise.'),
dict(id='P01',action='Replace only these final two sentences in Knowledge Gaps and the New Tools',comment='ajh26, 23 October 11:21 am: need a better word; second marker *',old='I developed the specific converted models, actuator interfaces, graphical-library elements, and circuit implementations below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.',new=r'I developed the model conversions, actuator-force calculations, Simulink blocks, and controller implementations described below. Table~\ref{tab:sns_capability} identifies the interfaces tested and the model integrations that remain incomplete.',status='Exact text read in current Overleaf. Table and label both present. Keep table while applying this patch; do not combine with the table-free local E wording.'),
dict(id='P02',action='Replace only item 4 wording; retain Ben’s block-programming description',comment='Ben 26 October 2:21 am tracked change and 5:15 am reply to huh?',old='I built a GUI for block-programming neural controllers that were attached mechanical models',new='I built a GUI for block-programming neural controllers coupled to mechanical models',status='Exact current text read live; grammar/coupling wording only.'),
dict(id='P03',action='Rename only the subsection title; preserve label',comment='Ben’s explicit heading instruction in chat',old=r'\subsection{Graphical Circuit Construction}\label{sec:sns_new_tools}',new=r'\subsection{Block Programming Neural Pathways}\label{sec:sns_new_tools}',status='Exact current heading read live. Preserve the existing label and comments.'),
]
(V/'patches.json').write_text(json.dumps(patches,indent=2),encoding='utf-8')
for row in patches:
 if row['id']!='P00':(V/(row['id']+'_replacement.tex')).write_text(row['new']+'\n',encoding='utf-8')
live={'read_date':'2026-10-26','project':'https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb','mode':'read-only DOM observation; show-more and comment-anchor navigation only; no source edits, accept/reject, replies, or resolves','deletion':{'author':'bbolen83','time':'26 October, 5:19 am','text':patches[0]['old'],'reply':'deleted','reply_time':'26 October, 5:20 am'},'terminology':{'author':'ajh26','time':'23 October, 11:21 am','comments':['need a better word','*'],'anchor_note':'Live character offsets 4587 and 4636; not line numbers or exact selected words. Ben identifies connections as the target. The paragraph contains those connections and connections demonstrated.'},'table':{'label':'tab:sns_capability','exists_live':True,'existing_explicit_reference':r'Table~\ref{tab:sns_capability} separates the implemented interfaces from incomplete mechanical integration.','uncited_mention':'The capability table identifies the connections demonstrated and those still incomplete.'},'versions':'The current live chapter differs substantially from local E. Local E is a review draft, not a drop-in replacement.'}
(V/'live_verification.json').write_text(json.dumps(live,indent=2),encoding='utf-8')
md='''# Targeted Overleaf edits — current coordination sheet

**Overleaf is the working document. Do not upload or paste entire local E chapters.** Local E/PDFs are reference drafts that contain proposed prose, newer figures, and supplemental evidence. They do not include the live review history or every author change. Whole-chapter replacement is not the application method.

## Coordination

- Ben owns the live document and author decisions. Preserve his tracked deletions and wording changes. A deleted passage stays deleted unless he explicitly reopens it.
- ChatGPT supplies bounded prose patches, linked to the advisor’s exact comment and existing text, with a reason and evidence limit.
- ZCode supplies source-verified modeling results and implementation details. Its ongoing experiments should arrive as one evidence packet, rather than parallel chapter rewrites. This file is a local handoff; no message was sent to ZCode.
- Apply each approved patch in place with tracking retained. Do not resolve comments or accept/reject existing changes as part of a prose patch. After application, check the associated thread still exists and its new location makes sense. Replacing a commented span can move its anchor; in-place editing is not a guarantee that an anchor remains unchanged.
- Match exact existing text, section label, and comment author/date. If the match is absent or changed, stop that patch and compare the new author text. Line numbers are supplementary navigation only.
- A large section can be a bounded patch when necessary, with explicit start/end text and a list of affected threads. A whole chapter is not the default patch unit.

## Three issues raised now

**Deletion:** the live tracked deletion is the sentence beginning “This suspended stepping result preceded…”. The remainder of the hand-tuning/outcome paragraph is still present. Ben replied “deleted.” The local E draft incorrectly put back a paraphrase, “This suspended result concerns…”. That paraphrase has now been removed locally. P00 records the instruction; nothing needs to be reinserted into Overleaf.

**Connections:** “framework” names the overall toolset. It is not a replacement for a particular signal exchange or mechanical coupling. Use “synaptic connections” for circuit edges, “interfaces” for exchanges between software/model components, and “coupling” for their operation together. In the criticized paragraph, replace the vague two concluding sentences with P01. Preserve the rest of the live paragraph for the next bounded review; this patch does not silently endorse or rewrite its five-system inventory.

**Capability table:** the current Overleaf source does contain a table labeled `tab:sns_capability`, and the Software Tools paragraph references it. The earlier paragraph calls it “the capability table” without a reference. P01 fixes that sentence while retaining the live table. Local E deliberately removed the table. Do not mix the table-free E paragraph with the table-present Overleaf passage. Removing or redesigning the table, if chosen later, must be a separate coordinated patch affecting its mentions and label together.

## Exact patches

'''
for row in patches:
 md+=f"### {row['id']} — {row['action']}\n\n**Comment:** {row['comment']}\n\n**Find exactly:**\n\n```tex\n{row['old']}\n```\n\n"
 if row['new']:md+=f"**Replace with:**\n\n```tex\n{row['new']}\n```\n\n"
 md+=row['status']+'\n\n'
md+='''## Previously requested edits still need bounded application

Figure 6.1’s (A)/(B) labels, Ben’s original shared-RG image, the later per-side RG chronology, and the conversion-repair appendix are prepared locally. They are not automatically present in Overleaf. Their source and evidence are in `review/figure61_correction/`. Apply the figure block, the chronology paragraph, and the appendix insertion as separate patches, preserving the existing document and review threads. The main-text appendix reference and appendix section must be applied together to avoid an unresolved reference. No chapter swap is required.

The larger clarity proposals in `review/clarity_pass/` remain proposals. A row marked “editorial treatment supplied” meant the local draft was edited, not that the live advisor comment was resolved. Those rows must not be used as a certificate of live compliance. This sheet takes precedence for application strategy.

No Overleaf source or review action was changed during this reconciliation. `live_verification.json` records what was actually observed. The replacement `.tex` files here contain fragments only, not complete chapters.
'''
(V/'OVERLEAF_TARGETED_PATCHES.md').write_text(md,encoding='utf-8')
# Put the application rule where each collaborator is likely to look first.
notice='''# Application rule — targeted patches, not chapter replacement

The live Overleaf document is the baseline. Do not upload/paste Revision E chapter files wholesale: they lack current tracked deletions, changes, and comment anchors. Use `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E/review/targeted_reconciliation/OVERLEAF_TARGETED_PATCHES.md`. It contains exact-text replacements, comment/date evidence, and the current table/deletion decisions. Ben’s deleted sentence must stay deleted. Supply modeling evidence separately from prose patches; do not launch parallel chapter rewrites. No Overleaf changes or inter-agent messages were made by this reconciliation.

---

'''
for p in [R/'CHATGPT_REPORT.md',R/'ZCODE_HANDOFF_20261023.md',E/'ZCODE_HANDOFF.md',E/'REPORT.md']:
 old=p.read_text(encoding='utf-8');p.write_text(notice+old,encoding='utf-8')
print('Prepared P00–P03, preserved live edits, and removed the local reintroduced qualifier.')
