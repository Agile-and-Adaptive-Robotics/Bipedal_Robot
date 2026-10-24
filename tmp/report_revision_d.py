from pathlib import Path
import json,shutil,re
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D';E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023'
changes=json.loads((D/'review/changes.json').read_text(encoding='utf-8'));matrix=json.loads((E/'review_matrix.json').read_text(encoding='utf-8')); build=json.loads((D/'review/build_result.json').read_text(encoding='utf-8'))
shutil.copy2(Path(build['build'])/'main.pdf',D/'Revision_D_local.pdf')
def link(p,line=None):return f'[{p.name}{":"+str(line) if line else ""}](<{p.as_posix()}{":"+str(line) if line else ""}>)'
intro=f'''# Local prose revision D — 23 October 2026

**18 targeted edits across Introduction, Background, Database, and SNS are complete locally. Overleaf is unchanged.** Ben clarified that ChatGPT Astra is his intended final prose model. No model was launched, no provider was configured, and no scheduled monitoring was restarted.

Working copy: {link(D/'Revision_D_local.pdf')}. Actionable code/model handoff: {link(D/'ZCODE_HANDOFF.md')}. Exact before/after changes: {link(D/'review/changes.diff')}; per-edit record: {link(D/'review/changes.json')}.

## Source and scope

The five captured Overleaf chapter snapshots supplied the reviewed source. Uncaptured chapters, main.tex, bibliography, and initial assets came from the existing Revision_C copy solely to provide local build support. This is **not a fully reconciled current Overleaf project** and is not an upload bundle. ZCode must merge the targeted changes into fresh source, not replace live chapters wholesale. ProofFinal and Revision_C were not edited by this pass.

The existing flexor panel Flx20mm_meetsHuman.pdf was missing from the Revision_C support assets. It was copied unchanged from tmp/overleaf_live_20261023 for the build. No figure was generated or replaced. The database pivot caption's “animal preparations” was explicitly changed to “experimental preparations”; all other figure blocks were preserved, and all figure inclusions, labels, and citation commands remain conserved.

## Changes and self-check

| Edit | File and final source line | Change / reason | Local self-score |
|---|---|---|---|
'''
for i,c in enumerate(changes,1):
 score=3 if c['file']=='10-introduction.tex' or c['reason'].startswith(('Define archived','Explain citation direction','Complete a sentence','Remove literal')) else 2
 c['local_score']=score
 # Anchor candidates are for cross-navigation, not proof that every nearby comment is answered.
 before=(D/'review/before'/c['file']).read_text(encoding='utf-8'); at=before.find(c['before']);start=before[:at].count('\n')+1 if at>=0 else None;end=start+c['before'].count('\n') if start else None
 c['matrix_context_ids']=[r['id'] for r in matrix if r['file']==c['file'] and not r['collapsed'] and start and start<=r['line']<=end]
 intro+=f"| D-{i:02d} | {link(D/'chapters'/c['file'],c['new_line'])} | {c['reason']} | {score}/4; acceptance unverified |\n"
intro+='''
Scores of 3 concern a directly checked local clarification; scores of 2 retain an unresolved contextual, reader-understanding, or implementation-evidence dependency. No score of 4 was assigned. Earlier 0/4 clarity failures remain in the historical matrix; these local revisions do not erase them. Context IDs in changes.json link nearby matrix entries but do not claim every nearby short comment is addressed.

## Verification

- Fresh pdflatex → bibtex → pdflatex ×2: **257 pages, zero errors, zero undefined references, zero missing figures.** Full logs and build_result.json are in review/ and the named build directory.
- Figure paths/counts, label lists, and citation command lists are unchanged in all four edited files. Bibliography was copied byte-for-byte and never re-encoded.
- All numerical measured/model results are preserved. In SNS, duplicated mentions of the 92-actuator count and Gait2392 name were shortened; their substantive values remain. Literal list ordinals 1–5 became an enumerate list. Numeric-token differences are logged in conservation.json; no result value was silently changed.
- Background's unsupported claims about thousands of overnight combinations, hundreds of cycles per condition, unmatched throughput, and statistical power were explicitly replaced by a qualitative statement of repeatable simulation tests. These were unverified performance claims, not benchmark results. The assertion that identical models are independent biological subjects was explicitly removed. Old wording is retained in changes.json/diff.
- Citation directions checked against SADb_audit/build_citation_graph.py:8 and SADb_audit/app/build_app.py:320–322,725–744. App-only archive status checked against the app header/archive control. These are local code checks, not claims of an exhaustive live-base schema audit.
- Database opening and SNS task-list page renders were inspected: text was legible, with no clipping on those sampled pages. **26 overfull boxes remain in the full build.** The document is not declared layout-clean.

## Rubric disposition and unresolved work

Gate 1: no new external prose copied; existing citation commands and publication notes preserved. Gate 2: changed prose checked for banned boilerplate, unexplained labels, and reviewer-comment leakage; no new reviewer text or agent comments placed in .tex. Gate 3: numerical results conserved, scoped claims checked as above; detailed implementation claims still require ZCode's artifact audit. Gate 4: direct clarifications made; broad clarity concerns remain provisional until Ben/advisor reads the result. **No overall scientific-validation GO is claimed.**

ZCode's handoff specifies evidence to retrieve for actual conversion defects, five-system integration status, bilateral organization/authorship, source-model replication, earlier-library comparison, and database links. It requires existing evidence before new runs, source-backed edge tables before optimization, and one baseline plus one published stimulus/deletion protocol before expanding a replication campaign. GLM-5.3 / Flash are Ben's preferred execution models there. ChatGPT Astra is reserved for bounded final prose polishing after evidence is settled.

All source edits are local and reviewable. No Overleaf changes, uploads, replies, accept/reject actions, resolutions, commits, or pushes were performed.

---

'''
(D/'REPORT.md').write_text(intro,encoding='utf-8')
(D/'review/changes.json').write_text(json.dumps(changes,indent=2,ensure_ascii=False),encoding='utf-8')
prior=(R/'CHATGPT_REPORT.md').read_text(encoding='utf-8');(R/'CHATGPT_REPORT.md').write_text(intro+prior,encoding='utf-8')
# Keep original feedback/score intact; attach the actual local attempt and evidence separately.
for row in matrix:
 attempts=[dict(edit=f'D-{i:02d}',file=str(D/'chapters'/c['file']),line=c['new_line'],score=c['local_score'],status='Local revision; reader acceptance unverified',evidence=str(D/'review/changes.json')) for i,c in enumerate(changes,1) if row['id'] in c['matrix_context_ids']]
 if attempts:row['local_revision_D']=attempts
(E/'review_matrix.json').write_text(json.dumps(matrix,indent=2,ensure_ascii=False),encoding='utf-8')
assert (R/'CHATGPT_REPORT.md').read_text(encoding='utf-8').endswith(prior)
print(json.dumps({'changes':len(changes),'pdf':str(D/'Revision_D_local.pdf'),'handoff':str(D/'ZCODE_HANDOFF.md'),'matrix_rows_with_local_attempt':sum('local_revision_D' in x for x in matrix)},indent=2))
