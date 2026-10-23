from pathlib import Path
import re,json,difflib
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/final_sweep_20261023'
entries=json.loads((S/'recorded_comments.json').read_text());assert len(entries)==124
def source(file,phrase):
    p=C/file
    for n,line in enumerate(p.read_text(encoding='utf-8-sig').splitlines(),1):
        if phrase in line:return f'[{p.name}:{n}](<{p.as_posix()}:{n}>)'
    raise ValueError((file,phrase))
mapping={}
def setids(prefix,nums,file,phrase,disposition):
    for n in nums:mapping[f'{prefix}-{n:02}']=(file,phrase,disposition)
setids('INTRO',range(1,13),'chapters/10-introduction.tex','Neural control exposes','The introduction names the evidence and model-coupling gaps directly; quantified corpus coverage is scoped to the database chapter. Undefined and awkward gap labels are removed.')
setids('INTRO',[13],'chapters/10-introduction.tex','Four connected' if 'Four connected' in (C/'chapters/10-introduction.tex').read_text() else 'four connected needs','Force, torque, afferent evidence, and simulation gaps appear in chapter order, matching the objectives and organization.')
setids('INTRO',[14],'chapters/10-introduction.tex','The present work builds','The humanoid-project paragraph connects the earlier robot, knee, muscle-routing work, and collaborators to this dissertation.')
setids('INTRO',range(15,35),'chapters/10-introduction.tex','Dissertation Organization','Replaced the detailed roadmap and numbered-gap shorthand with named chapters, brief purposes, and appendix pointers. Tracked wording is incorporated or its obsolete sentence is removed.')
setids('BG',[1,2],'chapters/15-background.tex','Simulation Tooling','The rewritten background incorporates the earlier tracked-edit intent and removes obsolete wording; simulation tools and their dissertation role are described briefly.')
setids('BG',[1],'chapters/15-background.tex','first improvement','The advisor’s replacement of “requirement” with “improvement” is retained in the BPA-force-model bridge.')
setids('BG',[3],'chapters/15-background.tex','fig:actuator_technology_comparison','Verified the replacement schematic: rotor and stator are concentric. The existing vector asset is unchanged in this sweep.')
setids('BG',[4],'chapters/15-background.tex','Braided Pneumatic Actuators','The confusing graphical abstract is absent from the background; the protected actuator-paper chapter is untouched.')
setids('BG',range(5,10),'chapters/15-background.tex','23 degrees','Gait2392 is identified by 23 degrees of freedom and 92 musculotendon actuators. Its moment arms are explicitly model-derived, not collected robot measurements; jargon is removed.')
setids('BG',range(10,15),'chapters/15-background.tex','Moment arms relate','The signed fixed-axis moment-arm definition and OpenSim length-sensitivity method are distinguished; Gait2392 planar knee kinematics are specified. The repository digression is removed.')
setids('BG',range(15,19),'chapters/15-background.tex',"In my master's work",'The master’s work, its comparison, and the subsequent placement optimization are described in separate sentences with citations.')
setids('BG',range(19,24),'chapters/15-background.tex','fig:masters_routing_and_arms','The figures have stated comparison purposes; the geometry explanation and caption are shortened. No unsupported cross-product definition remains.')
setids('BG',range(24,31),'chapters/15-background.tex','I adjusted BPA origins','Removed file-saving chronology; named origins, insertions, and via points; explained contraction/clearance constraints and the human-versus-robot comparison.')
setids('BG',[31,32,33],'chapters/15-background.tex','fig:rybak_cpg_circuit','The deleted literature-synthesis and three-connectome catalogue remain absent. The compact Rybak/Shevtsova pair is retained and placed together under Ben’s latest instruction.')
setids('BG',range(34,51),'chapters/15-background.tex','Computational neuroscience expresses','Replaced the awkward paragraph and its tracked fragments with a general account of computational neuroscience. CPG-specific inventories and the misplaced dissertation-controller sentence are removed.')
setids('DB',[1,2,3,4,5,6,7],'chapters/50-database.tex','The knowledge gap was','The resource is introduced first, followed immediately by the evidence gap and what SADb supplies. The binary laterality statement is removed; remaining unknowns are acknowledged without claiming a complete physiological specification.')
setids('DB',[8],'chapters/50-database.tex','The collection began','The review-led starting point is in Curation Methods; no import chronology is restored.')
setids('DB',[9],'chapters/50-database.tex','AI-Assisted Triage','Both Neuroscience presentations remain identified: the database in 2023 and AI-assisted triage in 2024, with their existing citations.')
setids('DB',[10,11,12,13,14],'chapters/50-database.tex','The working base in Airtable','Defines curation notes and feedback records with an example; removes first-class/canonical jargon; distinguishes editable Airtable from the current full-record app distribution.')
setids('DB',[15,16],'chapters/50-database.tex','Use in Neural-Controller Design','The controller-use implications are at the end, rather than interrupting the resource introduction.')
setids('DB',[17,18,19,20],'chapters/50-database.tex','Curation Methods','Sources, identity matching, relevance criteria, source-text grounding, and human decision authority are explicit. The resource is not misrepresented as a registered or exhaustive Cochrane review.')
setids('DB',[21,22],'chapters/50-database.tex','The working base in Airtable','The opening explains the finished resource, five record types, organization, and what curation notes contain. The bubble screenshot is in the main chapter.')
setids('DB',[23,24,25,26],'chapters/50-database.tex','I retained responsibility','Ben receives credit and decision authority; inclusion scope is explicit. Working curation and read-only distribution are separated, without the confusing “curation tool, not a reader” phrase.')
setids('SNS',[1,2,3,6,7,9,10,14],'chapters/60-sns.tex','Within this project, five','The causal opening leads from Ben’s bilateral walker to all five combinations in his reply. It states why those missing connections matter, what tools were built, complementary simulator limits, and which combinations remain unverified. “Build order” and “tooling” shorthand are removed.')
setids('SNS',[4,5,8],'chapters/60-sns.tex','The AnimatLab Bipedal Walker','The account clearly concerns a simulated model. The obsolete “Preliminary Results” heading is gone, and the rewritten passage incorporates the tracked-edit intent.')
setids('SNS',[11],'chapters/60-sns.tex','fig:pipeline','The existing replacement vector pipeline schematic has separated columns and labels; the main chapter uses that corrected asset.')
setids('SNS',[12,13],'chapters/60-sns.tex','Model Conversion and Actuator Coupling','The “laboratory stands on tools” sentence and the unnecessary SADb inventory are deleted.')
setids('SNS',[15],'chapters/60-sns.tex','Compared with the earlier SNS models','A functional comparison with Hunt/Young is provided: shared neural formulation versus added actuator, encoding, and mechanical interfaces. No unperformed historical-library performance benchmark is claimed.')
rows=[];inventory={}
for e in entries:
    assert e['id'] in mapping,e['id']
    file,phrase,disposition=mapping[e['id']];link=source(file,phrase)
    q=next((q['text'] for q in e['quotes'] if q['speaker'].startswith('Advisor')),'Tracked edit (original text retained in the ledger)')
    q=re.sub(r'\s+',' ',q).lstrip('> ')[:160].replace('|','/')
    rows.append(f"| {e['id']} | {e['kind']} | {q} | {disposition} | {link} |")
    inventory[e['id']]={'kind':e['kind'],'disposition':disposition,'source':file,'evidence':phrase}
out='''# Advisor-comment disposition sweep — 2026-10-23

Read the expanded live Overleaf Review Overview and compared it with CHATGPT_REPORT.md. Inventory matches: Introduction 34, Background 49, SADb 26, SNS 15; total **124 advisor items = 73 comment threads + 51 pending tracked edits**. The ledger also retains Ben’s replies separately. No new open advisor items were detected in this sweep. The three Ben-only appendix/main tracked deletions and the general-knee-transform comment were checked separately: obsolete text is absent, and Appendix C gives the angle-dependent transform.

“Addressed” here describes the local edited text. No Overleaf comment was resolved and no tracked edit was accepted/rejected. It does not assert advisor approval. Historical resolved threads remain in the existing saved archive; they were not reopened in this sweep. Deleted source context and the full original review entries are retained in the before snapshot and recorded_comments.json.

SNS-15 is a scoped functional comparison, not a performance benchmark of all earlier libraries. The five requested SNS combinations are not all declared fully achieved; remaining integration limits are explicit. Those are scientific scope limits, not silently claimed completed experiments.

| ID | Kind | Recorded request | Local disposition | Evidence |
|---|---|---|---|---|
'''+ '\n'.join(rows)+'\n'
(C/'COMMENT_DISPOSITIONS.md').write_text(out,encoding='utf-8');(S/'comment_dispositions.json').write_text(json.dumps(inventory,indent=2))
diff=[]
for old in sorted((S/'before').glob('*.tex')):
    new=C/'chapters'/old.name
    diff.extend(difflib.unified_diff(old.read_text(encoding='utf-8-sig').splitlines(True),new.read_text(encoding='utf-8-sig').splitlines(True),fromfile=str(old),tofile=str(new)))
(S/'final_sweep.diff').write_text(''.join(diff),encoding='utf-8')
root=R/'CHATGPT_REPORT.md';text=root.read_text(encoding='utf-8-sig')
header='## Final Revision C sweep — 2026-10-23'
if header in text:
    print('Updated all 124 dispositions; the root report already contains the final-sweep record')
    raise SystemExit(0)
text+='\n\n'+header+'''

The expanded live Overleaf Review Overview was checked against this ledger. The inventory remains **73 advisor comment threads and 51 pending advisor tracked edits**: Introduction 34 items, Background 49, SADb 26, SNS 15. Ben’s replies and the four Ben-only appendix/main items remain separate. No new open advisor items were found. No Overleaf source, comment resolution, or tracked-edit disposition was changed.

Every advisor ID now has a local disposition and source-line evidence in `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C/COMMENT_DISPOSITIONS.md`. The original expanded entries are retained in `review/final_sweep_20261023/recorded_comments.json`. These are editorial dispositions, not claims of advisor approval.

This pass makes the SNS knowledge gap match all five combinations in Ben’s reply, keeps the interconnected bilateral RG construction in the causal opening, adds an early SADb evidence-gap statement, and puts the bubble chart in the SADb main chapter. The current app payload was independently read: it now includes all 1,393 records, so the old paper-only app limitation is superseded. The separate monthly scraper remains under development.

The completed physical-test plots from the latest commits are inserted; no spreadsheet sanity-check torque values are used. Four full fit tables move to Appendix C, with a compact mean-error summary in Chapter 4. The cat caption is shortened, and Figures 2.8 and 2.9 share one page. The active bipulley run is not reported as completed. Original sources and the protected force chapter/figures remain unchanged by this pass.
'''
root.write_text(text,encoding='utf-8')
print('Wrote all 124 dispositions and appended the fresh review audit to CHATGPT_REPORT.md')
