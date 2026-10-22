from pathlib import Path
import json,re,shutil,datetime,subprocess
R=Path(r'D:\Github\Bipedal_Robot');H=Path(__file__).resolve().parent
N=R/'Notes/advisor_review_20261022';N.mkdir(exist_ok=True)
for src in ['review_snapshot.json','review_anchored.json','highlights.json','source_comparison.json','applied_figure_patches.json','layout_patches.json']:
    if (H/src).exists():shutil.copy2(H/src,N/src)
report=R/'CHATGPT_REPORT.md';prompt=R/'Notes/new_session_prompt_20261022.md'
for p in [report,prompt]:
    backup=N/(p.stem+'_before_this_pass.md')
    if not backup.exists():shutil.copy2(p,backup)
previous=json.loads((N/'advisor_inventory.json').read_text(encoding='utf-8')) if (N/'advisor_inventory.json').exists() else []
def key(file,e):
    return (file,e['messages'][0]['id']) if e['type']=='comment' else (file,e['time'],e['change'])
priorids={key(g['file'],e):e['report_id'] for g in previous for e in g['entries'] if 'report_id' in e}
groups=json.loads((H/'review_final.json').read_text(encoding='utf-8'))
original=json.loads((H/'review_snapshot.json').read_text(encoding='utf-8'))
oldhi={(e['file'],e['position']):e['highlight'] for e in json.loads((H/'highlights.json').read_text(encoding='utf-8'))}
hiids={key(g['file'],e):oldhi.get((g['file'],e['position'])) for g in original for e in g['entries'] if e['type']=='comment'}
newhi={(e['file'],e['position']):e['highlight'] for e in json.loads((H/'highlights_final_additions.json').read_text(encoding='utf-8'))}
high={}
for g in groups:
    src=H/'background_final.tex' if g['file']=='15-background.tex' else H/'overleaf_before/chapters'/g['file']
    if not src.exists():src=H/'overleaf_before'/g['file']
    s=src.read_text(encoding='utf-8')
    for e in g['entries']:
        pos=e['position'];e['line']=s[:pos].count('\n')+1;e['context']=s[max(0,pos-100):pos+220]
        high[(g['file'],pos)]=newhi.get((g['file'],pos)) or hiids.get(key(g['file'],e))
counts=[]
for g in groups:
    entries=[e for e in g['entries'] if e['author']=='ajh26']
    if entries:counts.append((g['file'],sum(e['type']=='comment' for e in entries),sum(e['type']=='tracked edit' for e in entries)))
comments=sum(x[1] for x in counts);edits=sum(x[2] for x in counts)
introCommentMessages=sum(len([m for m in e['messages'] if m['author']=='ajh26']) for g in groups for e in g['entries'] if e['type']=='comment')
out=[f'''# Advisor review and revision handoff

Updated from the live **Bolen-Dissertation** Overleaf Review → Overview on 22 October 2026. This was a review of advisor contents and their anchors, not another full-document quality review. The advisor is actively reviewing; this is a snapshot, not a claim that no later comments exist.

## Current inventory

**{comments} open advisor threads ({introCommentMessages} advisor messages) and {edits} pending advisor tracked edits.** Ben's replies are explicitly separated below. Resolved comments and Ben-only appendix/main entries are excluded from the advisor totals. The previous 36-thread/8-edit inventory is superseded and preserved in Notes/advisor_review_20261022/CHATGPT_REPORT_before_this_pass.md. Existing SNS-01..15, DB-01..26, INTRO-01, BG-01, and BG-02 identifiers are retained.

| File | Advisor threads | Pending advisor edits |
|---|---:|---:|
''']
out.extend(f'| {f} | {c} | {t} |\n' for f,c,t in counts)
out.append('''
Line numbers and context refer to the fresh pre-upload source snapshot in tmp/advisor_review_20261022/overleaf_before/chapters, not the differing local ProofFinal chapters or PDF line numbers. A long paragraph occupies one source line. Positions are the rendered review panel's character offsets. Exact highlighted text was captured from the editor where available; a context window is explicitly labeled when no highlight was returned. Thread IDs are retained in the JSON evidence. Source offsets can drift after concurrent edits. Tracked edits remain pending; recording them does not accept them.

## Ben's current controlling instructions

These supersede conflicting older campaign plans and the earlier ZCode prompt:

1. Reduce the dissertation, especially SADb and SNS. Quality and a clear argument take precedence over more text. Do not bury a cut by adding an equally long appendix.
2. SNS does **not** require Methods/Results/Discussion headings. Lead with Ben's PhD-time combination of the Deng CPG and Li walker in AnimatLab: stepping-like movement in air, but no walking achieved by hand tuning. Use that limitation to motivate the knowledge/capability gaps in Ben's replies below. Follow with brief results of new software tools and concrete verification showing what works and what gap is closed. Future work belongs mainly in Chapter 7.
3. Preserve the advisor's active review. Supply narrowly targeted changes wherever possible; prefer a few paragraphs or one section to replacement of entire chapters. Log moves/deletions, retain labels and citations, and do not overwrite Overleaf from stale local chapters.
4. Bring Chapter 4 into line with verified recent torque-paper work from GitHub, distinguishing new measurements, uncorrected estimates, completed reductions, and predicted results.
5. Preserve Chapter 3's isometric-force figures and MATLAB figures. Reports and prompts stay local; **never upload CHATGPT_REPORT or handoff/administrative files to Overleaf.**

## Advisor comments and tracked changes
''')
for g in groups:
    file=g['file'];entries=[e for e in g['entries'] if e['author']=='ajh26']
    if not entries:continue
    out.append(f'\n## {file}\n')
    # Stable IDs for the older background edits; new material receives BG-03 onward.
    bgnew=max([int(v.split('-')[1]) for k,v in priorids.items() if k[0]=='15-background.tex']+[2])+1
    for i,e in enumerate(entries,1):
        if file=='60-sns.tex':id=f'SNS-{i:02}'
        elif file=='50-database.tex':id=f'DB-{i:02}'
        elif file=='10-introduction.tex':id=f'INTRO-{i:02}'
        elif key(file,e) in priorids:id=priorids[key(file,e)]
        elif 'Changed: requirement' in (e.get('change') or ''):id='BG-01'
        elif 'Deleted:' in (e.get('change') or '') and 'these tools were' in e['change']:id='BG-02'
        else:id=f'BG-{bgnew:02}';bgnew+=1
        e['report_id']=id
        out.append(f"\n### {id} · line {e['line']} · {e['type']} · {e['time']}\n\n")
        hi=high.get((file,e['position']))
        if hi:out.append('**Highlighted source:**\n\n```latex\n'+hi+'\n```\n\n')
        else:out.append('**Context near the review anchor (not an exact highlight):**\n\n```latex\n'+e['context']+'\n```\n\n')
        if e['type']=='comment':
            for m in e['messages']:
                who='Advisor' if m['author']=='ajh26' else "Ben's reply (author context, not an advisor instruction or independent verification)"
                out.append(f"**{who} — {m['time']}:**\n\n"+'\n'.join('> '+line for line in m['text'].splitlines())+'\n\n')
        else:
            change=re.sub(r'^(edit|add_circle|delete)\s*','',e.get('change') or '').replace('\nshow less','')
            out.append('**Pending advisor edit:**\n\n```text\n'+change+'\n```\n\n')
out.append('''
## What the advisor's new notes require

The introduction has newly annotated undefined references to numbered gaps, unexplained spiking-neuron mirrors, an overly specific roadmap, and a weak connection to the Biomimetic Humanoid Project. Name the gaps directly and present them in chapter order. Explain Ben's work as a continuation of the humanoid project, rather than inserting adjacent laboratory projects without a clear connection.

The new background notes target the motor drawing, the graphical abstract's usefulness, Gait2392's name, the moment-arm discussion, joint-center modeling, run-on sentences, and the chronology of the master's model. State what is compared and how it supports the present work. Replace abstractions such as “first-class,” “definitions are stressed,” and “written against” with concrete descriptions. Do not simply add explanatory paragraphs after every “huh?”; replace or remove the unclear statement and its unnecessary detail.

SADb should start with the finished resource: what it contains, its organization, and who can use it. Then give only the methods needed to reproduce the curation, adapting the advisor's Cochrane suggestion without pretending it was a registered systematic review. Verify the two conference presentations, AI-assisted workflow, decision authority, and actual inclusion criteria from evidence. Ben's role should be explicit. Drop the import chronology and operational record counts unless needed to interpret the final corpus.

For SNS, the old suggestion to compress AnimatLab to one generic paragraph is superseded by Ben's requested causal narrative. The account must identify what Ben combined, what the air-stepping result shows, what failed under hand tuning, and why the new tools were needed. A compact capability/evidence table can replace repeated inventories. Verification should identify the expected behavior, the test or comparison, the observed outcome, and its limit. Do not expand this into a catalogue of all experiments.

## GitHub torque update checked in this pass

Fetch succeeded. The active branch is KneeTestSetup_BenBo_stw. During the pass HEAD and origin both reached **65b05aae161d6a4e87f561a38e57dbae90ec6f95** (22 October, 12:00 Pacific), committed by another session/user; this agent did not commit or pull over local changes. Re-fetch at the start of ZCode work because other machines remain active.

- **fd1c13f4**, 22 October 02:22: updates Knee_Torque_Test/ProofFinal/Main.tex:540–542 with the six-point extensor cam-route sweep, explicitly pending reaction-point geometry and adjoint torque reduction. Read as manuscript status, not independently verified final torque values. Chapter 4 should not continue to say no full angle sweep exists if the test record confirms this sweep.
- The **same manuscript still says “full sweep ... pending” at line 556 and “first ... points ... near full extension” at line 560**. Resolve this internal inconsistency against the actual test sheets before adapting it. Do not copy both claims.
- **7862e4e2**, 22 October 04:31: adds ExtTest_3_01..03 records, updates Results_table_20mm.xlsx, Knee_Extensor_20mm.m, and Ext20mm_compare_20261022.png, and adds the new BiPam_pulley2_20261022.m plus sweep33/escalation records. These are newer than the six-point paragraph; inspect their provenance before selecting the current validation status.
- The 21 October commits also update torque and margin surfaces and the transmission doctrine. PULLEY_SWEEP_DIRECTIVE_20261021.md is the current scope: nBPA and tackle parts each span 1–3, re-solve actuator sizes, account for moving-via load and Bowden compliance, and audit OpenSim path artifacts. The existence of CSVs/mats is not proof that every configuration is feasible or implemented in the simulation model.

The dissertation's 30-torque.tex also differs between live Overleaf and local ProofFinal. Neither was overwritten in this pass. Integrate in a small evidence-backed patch against fresh Overleaf text, not by replacing the chapter with the torque manuscript.

## Figure synchronization in this pass

Uploaded the 11 revised PDF assets to their existing figs/10-introduction, figs/15-background, figs/60-sns, and figs/70-overall folders. Kept editable PPTX, EPS, alt text, and generators local. Chapter 3 assets were not uploaded or edited.

Applied seven exact-fragment changes to figure references, sizing, captions, and the sentence introducing the functional neuron view. No whole chapter was uploaded. A subsequent layout check led to six further exact-fragment caption/sizing changes to the two combined figures: shorter captions and 0.94-textwidth panels retain figure text above 10 pt. The separate patches are recorded in the JSON evidence. Advisor tracked edits and comment status were not accepted, rejected, or resolved.

No report, prompt, or administrative file was uploaded to Overleaf; the fresh pre-upload ZIP contains no chatgpt_report file. The source snapshot revealed six divergent chapters (introduction, background, torque, SADb, SNS, and Discussion). Local drafts and live advisor edits must be reconciled, not treated as interchangeable. A fresh source ZIP is evidence at its download time, not a permanent canonical version.

## ZCode execution contract

Use Notes/new_session_prompt_20261022.md. It replaces the stale broad task list with a focused reduction/revision pass, explicit prose and caption limits, the new SNS structure, small patch batches, and evidence requirements. The old prompt and report are preserved locally in Notes/advisor_review_20261022. No chapter reduction has been claimed completed in this pass: ZCode is assigned that heavy revision work.
''')
report.write_text(''.join(out),encoding='utf-8')
(N/'advisor_inventory.json').write_text(json.dumps(groups,ensure_ascii=False,indent=2),encoding='utf-8')
print('COUNTS',counts,comments,edits,'MESSAGES',introCommentMessages)
