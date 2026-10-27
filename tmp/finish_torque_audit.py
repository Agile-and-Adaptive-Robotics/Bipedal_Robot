from pathlib import Path
import json,re,difflib,datetime,hashlib
root=Path(r'D:\Github\Bipedal_Robot')
folder=root/'Documentation/Reports and Papers/Knee_Torque_Test'
w=folder/'review_audit_20261026'
p=(w/'paper_Main_live.tex').read_text(encoding='utf-8')
d=(w/'dissertation_30-torque_live.tex').read_text(encoding='utf-8')
pa=(w/'paper_Main_after_GoF.tex').read_text(encoding='utf-8')
da=(w/'dissertation_30-torque_after_GoF.tex').read_text(encoding='utf-8')
entries=json.loads((w/'paper_review_entries_with_positions.json').read_text())
dentries=json.loads((w/'dissertation_torque_review_entries.json').read_text())
archive=json.loads((w/'paper_resolved_comments.json').read_text())
def clean(s): return re.sub(r'\s+',' ',s).strip()
def line_at(s,pos): return s[:max(0,int(pos or 0))].count('\n')+1
def line_of(s,needle):
    i=s.find(needle); return s[:i].count('\n')+1 if i>=0 else 'absent'
def section_at(s,pos):
    names=re.findall(r'\\(?:sub)*section\*?\{([^}]+)\}',s[:max(0,int(pos or 0))])
    return names[-1] if names else 'Preamble/front matter'
def current_line(s,pos):
    n=max(0,int(pos or 0)); start=s.rfind('\n',0,n)+1; end=s.find('\n',n)
    return s[start:end if end>=0 else len(s)]
def item_type(e):
    if 'entry-comment' in e['classes']: return 'Comment'
    if 'entry-insert' in e['classes']: return 'Insertion'
    if 'entry-delete' in e['classes']: return 'Deletion'
    return 'Replacement'
def tidy_item(s):
    return '\n'.join(l for l in s.splitlines() if l.strip() not in ['check','close','edit','delete','add_circle','more_vert','Accept change','Reject change','Resolve comment','More options','show less','Reply'])
def figure_blocks(s):
    found=[]
    for m in re.finditer(r'\\begin\{(figure\*?|table\*?|sidewaystable)\}',s):
        end=s.find('\\end{'+m.group(1)+'}',m.end())
        block=s[m.start():end+len('\\end{'+m.group(1)+'}')]
        found.append({'type':m.group(1),'line':line_at(s,m.start()),'labels':re.findall(r'\\label\{([^}]+)\}',block),'images':re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',block),'text':block})
    return found
manifest={'paper':figure_blocks(p),'dissertation':figure_blocks(d)}
(w/'figure_table_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
(w/'torque_documents_raw.diff').write_text(''.join(difflib.unified_diff(p.splitlines(True),d.splitlines(True),fromfile='paper_Main_live.tex',tofile='dissertation_30-torque_live.tex')),encoding='utf-8')
norm=lambda s:re.sub(r'\s+',' ',re.sub(r'%[^\n]*','',s)).strip()
paragraphs=lambda s:[x.strip() for x in s.split('\n\n') if len(x.strip())>100]
dp=paragraphs(d); align=[]
wordset=lambda s:set(re.findall(r'[A-Za-z]{3,}',s.lower()))
dsets=[wordset(q) for q in dp]
for paragraph in paragraphs(p):
    ws=wordset(paragraph)
    candidates=sorted(range(len(dp)),key=lambda i:len(ws&dsets[i])/max(1,len(ws|dsets[i])),reverse=True)[:3]
    score,best=max((difflib.SequenceMatcher(None,norm(paragraph),norm(dp[i]),autojunk=True).ratio(),dp[i]) for i in candidates)
    align.append({'paper_line':line_of(p,paragraph),'dissertation_line':line_of(d,best),'similarity':round(score,3),'paper_excerpt':paragraph[:180],'dissertation_excerpt':best[:180]})
(w/'paragraph_alignment.json').write_text(json.dumps(align,indent=2),encoding='utf-8')
gofdiff=''.join(difflib.unified_diff(p.splitlines(True),pa.splitlines(True),fromfile='paper_before_GoF',tofile='paper_after_GoF'))+'\n'+''.join(difflib.unified_diff(d.splitlines(True),da.splitlines(True),fromfile='dissertation_before_GoF',tofile='dissertation_after_GoF'))
(w/'applied_GoF_prose.diff').write_text(gofdiff,encoding='utf-8')
historical=['# Torque-paper resolved review archive','', 'Read-only export. These threads are resolved history, not a new action queue. Original UI date strings are retained without inventing missing years.','']
ben_history=[]
for i,t in enumerate(archive,1):
    historical += [f'## R{i:03d} — {t["file"]}','',f'Quoted text: {t["quote"] or "[no surviving quote in UI]"}','']
    for msg in t['messages']:
        historical += [f'**{msg["author"]} — {msg["date"]}**', '',msg['text'] or '[empty message]','']
        if msg['author']=='bbolen83' and msg['text']!='Marked as resolved': ben_history.append((i,msg))
(w/'historical_review_archive.md').write_text('\n'.join(historical),encoding='utf-8')

rows=['# ChatGPT torque-paper / dissertation report','',f'Completed {datetime.datetime.now().astimezone().isoformat(timespec="seconds")}. Machine: EB475WS4.','', '**Main findings:** the two documents contain materially different figures, data counts, historical-fit descriptions, and presentation of results. Several of Ben\'s current paper comments also apply to the dissertation. A blind chapter replacement or automatic transfer from either document would propagate errors.','', 'Paper: [Design Factors for Isometric Robotic Knee Torque](https://www.overleaf.com/project/622e9536b8ea1e0a7f167cf5). Dissertation: [Bolen-Dissertation](https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb), `30-torque.tex`.','', '## Scope and evidence','', '- Captured live source from both Overleaf editors. Line numbers below refer to the saved **pre-GoF-edit snapshots**, not PDF line numbers or character positions. Use section/figure labels to locate the current live passages.', '- Paper review overview: **31 active items** before this audit\'s GoF edit: 29 in Main.tex, one in BolenFrontiers22.bib, and one in latexmkrc. Main.tex has **11 comment threads** (nine Ben, two older advisor) and **18 tracked changes** (17 Ben, one advisor). The two other files each have a tracked insertion by Ben.', '- Dissertation torque overview: **13 active items** before this audit\'s GoF edit: four comment threads and nine tracked changes, all attributed to Ben.', '- Resolved paper history: **272 threads**, fully expanded and preserved, including **74 substantive messages by Ben**. Resolution events are not counted as substantive replies.', '- Only the subsequently authorized GoF prose was edited during this audit. Figures, tables, numerical results, pending review dispositions, and comments were not changed. No commits, pushes, model runs, or messages to ZCode.','', '## GoF: does the paper explain it?','', '**Yes, but briefly.** The paper\'s Identification Procedure (snapshot line 300) defined RMSE, FVU, and maximum absolute residual, and named the solver. The dissertation had a longer account at lines 310–314. Its one-to-one mapping of metrics onto offset, shape/amplitude, and phase was too strong.','', '**Ben\'s clarification is applied in both documents:** the metrics are interpreted together, rather than assigning one number to each type of curve disagreement. The revised explanation states what each metric computes, how their combined interpretation assesses torque magnitude and variation, and why multiobjective optimization retains trade-offs between overall agreement, comparison across tests, and local error. The Pareto front and held-out selection are explained. Existing algorithm settings and selection procedures were retained.','', 'For a fixed test, RMSE and FVU are transformations of the same squared residuals; their aggregation across tests of different scales can weight those tests differently. FVU does not exclude mean bias, and maximum residual is not insensitive to a uniform offset. This mathematical point is why the combined interpretation is appropriate.','', 'Exact applied prose: [gof_prose_patches.json](review_audit_20261026/gof_prose_patches.json). Whole-file readbacks matched precisely one bounded replacement in each document. [Applied diff](review_audit_20261026/applied_GoF_prose.diff).','', '## Current paper comments and corresponding dissertation actions','', '| ID / Ben\'s comment | Paper anchor | Dissertation finding / required action |','|---|---|---|']
comments=[
('P-C01','This figure got replaced with an incorrect one. The blue dots are where the back-calculated yellow dots should be.','Figure environment `fig:Ext20mm`, lines 506–510; `Figures/Ext20mm_cam_20261023.eps`.','Compare with `fig:extensor_physical_final`, lines 599–604, `extensor_20mm_physical.png`. The paper plots six blue measurements; the dissertation describes nine measurements plus yellow pressure-specific predictions. Audit data provenance before changing colors or transferring either result.'),
('P-C02','Don\'t use "holds". Call it data measurements.','Physical extensor results, especially line 504 and caption line 510.','Still applies: `sec:redesign_physical`, line 595, says nine/fifteen “holds,” and later text says “between holds.” Search both documents and captions for measurement uses; replace with data measurements. Do not alter the mechanical meaning of holding a load when describing the procedure.'),
('P-C03','Too big for page. Should probably have a femur and tibia shape in there to make it less confusing. We could reduce the amount of panels shown.','`fig:ExtensorRoute_redesigned`, lines 485–490; eight poses plus legend.','Dissertation uses a different composite under the same label, lines 551–567: cropped extensor_cam_route.pdf panels and a hardware photograph. Review leg outlines and final-page readability in each version; the paper\'s eight-pose replot should not silently replace the dissertation composite.'),
('P-C04','Candidate to move to supplementary material if we are over the Figure/Table limit.','Blank-line anchor directly before `fig:RouteVariants`, lines 474–480; confirmed by comment navigation and surrounding source.','No equivalent RouteVariants figure in the torque chapter. Its route/pulley material is distributed between the chapter and Appendix C. Treat this as a conditional paper-length decision, not an instruction to delete dissertation results. Main paper currently has 11 figure floats and four table floats; journal limits were not assumed.'),
('P-C05',"There\'s no /Chi_{4} lol XD. Reply: Make a 2x2 panel figure with the legend in the last panel. Title panels 10 mm, 20 mm, 20 mm.",'`tab:AdoptedXi` header, line 457, has a χ₄ column.','Dissertation table at lines 528–540 already has χ₀–χ₃ only. Remove the empty paper χ₄ column in a future table edit. **Ben explicitly confirmed that the reply\'s 2×2 instruction applies to fig:FlxGrp**, not this table. Preserve these as two separate actions.'),
('P-C06','Leaves the page and is rotated the wrong way.','`tab:FlxPinGof`, especially line 356; 14 columns in a sidewaystable.','Dissertation separates pre-optimization and adopted-term tables and already carries Ben\'s tracked table-layout changes (lines 350–407). Fit and orientation must be checked in the paper\'s article layout separately. Do not rotate until the final reader-facing page is inspected.'),
('P-C07',"I don\'t think Hunt Equation was calculated correctly here.",'`fig:FlxPin`, caption line 334; `FlxPin_allmodels_20261023.eps`.','Dissertation pre-optimization Figure `fig:FlxPin_preopt` and its GoF table, lines 341–368, use the same Hunt comparison. Audit the calculation and score mask in both documents. Existing low-strain undefined-point exclusions do not establish that the equation itself was evaluated correctly.'),
('P-C08','This is all incorrect','Initial extensor paragraph, line 368.','Corresponding historical identification paragraph remains in the dissertation at line 425. Reconcile with saved fit records, intended parameter stage, held-out fold, and the terms actually used to generate the displayed curves. Do not relabel old single-bracket numbers as adopted two-stage results.'),
('P-C09','This is incorrect in the values and which tests were training and validation. The whole paragraph is incorrect.','Initial flexor paragraph, line 328.','Near-verbatim corresponding paragraph remains in the dissertation at line 375, including 11.9 mm, 102/10.6 kN/m, initial-fit errors, and a single validation test. Both need the same record-level reconciliation. Adopted two-bracket fold is two training/two validation; identification history must be distinguished from current figure membership.'),
('P-C10','Advisor, 6 May: how does the shorter length result in this?','Old thread currently navigates to the conclusion near line 556. Original highlighted phrase no longer survives.','Anchor is historically uncertain; retain the thread. The force/resting-length and wrapping explanations may be relevant, but do not claim the original question answered without the original quoted context.'),
('P-C11',"Advisor, 6 May: Can this be connected with Jon\'s work?",'Old thread also navigates to conclusion line 556; original anchor no longer survives.','Identity of Jon\'s specific work is not established by this thread. Preserve as unresolved historical context; do not invent a reference or retrofit a claim.')]
for c in comments: rows.append('| '+' | '.join(x.replace('|','\\|') for x in c)+' |')
rows += ['', '## Differences beyond the comment anchors','', '| Issue | Paper | Dissertation / disposition |','|---|---|---|',
'| Cam-extensor data count and result | Six measurements; only one exceeds the human target at the measured, lower pressures. | Nine measurements and a broader “most of the tested range” summary. Current `Knee_Extensor_20mm.m` contains nine entries. This needs a shared data export and recalculated comparisons, not a prose-only merge. |',
'| Biomimetic extensor GoF | Line 407 reports adopted RMSE 1.1702, FVU 0.6681, max 2.2088; line 468 also distinguishes older terms. | Line 465 reports 0.8920 / 0.3396 / 1.4188, while line 540 separately reports about 1.17 / 0.67. The dissertation is internally inconsistent about the parameter/evaluation version. |',
'| Flexor panel C | `FlxGrp_20261023.eps` is a single merged three-panel asset; caption fully describes the redesigned route. | A/B in FlxGrp.eps plus a small separate C from Flx20mm_meetsHuman.pdf. C caption ends at “The redesigned two-BPA flexor route.” Carry over a complete explanation with the confirmed 2×2 replot. |',
'| Original 10 mm extensor photograph/curve | One A/B figure, photo A and torque plot B. | Separate torque plot and pass-through photograph with two labels. A future merge must retain reference aliases and update panel references. |',
'| Adopted χ₃ scope | Flexor row is a dash, extensor 0.158; table describes identification. | Both rows show 0.158, caption says it was identified on extensor and used in flexor redesign, not fitted to pinned flexor. Explain identification versus deployment instead of silently harmonizing cells. |',
'| Pulley claims | Route procedure says mono is one tackle and bi one tackle per joint; says all record runs were unity transmission. | Chapter describes a later size-resolved transmission study; Appendix C has pulley torque/margin surfaces. Distinguish the older unity-routing runs from the later study. Earlier reviewed rulings describe mono/bi by the dominant DoF used in the route calculation; verify implementation before porting the paper explanation. |',
'| Deleted unity-gain paragraph | Ben deleted the long “The knee-flexor run quantifies…” paragraph at 2:20 pm. | The first paragraph of `sec:disc_pulley`, line 589, still leads with unity transmission and 81.23→85.74 N m / 0.803 N m. Mark for editorial reconciliation with Ben\'s deletion; do not restore the long paragraph from history. |',
'| Design-factor wording | Ben changed budget→allowance, added the biomimicry trade-off, and added isokinetic/isotonic future-testing discussion. | Chapter contains related anatomical-bounds and routing trade-offs, but not the same explicit design-factor synthesis. Selectively carry over meaning, correcting “biomimcry” and the repeated “in in”; do not copy tracked prose mechanically. |',
'| Removed caveat block | Ben deleted a long limitations paragraph, then scoped the remaining statement to this stand/configurations/actuator families. | Related limits are distributed in route and physical-test sections. Preserve the deletion as editorial history; retain necessary, accurately scoped limitations without importing the deleted block wholesale. |',
'| Whole-body capability survey | Section `sec:capability_survey` reports a 19-actuator 5×5-domain comparison (66.32→75.30% torque coverage; 97.05→95.23% active-force coverage). | Not a corresponding standalone torque-chapter section. Appendix C has separate route-search/capability-support material. Verify whether the same campaign is represented before cross-referencing; do not conflate these coverage statistics with route feasibility counts. |',
'| Torque/knee-axis notation | Existing paper equations retain several lowercase moment/torque forms and plain knee subscripts. | Ben explicitly requests uppercase torque M and bold knee-axis subscript k in dissertation comments. Apply consistently to both manuscripts in a separate notation pass, distinguishing torque from moment-arm variables and knee subscripts from summation indices. |',
'| Composite mechanics illustration | Paper has the frame figures but no `fig:xiProjection` composite. | Ben asks to split `fig:xiProjection` and put alpha on the circle. This is a dissertation-specific replot, not a paper change to force through without checking its different figure set. |',
'| GoF exposition | Shorter before this audit. | Longer but overassigned metrics to specific error causes. **Both revised live during this audit** to describe their combined effect and multiobjective selection. |',
'| Paper conclusion grammar | Ben\'s tracked opening currently reads “started with an improved … model and actually delivers about a joint.” | Preserve the recorded edit, but polish the resulting sentence before submission; it is not a suitable sentence to transplant into the dissertation. |',
'| Compile warnings | GoF revision compiles with 0 errors and 20 warnings: duplicate fig:extensor_misalign, small-font substitutions, and their summary. | Record/reconcile these remaining warnings separately from scientific revisions. No new labels or citation keys were introduced by the GoF edits. |',
'', '## Every active paper item: exact text and source location','', 'The original UI character positions are preserved in JSON as anchor evidence. They are converted to snapshot line numbers here; they are not presented as line numbers themselves. The excerpt is the containing source line, not a claim to recover the complete highlighted range. Short tracked fragments are grouped by their actual location, so “is→was” is not treated as a free-standing scientific instruction.','']
for i,e in enumerate(entries,1):
    file=e['file'] or 'unknown file'; ismain=file.startswith('Main.tex')
    rows += [f'### P-I{i:02d} — {item_type(e)} — {file}', '', '```text',tidy_item(e['text']),'```','']
    if ismain:
        ln=line_at(p,e['position']); sec=section_at(p,e['position']); anchor=current_line(p,e['position'])
        rows += [f'Anchor: paper snapshot line **{ln}**, **{sec}**.','',f'> {anchor[:1800] or "[blank-line anchor; surrounding figure context is recorded in paper_comment_anchors.json]"}','']
    else: rows += ['File-level item; source line not established from the current Main.tex editor. Full inserted content is preserved in the overview JSON. The bibliography insertion is empty in the displayed UI; do not infer an entry.','']
rows+=['## Every active dissertation torque item','']
for i,e in enumerate(dentries,1):
    rows += [f'### D-I{i:02d} — {item_type(e)} — line {line_at(d,e["position"])}','', '```text',tidy_item(e['text']),'```','',f'Containing section: **{section_at(d,e["position"])}**.','',f'> {current_line(d,e["position"])[:1800]}','']
rows += ['## Historical comments and Ben replies','',f'The full, untruncated archive contains {len(archive)} resolved threads and {len(ben_history)} substantive messages by Ben: [historical_review_archive.md](review_audit_20261026/historical_review_archive.md). It preserves each surviving quoted anchor, author, date string, reply, and resolution event. The following index identifies every substantive Ben message; these are historical records rather than new requirements.','']
for ri,msg in ben_history: rows.append(f'- **R{ri:03d}, {msg["date"]}:** {clean(msg["text"])}')
rows += ['', '## ZCode delivery and verification','', 'Use [ZCODE_REPLOT_REQUESTS.md](review_audit_20261026/ZCODE_REPLOT_REQUESTS.md). Return bounded `.tex`/figure changes under the existing machine/batch review inbox. Include exact generating script, source-data identifiers, parameter row/fold, plotted series export, and a rendered preview. Do not overwrite live chapters or accept/reject tracked changes.','', 'Evidence files: [raw document diff](review_audit_20261026/torque_documents_raw.diff), [paragraph alignment](review_audit_20261026/paragraph_alignment.json), [figure/table manifest](review_audit_20261026/figure_table_manifest.json), [paper review items](review_audit_20261026/paper_review_entries_with_positions.json), [dissertation review items](review_audit_20261026/dissertation_torque_review_entries.json).','', 'The diff is a complete source diff; structural/template differences are expected because one file is a standalone article and the other a dissertation chapter. The curated differences above identify the actionable scientific and editorial issues. No current Overleaf asset hash comparison was obtained; matching local filenames are generator leads, not proof of byte-identical live images.','']
(folder/'ChatGPT_reportTorque.md').write_text('\n'.join(rows),encoding='utf-8')
replots=r'''# ZCode: torque replot and reconciliation requests

Ben requested these changes in the live torque-paper/dissertation reviews. Read `../ChatGPT_reportTorque.md` and the saved live sources first. Do not alter the published isometric-force chapter. Do not upload whole chapters, resolve comments, accept/reject tracked changes, or launch new optimization campaigns. Replots should use existing data and the adopted correction terms. Ben commits.

## P1 — Correct the physical torque data/series before styling

Paper `fig:Ext20mm` uses `Figures/Ext20mm_cam_20261023.eps`; Ben says blue measurements occupy positions expected for the yellow back-calculated series. The dissertation `fig:extensor_physical_final` uses `extensor_20mm_physical.png` and has nine measurements, versus the paper's six. Do not fix this by swapping colors. Export a table containing angle, load, load-cell angle convention, pressure, reaction-point geometry, measured torque, predicted torque at measured pressure, 620 kPa prediction, and human target. Classify pressure-specific force-model predictions separately from a hybrid calculation using measured muscle length/moment arm.

Generator lead: `tmp/paper_bridge/figscripts/gen_Ext20mm_cam_20261023.m`, especially lines 41–137 (six-point arrays/reduction) and 140–159 (blue scatter). Current source lead: `Code/Matlab/Mesh_Optimization/Results/Knee_Extensor_20mm.m`, especially lines 582–653: `K_deg` contains nine values. Data: `Testing_Data/2026_06_Festo/Results_table_20mm.xlsx`, `ExtTest20mm_2`. Follow AGENTS.md's load-cell convention and the current measured rig constants. Audit old flexor-template sign-flip machinery rather than assuming it is authoritative. Recalculate both documents' claims from the same declared dataset. Replace “holds” with data measurements in descriptive prose/captions.

## P1 — Hunt equation, fit stages, and GoF tables

Ben questions the Hunt equation in paper `fig:FlxPin`; the dissertation `fig:FlxPin_preopt`/pre-fit table carry the same issue. Generator lead: `tmp/paper_bridge/figscripts/gen_FlxPin_allmodels_20261023.m`, lines 55–72. It interpolates `ForceStrainForFit.mat` after replacing NaNs in z with zero, rather than directly calling the equation. Cross-check against `Code/Matlab/Functions/HuntEq.m`, correct diameter, pressure units, strain versus normalized strain, maximum contraction, hysteresis setting, force direction, and joint moment sign. Compare the source equation with the lookup table before declaring either wrong. Export the valid-point mask and report both per-model-domain and common-domain scores if the domains differ.

Ben rejects the initial flexor paragraph (paper line 328 / dissertation line 375) and initial extensor paragraph (paper line 368 / dissertation line 425). Use `minimizeFlxPin10_results.mat` row 77 and `minimizeExt10mmX3_results.mat` row 32 as adopted records; distinguish them from explicitly historical single-bracket fits. Follow AGENTS.md's fold bookkeeping: two training/two validation/one unused flexor test, and three training/three validation/three unused extensor tests. Map nominal test IDs to the actual resting lengths printed in the manuscript. Recompute scores and memberships before replacing prose, not the reverse.

## P1 — Confirmed fig:FlxGrp layout

Ben explicitly confirmed in this chat: use a **2×2 tiled layout**, three data panels titled **10 mm**, **20 mm**, **20 mm**, and a dedicated legend in tile 4. Preserve A/B/C margin letters. The last panel is the redesigned two-BPA route; state that distinction in the caption despite the repeated 20 mm title. Keep measured points, corrected predictions, original/improved-force baselines, reduced-pressure points, and human targets correctly identified. Use shared axes where comparison is meaningful, but do not clip data merely to share a range. Deliver one coherent asset for each document, not the dissertation's tiny pasted C inset.

Generator lead: `tmp/paper_bridge/figscripts/gen_FlxGrp_3panel_20261023.m`, especially layout lines 184–277. Paper asset: `FlxGrp_20261023.eps`; dissertation currently combines `FlxGrp.eps` and `Flx20mm_meetsHuman.pdf`. Preserve all fifteen verified redesigned-flexor measurements and their actual source/pressure metadata. The 2×2 instruction was posted on the χ₄-table thread, but its target is confirmed as this figure.

## P2 — Routing geometry and supplementary candidate

Paper `fig:ExtensorRoute_redesigned`: reduce the oversized eight-pose diagram; add recognizable femur/tibia outlines; Ben permits fewer representative poses. Generator lead: `tmp/paper_bridge/figscripts/gen_ExtensorRoute_A4.m`. Select poses that explain contact/release behavior, rather than retaining eight merely to fill a grid. Maintain ≥10 pt final text, meaningful dimensions, and a clear legend. The dissertation's same-named figure uses a different composite/hardware photograph; reconcile purpose before replacing it.

Ben's conditional supplementary-material comment anchors to `fig:RouteVariants`, not the routing-geometry figure. Count the final floats against the actual target-journal rules before moving it. Current Main.tex has eleven figure floats and four table floats; no journal limit was presumed by this audit.

## P2 — Table fit, χ₄, and notation

The paper's fourteen-column `tab:FlxPinGof` leaves the page and has the wrong orientation. Reduce unnecessary repetition or split into readable tables before shrinking text. Check the final page in its intended reading orientation. Preserve values and memberships from the numerical audit.

Remove the paper's empty χ₄ column in `tab:AdoptedXi`; this model declares χ₀–χ₃. Explain why χ₃ is not fitted to pinned-flexor identification but is used in later flexor redesign, rather than reconciling the paper dash and dissertation 0.158 silently.

Ben's dissertation requests also apply to shared equations: uppercase torque M, and `\mathbf{k}` for knee-axis subscripts, distinguishing them from basis/summation indices. Apply consistently to axes, captions, equations, and generated figures. Preserve moment-arm variables.

## P2 — Dissertation mechanics schematic

At `fig:xiProjection`, split the existing compliance/wrapping composite into two figures as Ben requested and put alpha on the circle. The follow-up “Then we don't need to mention” has an incomplete instruction; inspect its saved source anchor before deleting an explanation. This request is dissertation-specific; the paper has a different mechanics-figure set.

## Return package

Use the established `Dissertation/Notes/advisor_revision_20261023/REVIEW_INBOX/<machine>/<new-batch>/` with `REVIEW.md`, bounded `patches/`, and `figures/`. Include original and replacement filenames, labels, numerical changes, and a series-level CSV for audited plots. Keep the current GoF prose edits from this audit. Use the dissertation, MATLAB, mechanics, and relevant plotting skills. No new optimization is needed for these repairs.
'''
(w/'ZCODE_REPLOT_REQUESTS.md').write_text(replots,encoding='utf-8')
summary='''# Torque-paper comparison and GoF follow-up

Created `Documentation/Reports and Papers/Knee_Torque_Test/ChatGPT_reportTorque.md`. It inventories all 31 initial active paper review items, all 13 dissertation torque review items, nine current Ben paper comment threads, and the complete 272-thread resolved archive (74 substantive Ben messages). It maps comments to snapshot source lines/labels and records the actionable differences, including six versus nine cam-extensor measurements and conflicting biomimetic-extensor GoF values.

Ben confirmed `fig:FlxGrp` should use a 2×2 grid: 10 mm, 20 mm, 20 mm, and legend tile. Detailed local ZCode instructions are in `Knee_Torque_Test/review_audit_20261026/ZCODE_REPLOT_REQUESTS.md`; no message was sent and no plot was regenerated.

Ben additionally authorized GoF prose improvements in both documents and clarified that the metrics have a combined effect rather than a one-to-one mapping to error types. Those two bounded passages were revised live, with exact whole-file readback checks and a saved diff. No scientific values, figures, or review dispositions were changed. Paper live compile: 0 errors, 20 remaining duplicate-label/font warnings. Dissertation live compile status is recorded with the audit evidence.

---

'''
r=root/'CHATGPT_REPORT.md'; r.write_bytes(summary.encode('utf-8')+r.read_bytes())
print(json.dumps({'report':str(folder/'ChatGPT_reportTorque.md'),'paper_active_items':len(entries),'dissertation_active_items':len(dentries),'historical_threads':len(archive),'ben_historical_messages':len(ben_history),'report_lines':len(rows),'figures':len([x for x in manifest['paper'] if x['type'].startswith('figure')]),'tables':len([x for x in manifest['paper'] if not x['type'].startswith('figure')])},indent=2))
