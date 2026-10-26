from pathlib import Path
import json,datetime
root=Path(r'D:\Github\Bipedal_Robot')
work=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026'
patches=json.loads((work/'patches.json').read_text(encoding='utf-8'))
applied=json.loads((work/'applied.json').read_text(encoding='utf-8'))
begin='<!-- LIVE_ASTRA_20261026_START -->'; end='<!-- LIVE_ASTRA_20261026_END -->'
rows=[begin,'# Live Overleaf revision — Astra, 26 October 2026','',f'Updated {datetime.datetime.now().astimezone().isoformat(timespec="seconds")}. Status: in progress. Ben explicitly authorized live edits and asked that this report be updated as work proceeds.','', 'Purpose: describe the knowledge gap and how the created tools address it; repair unclear SNS prose and align the introduction, discussion, and conclusion. The peer-reviewed isometric-force chapter is excluded.','', 'Live project: https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb','', 'Evidence and local before/after source: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/LIVE_ASTRA_20261026/`. A fresh Overleaf source archive was retrieved before substantive edits. Changes are made through the live editor in Reviewing mode; no advisor comment has been resolved and no existing tracked change has been accepted or rejected.','', '| Patch | File | Applied change |','|---|---|---|']
for i in applied: rows.append(f'| A{i:02d} | `{patches[i]["file"]}` | {patches[i]["reason"]} |')
rows+=['',f'Applied so far: {len(applied)} of {len(patches)} prepared text patches. Each replacement is checked against the exact selected live text before insertion. The four main revised chapters matched their prepared text on whole-file readback (line-ending normalization only).','', 'Additional completed actions: inserted the detailed conversion-repair section in Appendix D; uploaded the correct Connor–Ben shared-RG diagram and the reviewed contact/Ib ablation figure into `figs/60-sns/`. First complete local build: zero errors, missing figures, or undefined references. Overleaf also reported Errors 0 and Warnings 0; box-layout messages remain. Final build following the background and Rybak additions is in progress.','', 'The live comments and Ben replies visible around the edited passages were captured in `before/review-panel.txt`. Replacing a passage addresses its content; it does not declare the advisor thread resolved. The original deleted suspended-stepping qualifier remains deleted.','',end,'']
path=root/'CHATGPT_REPORT.md'; old=path.read_text(encoding='utf-8-sig')
if begin in old:
    a=old.index(begin); b=old.index(end,a)+len(end)
    old=old[:a]+old[b:].lstrip('\r\n')
    while old.startswith('---\n\n'): old=old[5:]
path.write_text('\n'.join(rows)+'\n---\n\n'+old,encoding='utf-8')
print('Updated CHATGPT_REPORT.md:',len(applied),'applied patches')
