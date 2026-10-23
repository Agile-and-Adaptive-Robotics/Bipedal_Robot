from pathlib import Path
import json,re
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/user_sweep_20261023'
p=S/'checklist.json';checks=json.loads(p.read_text(encoding='utf-8'))
checks['22']['status']='DONE';checks['22']['evidence']='Unsupported comparison removed from manuscript. Earlier Overleaf reply edited in place: it no longer claims a Hunt/Young library comparison and leaves comparison with the actual libraries pending verification.'
checks['25']['status']='DONE';checks['25']['evidence']='Fresh 5-pass local build passed. All 36 Table3.3 results and all Table4.1–4.3 data bodies conserved; Chapter3 source plots byte-identical. Changed pages rendered and visually inspected. Uploaded six changed chapters and eight changed/new figure assets individually to Overleaf. Live compilation: Errors0, Warnings0; existing typesetting warnings remain. Live photo and Fig3.3 page inspected; plot and caption fit. Prior misleading reply corrected in place. Reports not uploaded. ZCode-owned U14/U21 remain outstanding.'
p.write_text(json.dumps(checks,indent=2,ensure_ascii=False),encoding='utf-8')
root=R/'CHATGPT_REPORT.md';s=root.read_text(encoding='utf-8')
for k,v in checks.items():
 line=f'- [{"x" if v["status"]=="DONE" else " "}] U{k} [{v["status"]}] {v["request"]} Evidence/status: {v["evidence"]}'
 s=re.sub(r'^- \[[ x]\] U'+k+r' \[[^\n]*',lambda m:line,s,flags=re.M)
s+='''\n\n### User correction sweep: final uploaded status

The pre-upload checkpoint above is superseded by this status. All six revised chapter sources and eight figure assets are now on Overleaf. Compile panel shows Errors 0 and Warnings 0; inherited typesetting warnings remain. The live BPA photograph and resized Figure 3.3 were visually checked. The earlier Hunt/Young advisor reply was corrected in place. No new comments were resolved in this sweep; substantive threads remain open.

Completed: U01–U13, U15–U20, U22–U25. ZCode-owned: U14 (Figure4.11 panel C) and factual verification under U21 (Table6.1 attribution, MuJoCo–Simulink bridge direction, actual imported mechanism/actuator integration, and scope of demos/curriculum). U21 layout and the clearly false welded-knee statement were corrected; do not infer that the remaining technical questions have been answered. No performance benchmark was invented.

Latest local review PDF: `Revision_C/Revision_C_user_sweep_20261023.pdf`. Clean upload folder: `Revision_C/OVERLEAF_UPLOAD`. Latest ZIP: `Revision_C/Revision_C_user_sweep_Overleaf.zip`. Evidence, before snapshots, complete diff, checklist, numeric conservation, and screenshots: `Revision_C/review/user_sweep_20261023/`. Original phone photos and MATLAB plot files are preserved. Photo copies have descriptive names. No commit/push or ProofFinal edits.
'''
root.write_text(s,encoding='utf-8')
patterns={'01':('30-torque.tex','What the correction terms'),'02':('15-background.tex','festo_bpas_inflated_deflated'),'03':('15-background.tex','Steele biomimetic knee mechanism'),'04':('15-background.tex','fig:aarl_opensim_comparison'),'05':('15-background.tex','robot_pam_routes'),'06':('20-muscles.tex','height=0.29'),'07':('20-muscles.tex','BPA force-model comparison'),'08':('15-background.tex','moment_arm_methods'),'09':('30-torque.tex','bktFrame.pdf'),'10':('30-torque.tex','bpa_wrapping_loss'),'11':('30-torque.tex','Scoring-point restrictions'),'12':('30-torque.tex','0.00394'),'13':('30-torque.tex','initial route searches bounded'),'15':('30-torque.tex','cam_extensor_hardware'),'16':('50-database.tex','structured database'),'17':('50-database.tex','Zotero libraries'),'18':('50-database.tex','ResearchRabbit'),'19':('60-sns.tex','one rhythm generator'),'20':('60-sns.tex','unsupported, relies'),'21':('60-sns.tex','Verification} &'),'22':('60-sns.tex','Verification uses selected circuits'),'23':('60-sns.tex','Walker Two-Layer'),'24':('60-sns.tex','many joint degrees')}
rows=[]
for k,v in checks.items():
 ev=v['evidence']
 if k in patterns:
  f,term=patterns[k];lines=(C/'chapters'/f).read_text(encoding='utf-8').splitlines();n=next((i for i,l in enumerate(lines,1) if term in l),None);assert n,(k,term)
  ev=f'chapters/{f}:{n}. '+ev
 rows.append(f'| U{k} | {v["status"]} | {ev} |')
(S/'COMPLIANCE.md').write_text('# Per-directive completion\n\n| Item | Status | Evidence |\n|---|---|---|\n'+'\n'.join(rows)+'\n',encoding='utf-8')
handoff='''# Current handoff after upload

The completed local sweep is uploaded to the existing dissertation project. Overleaf: Errors0, Warnings0; inherited typesetting warnings remain. Use current Revision_C source and figures as the baseline, not the earlier uploaded-table PDF. Before snapshots, diff, and numerical checks are in this folder. All 25 directives are in checklist.json and the root CHATGPT_REPORT.md with statuses.

ZCode remaining work:
1. U14: Figure4.11 panel C (Ben says ZCode is doing the plot). Do not overwrite the other current torque edits when integrating it.
2. U21: verify Table6.1 and related prose against actual code/results. Identify authored models/interfaces versus external packages; attribute the external MuJoCo–Simulink bridge; establish one- versus bidirectional operation; distinguish circuit ports, demo tests, and curriculum; verify actual SolidWorks/actuator integration. Ben says welded-knee issue was resolved. Do not resurrect it. Verification/limit columns are already separated.
3. The specific libraries referenced by the advisor remain unidentified; do not invent a Hunt/Young Simulink library or benchmark. The false manuscript comparison and old advisor reply were corrected.

Only those factual checks and Figure4.11 panel C remain from this request. Source MATLAB plots were not edited. Ben explicitly authorized Chapter3 layout/table changes. Original photographs remain under Figures/components; renamed copies are in Revision_C/figs. No changes to ProofFinal, no commit/push. Do not run old one-shot revision scripts over these files. Use UTF-8 explicitly. For deep paths, Python writes may require the Windows extended-length prefix.
'''
(S/'HANDOFF.md').write_text(handoff,encoding='utf-8');(S/'REPORT.md').write_text(handoff+'\nSee COMPLIANCE.md for source-line evidence and verification.json for build/numeric conservation.\n',encoding='utf-8')
for name in ['REPORT.md','README.md']:
 p=C/name;s=p.read_text(encoding='utf-8');note='> Latest user correction sweep is uploaded. Use Revision_C_user_sweep_20261023.pdf locally; Overleaf compiled with Errors0 and Warnings0. See review/user_sweep_20261023/REPORT.md and COMPLIANCE.md. Only ZCode-owned Figure4.11 panel C and Table6.1 factual checks remain. Earlier pass descriptions below are historical.\n\n'
 p.write_text(note+s,encoding='utf-8')
print('Updated root report and final handoff: 23 completed, 2 ZCode-owned items')
