from pathlib import Path
import json,re
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/user_sweep_20261023';p=S/'checklist.json';data=json.loads(p.read_text(encoding='utf-8'))
evidence={
'01':'30-torque.tex: quoted sentence removed; matching Appendix C selection narrative removed. Before snapshot preserves it.',
'02':'Copied original JPEGs without resampling to figs/15-background/festo_bpas_inflated_deflated.jpg and figs/30-torque/cam_extensor_hardware.jpg; original phone files untouched. Photo wired into Fig.2.3; camera EXIF rotation handled in TeX. Visual page check pending.',
'03':'steeleknee.pdf: uppercase A/B native PDF overlays; caption uppercase too.',
'04':'15-background.tex: three top-aligned minipages, letters at upper left.',
'05':'Fig.2.5 now skeleton / robot PAM routes / human Gait2392 paths; split original masters_leg_model.pdf into robot_pam_routes.pdf and human_gait2392_paths.pdf. Fig.2.6 replaced duplicate renderings with moment_arm_methods.pdf.',
'06':'20-muscles.tex: 3 figures use page floats, centered, bounded panel width/height. Source EPS/MATLAB plots unchanged. Full page verification pending.',
'07':'Table model_compare unpivoted into 12 readable rows (model per row), centered, small rather than tiny, no negative adjustwidth. Numeric conservation/page check pending.',
'08':'figs/15-background/moment_arm_methods.pdf: cross-product torque and dot-product moment-arm projection; Fig.2.6 label preserved and cross-referenced near Eq.4.2. Fixed prior self-orthogonal definitions of p_proj versus p_f in nearby prose; equation unchanged. Diagram visually checked.',
'09':'bktFrame.pdf: native text A–D replaced with bold labels; CAD/image content preserved.',
'10':'No existing BPA-wrapping schematic was located in searched dissertation/paper assets (bracket-bending assets are a different effect). Added native vector bpa_wrapping_loss.pdf as panel B of existing Fig.4.5, referenced at Eq.4.24, preserving later figure numbering. Diagram states fitted loss is not arc-minus-chord length.',
'11':'Tables4.1–4.3 captions shortened; point-set restrictions and counts moved before Table4.1; fold membership already in rows; excluded flexor test retained in prose. Numbers unchanged.',
'12':'Table4.4 header meters; +3.94mm -> +0.00394m, -6.38mm -> -0.00638m. Other parameters unchanged.',
'13':'Added paragraph before route redesign: anatomical bounds, relaxation for RoM, subsequent lower-force hardware tradeoff; cam paragraph simplified. No new bounds or optimizer outcomes invented.',
'14':'ZCODE: user says Figure4.11 panel C is being changed by ZCode. No concurrent plot editing performed.',
'15':'cam_extensor_hardware.jpg copied byte-identically and included below cam route panels in Fig.4.12. Needs full page visual check.',
'16':'50-database.tex: structured account -> structured database.',
'17':'50-database.tex: Zotero libraries of the AARL and myself, and papers already curated in Airtable.',
'18':'Read SfN24_Bolen.pdf, extracted text and visually inspected poster. Paragraph names ResearchRabbit, Litmaps, txyz, SciSpace, MirrorThink; categorization: none good, txyz adequate, others insufficient. Limited benefits and substantial correction stated. Existing bolen_using_2024 citation retained.',
'19':'60-sns.tex: initial Connor/Ben model ONE RG both sides; air stepping at NIH2021; distinguishes Ben feedback/contact additions and later bilateral interconnected RG model.',
'20':'60-sns.tex: AnimatLab unsupported, older physics engine, slow in this workflow.',
'21':'PARTIAL / ZCODE: split verification/limit columns; removed known-false welded-knee limitation and unsupported incomplete-integration account; row names distinguish authored models/circuit implementations from external packages. ZCode must verify actual SolidWorks integration, bridge attribution/direction, controller-port/demos/curriculum scope, and all table capability/limit claims. Do NOT revive welded-joint claim. No answers guessed.',
'22':'60-sns.tex: unsupported Hunt/Young Simulink-library comparison removed. IMPORTANT: earlier Overleaf reply claiming a comparison still needs correction; do not leave it uncorrected. Actual references may be legitimate papers but do not establish a shared Simulink library.',
'23':'Defined Walker Two-Layer (W2L) at first mention; grounded in spinal/w2l_cpg/README.md / Walker_2_Layer_CPG_BilateralRG model naming.',
'24':'60-sns.tex gap now explicitly includes many degrees of freedom and 92 musculotendon actuators in Gait2392. Does not claim no such models exist anywhere.',
'25':'IN PROGRESS. Fresh build running from tmp/revision_c_build.py; session 27809. Passes1–2 succeeded at last check. Need final build, rendered changed-page checks, numeric conservation for reformatted Table3.3, updated source integrity allowing Ben-authorized Ch3 layout edits, diff/line evidence, refreshed OVERLEAF_UPLOAD/ZIP, upload finished changes and correct misleading prior Hunt/Young reply. NOTHING FROM THIS NEW SWEEP UPLOADED YET.'}
for k,v in data.items():
 v['status']='ZCODE' if k=='14' else 'PARTIAL / ZCODE' if k=='21' else 'IN PROGRESS' if k=='25' else 'EDITED / QA PENDING' if k in ['02','03','04','05','06','07','09','10','15','22'] else 'DONE'
 v['evidence']=evidence[k]
p.write_text(json.dumps(data,indent=2,ensure_ascii=False),encoding='utf-8')
root=R/'CHATGPT_REPORT.md';s=root.read_text(encoding='utf-8')
for k,v in data.items():
 pattern=r'^- \[[ x]\] U'+k+r' \[[^\n]*'
 line=f'- [{"x" if v["status"]=="DONE" else " "}] U{k} [{v["status"]}] {v["request"]} Evidence/status: {v["evidence"]}'
 s=re.sub(pattern,lambda m:line,s,flags=re.M)
heading='### Handoff checkpoint: user correction sweep'
note='''\n\n### Handoff checkpoint: user correction sweep

All edits are in Revision_C, NOT ProofFinal. Nothing from this new sweep has been uploaded to Overleaf yet. Do not run older one-shot revision scripts or packaging scripts: they would overwrite current work or assert the old table/protection state. Current sources are UTF-8; read/write explicitly with UTF-8. Bibliography remains untouched. Source MATLAB/experimental plots remain untouched; Chapter3 layout and table reformatting were explicitly authorized by Ben in this request.

Before snapshot: `Revision_C/review/user_sweep_20261023/before/` (all chapters and figs). Checklist: same folder `checklist.json`; evidence includes poster text/render, asset renders, and original copies. New scripts in root tmp: `revision_c_user_edit.py` (ONE-SHOT, already applied), `revision_c_user_assets.py` (vector PDF authoring, do not rerun on modified assets without restoring snapshot), `revision_c_user_layout2.py` (already applied), `revision_c_user_progress.py` (this checkpoint). Generating vector diagram code is in the assets script. Existing full build helper: `tmp/revision_c_build.py`.

ZCode priorities: Figure4.11 panel C, Table6.1 factual verification and attribution/direction questions, final layout QA, and correcting the prior Overleaf Hunt/Young reply. Do not invent a benchmark or claim Hunt/Young shared a neuromechanical Simulink library. Ben's original technical questions remain in U21 above. Fig4.5 now has wrapping schematic panel B to avoid renumbering the downstream torque figures.
'''
if heading not in s:s+=note
root.write_text(s,encoding='utf-8');(S/'HANDOFF.md').write_text(note.strip()+'\n',encoding='utf-8')
print('CHATGPT_REPORT.md updated: all 25 items have explicit progress and handoff evidence')
