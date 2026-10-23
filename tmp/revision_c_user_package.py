from pathlib import Path
import json,re,shutil,zipfile,hashlib
C=Path(r'\\?\D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');S=C/'review/user_sweep_20261023';U=C/'OVERLEAF_UPLOAD';B=Path((C/'review/latest_build.txt').read_text(encoding='utf-8'))
# Detect concurrent Overleaf source edits before replacing anything.
for p in (S/'overleaf_before').glob('*.tex'):
 old=(S/'before/chapters'/p.name).read_text(encoding='utf-8-sig');live=p.read_text(encoding='utf-8-sig')
 assert old.rstrip()==live.rstrip(),f'Concurrent Overleaf edits: {p.name}'
for p in (C/'chapters').glob('*.tex'):shutil.copy2(p,U/'chapters'/p.name)
for chapter,names in {'15-background':['festo_bpas_inflated_deflated.jpg','steeleknee.pdf','robot_pam_routes.pdf','human_gait2392_paths.pdf','moment_arm_methods.pdf'],'30-torque':['cam_extensor_hardware.jpg','bktFrame.pdf','bpa_wrapping_loss.pdf']}.items():
 for name in names:shutil.copy2(C/'figs'/chapter/name,U/'figs'/chapter/name)
manifest={p.relative_to(U).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in U.rglob('*') if p.is_file()}
(S/'upload_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
with zipfile.ZipFile(C/'Revision_C_user_sweep_Overleaf.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(U.rglob('*')):
  if p.is_file():z.write(p,p.relative_to(U).as_posix())
shutil.copy2(B/'main.log',C/'review/build.log')
for p in B.glob('pass_*.log'):shutil.copy2(p,C/'review'/p.name)
checks=json.loads((S/'checklist.json').read_text(encoding='utf-8'))
for k,v in checks.items():
 if k not in ['14','21','22','25']:v['status']='DONE'
checks['06']['evidence']+=' Full 5-pass build passed; figures3.2–3.4 rendered and visually checked within page bounds.'
checks['07']['evidence']+=' All 36 numeric results checked against original; centered readable table visually verified.'
checks['25']['evidence']='Full fresh build succeeded (all 5 passes), no undefined refs/cites, missing figures, or oversized floats. Selected pages rendered and inspected. Chapter3 source plots byte-identical; Table3.3 all36 results and Tables4.1–4.3 bodies preserved exactly. Current PDF Revision_C_user_sweep_20261023.pdf. Clean upload bundle refreshed; six live Overleaf sources checked against before snapshot (no intervening edits). Upload pending; ZCode items14/21 remain open.'
(S/'checklist.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False),encoding='utf-8')
root=C.parents[5]/'CHATGPT_REPORT.md'
# Resolve repo explicitly rather than depending on folder depth.
root=Path(r'D:\Github\Bipedal_Robot\CHATGPT_REPORT.md');s=root.read_text(encoding='utf-8')
for k,v in checks.items():
 line=f'- [{"x" if v["status"]=="DONE" else " "}] U{k} [{v["status"]}] {v["request"]} Evidence/status: {v["evidence"]}'
 s=re.sub(r'^- \[[ x]\] U'+k+r' \[[^\n]*',lambda m:line,s,flags=re.M)
root.write_text(s,encoding='utf-8')
print('Bundle refreshed; live sources match saved baseline; QA recorded')


