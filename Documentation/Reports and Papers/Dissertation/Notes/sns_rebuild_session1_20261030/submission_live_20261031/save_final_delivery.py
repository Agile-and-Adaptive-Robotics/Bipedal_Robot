from pathlib import Path
import zipfile, hashlib, json, sys
import pymupdf as fitz
sys.stdout.reconfigure(encoding='utf-8')
repo=Path(r'D:\GitHub\Bipedal_Robot')
dis=repo/'Documentation/Reports and Papers/Dissertation'
notes=dis/'Notes/sns_rebuild_session1_20261030/submission_live_20261031'
archive=dis/'Archive/session1_submission_20261031'
archive.mkdir(parents=True, exist_ok=True)
source=Path(r'C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.zip')
pdf=notes/'Bolen_Dissertation_SESSION1.pdf'
(archive/'Overleaf_source.zip').write_bytes(source.read_bytes())
(archive/'Bolen_Dissertation_SESSION1.pdf').write_bytes(pdf.read_bytes())
files={}
with zipfile.ZipFile(source) as z:
 for basename in ['60-sns.tex','95-AppendixD.tex']:
  matches=[n for n in z.namelist() if n.endswith('/'+basename) or n==basename]
  assert len(matches)==1,matches
  b=z.read(matches[0]); target=archive/'source'/basename;target.parent.mkdir(exist_ok=True);target.write_bytes(b)
  files[basename]=b.decode('utf-8')
 c=files['60-sns.tex'];a=files['95-AppendixD.tex']
 for s in ['LI_DENG_FALL_FRAME.png','W2LNEW_TREADMILL_STILL.png','fig:animatlab_shared_rg','tab:sns_runtime','DENG_RAT_NATIVE_SUBNETWORKS.pdf','W2L108_MAIN_SUBNETWORKS.pdf']:
  assert s in c,s
 assert '\\ContinuedFloat' not in c
 assert 'synthetic hysteresis term' in c
 assert 'Additional execution checks.' in a
print('Verified source archive matches final live edits.')
d=fitz.open(pdf)
assert len(d)==310
assert not [i+1 for i,p in enumerate(d) if '??' in p.get_text()]
assert 'Figure 6.1:' in d[152].get_text()
assert 'Figure 6.2:' in d[153].get_text()
assert '122.120' in d[180].get_text() and '17.75' in d[180].get_text()
d[152].get_pixmap(matrix=fitz.Matrix(1.7,1.7),clip=fitz.Rect(55,150,565,325)).save(notes/'Figure_6_1_preview.png')
evidence=[]
for name,terms in [('60-sns.tex',['fig:animatlab_model_architecture','fig:animatlab_shared_rg','Headless execution','fig:w2l108_subnetworks','fig:rat_native_subnetworks','tab:sns_runtime','synthetic hysteresis term','compliant-tendon prototype']),('95-AppendixD.tex',['fig:w2l_crossengine','fig:gait2392_ik_overlay','Additional execution checks.'])]:
 lines=files[name].splitlines()
 for term in terms:
  evidence.append(f'- `{name}:{next(i+1 for i,x in enumerate(lines) if term in x)}` — {term}.')
record=f'''# Session 1: final live SNS revision

Machine: **easteregg2**. Ben explicitly authorized live Overleaf edits during the final submission pass. This supersedes the earlier editor-only manuscript handoff instruction for this pass.

## Delivery

- Live project: https://www.overleaf.com/project/644dc82bcd6a5481e3c3a8bb
- Final PDF: `{pdf}`
- Complete downloaded live source archive: `{archive/'Overleaf_source.zip'}`
- Exact edited chapter copies: `{archive/'source'}`. These are archival snapshots; the local working ProofFinal was not overwritten with an older source.
- Browser-download PDF: `C:\\Users\\Ben Bolen\\Documents\\ME603 Dissertation\\Defense\\Bolen_Dissertation (2).pdf`.
- Final downloaded PDF: **310 pages**. Live compiler: **Errors 0, Warnings 0**. PDF text has no unresolved `??` references. No edits to the published muscle chapter or torque chapter.

## Requested changes completed

- Figure 6.1 now pairs the original walker (A) with Ben's supplied ground-test fall frame (B), in a compact float.
- The old shared-RG controller is standalone **Figure 6.2**, on its own sideways page (PDF page 154, printed page 133). It is no longer a continued panel.
- The newer supported treadmill still is separate Figure 6.4; 9 left / 8 right qualifying swings over 10 s, with 40.5% mean net body-weight assistance.
- Native W2L108 subnetwork schematic includes contact and bilateral joint flexor/extensor type-II feedback (Figure 6.6).
- Native public rat subnetwork figure includes the biarticular subsystem. Caption identifies the public precursor rather than claiming the final 2023 model was reproduced.
- Added actual Gait2392 geometry with muscle routes. Preserved the CAD/SNS-Simscape/BPA standing example and constructive future-work wording.
- Moved the matched-input engine plot and Gait2392 IK diagnostic plot to Appendix D while retaining labels and references. Existing editor D4 replacement was preserved.
- Added concise GUI/headless/API/SDK explanation and acceleration procedure: stored topology, sparse connections, preallocated arrays, compiled Numba calculations, and reference checks.
- Table 6.4 compares complete loops: Gait2392 8 s simulation, 122.120 to 13.253 s wall time (9.21x); W2L80 4 s simulation, 17.208 to 0.970 s wall time (17.75x). Different references are named. No unmatched AnimatLab/Simscape timing was invented.
- Retained compliant-tendon and damping prototype coverage; added the separate BPA pressure-dynamics and synthetic hysteresis test, with dynamic identification framed as future research.
- Additional old execution measurements were moved to Appendix D, preserving their numerical record rather than silently deleting them.

## Source line evidence

{chr(10).join(evidence)}

## Timing provenance

- `{repo/'Code/MuJoCo_SNS/spinal/hunt_training/runs/backend_benchmarks.json'}`
- `{repo/'Code/MuJoCo_SNS/spinal/virtual_walker_20261031/runs/BENCHMARK_FINAL.json'}`
- `{repo/'Code/MuJoCo_SNS/robot_torso/bpa_neural_coordination/EXPERIMENT_RECORD.md'}`

## Saved next modeling instruction

> Keep working on Gait2392. Fix stage 3 with the same controller, then work toward less weight support in stage 4. Don't let the knees hyperextend. Use the speed-ups and keep the separate training record updated.

No new controller training was run during this final manuscript pass. The separate experiment records remain authoritative for model performance.
'''
(notes/'SESSION1_FINAL_HANDOFF.md').write_text(record,encoding='utf-8')
for target in [repo/'CHATGPT_REPORT.md',repo/'SESSION1_HANDOFF_20261031.md']:
 with target.open('a',encoding='utf-8') as f:
  f.write('\n\n## Session 1 — final live submission revision, 2026-10-31\n\nBen authorized direct live Overleaf edits. Final PDF: 310 pages; Errors 0 / Warnings 0. Figure 6.1 pairs original walker and fall; original controller is standalone sideways Figure 6.2. Added W2L108 contact/type-II and rat subnetwork plates, supported treadmill still, human geometry, headless/GUI explanation, scoped speed table, and actuator fidelity paragraph. Detailed engine and IK plots moved to Appendix D. Current complete source ZIP and PDF preserved in Dissertation/Archive/session1_submission_20261031. Full change record: '+str(notes/'SESSION1_FINAL_HANDOFF.md')+'. No new training in this final publication pass.\n')
result={'machine':'easteregg2','pdf_pages':len(d),'live_errors':0,'live_warnings':0,'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest(),'source_zip_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_archive':str(archive/'Overleaf_source.zip'),'pdf':str(pdf),'source_checks':'pass','undefined_reference_markers':0}
(notes/'FINAL_VERIFICATION.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result,indent=2))
