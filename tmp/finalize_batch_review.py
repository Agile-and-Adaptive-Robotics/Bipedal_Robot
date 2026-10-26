from pathlib import Path
import json,shutil
import pymupdf
root=Path(r'D:\Github\Bipedal_Robot')
base=root/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023'
out=base/'REVIEW_INBOX/REVIEWED_20261026'
for name in ['OurWalker.png','connor_ben_shared_rg_original.png']:
    source=base/'Revision_E/figs/60-sns'/name
    if not source.exists(): source=base/'Revision_E/figs/60-sns/AnimatLab_reference'/name
    shutil.copy2(source,out/'figures'/name)
shutil.copy2(root/'tmp/review_batches_20261026.py',out/'evidence/review_packet_generator.py')
doc=pymupdf.open(out/'REVIEWED_FRAGMENTS.pdf')
selected=[]
for i,page in enumerate(doc):
    text=page.get_text()
    if ('Tool developed' in text or 'RG flexor: RG stimulus' in text or '51.8' in text):
        dest=out/f'proof_build_01/visual_page_{i+1}.png'
        page.get_pixmap(matrix=pymupdf.Matrix(1.3,1.3)).save(dest)
        selected.append({'page':i+1,'image':str(dest)})
(out/'proof_build_01/visual_pages.json').write_text(json.dumps(selected,indent=2),encoding='utf-8')
lines=['# Directive compliance and verification','','These line numbers refer to the delivered patch files, not stale Overleaf positions. Live advisor-thread mapping still requires checking the application anchors.','']
checks=[('Researchers and Gait2392 second','P01_replacement.tex','Researchers'),('Actual referenced tools table','P01_replacement.tex','Table~'),('Table Gait2392 second','P01_companion_table.tex','Gait2392/SNS'),('Correct Ib versus commissural claim','EG_controlled_tests.tex','Removing Ib'),('Mass count/height/tilt limits','EG_controlled_tests.tex','Across ten'),('Partial feedback and scored-window scope','EB_feedback_ablation.tex','At fixed'),('Original engine and original model scope','DENG_original_engine.tex','These are observations'),('Actual panel letters and preserved labels','EG_extensor_merged_figure.tex','makebox'),('Captionless figure counter','EG_physical_figures_captionless_OPTION.tex','refstepcounter'),('Two-part tackle','EG_pulley_discussion.tex','two-part reversed')]
for title,name,needle in checks:
    content=(out/'patches'/name).read_text(encoding='utf-8').splitlines()
    hits=[str(i+1) for i,l in enumerate(content) if needle in l]
    lines.append(f'- {title}: `patches/{name}`, line(s) {", ".join(hits)}.')
lines+=['','Proof build: pdflatex, BibTeX, and two further pdflatex passes all exited zero. Final log: no errors, undefined citations/references, or overfull boxes. Six references to dissertation material outside this fragment packet are explicitly mocked as X, as recorded in `proof_result.json`. This is a fragment proof, not the assembled dissertation.','',f'Proof length: {len(doc)} pages. The new Deng plot was visually inspected before and after layout repair. Figure/model simulations were not rerun.','', 'Overleaf, active chapter source, and existing bibliography were not edited. No commit, push, or outgoing ZCode message was made.']
(out/'COMPLIANCE.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
entry='''# Three-batch review: reconciled packet, 26 October 2026

Reviewed EB475WS4, easteregg2, and easteregg2_2. Delivery: `Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/REVIEW_INBOX/REVIEWED_20261026/REVIEW_DECISIONS.md`. Use its bounded patches and application map; do not apply the incoming batches verbatim. `ZCODE_REPAIR_NOTES.md` identifies held figures and exact repairs. `REVIEWED_FRAGMENTS.pdf` is the standalone proof, not a chapter replacement.

Approved tools framing retained (Researchers; Gait2392 second). Corrected stale W2L necessity/mass claims, scoped the contact/Ib ablation, repaired the original-engine Deng figure from preserved raw charts, completed a separate publisher-verified bibliography fragment, and corrected torque panel lettering/label handling. Knee CPG figure is held because the plotted interval does not show the claimed cycle. Knee-reflex and loaded-arm graphics need clipping/annotation repairs. Repeated failures and ownership are recorded in the packet's self-scoring matrix. Prior invalidated counts/no-op deletions remain invalidated.

Fragment proof compiles cleanly; external chapter references are explicitly mocked. No Overleaf edits, active chapter replacements, bibliography rewrites, model runs, commits, pushes, or messages to ZCode were performed. Application must reconcile the live text, advisor comments, and Ben's tracked deletions first.

---

'''
for path in [root/'CHATGPT_REPORT.md',root/'ZCODE_HANDOFF_20261023.md',base/'Revision_E/ZCODE_HANDOFF.md']:
    old=path.read_bytes()
    if not old.startswith(entry.encode('utf-8')): path.write_bytes(entry.encode('utf-8')+old)
print(json.dumps(selected,indent=2))
