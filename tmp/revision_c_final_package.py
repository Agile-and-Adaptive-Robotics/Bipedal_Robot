from pathlib import Path
import re,json,hashlib,shutil,zipfile,fitz
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';C=D/'Notes/advisor_revision_20261022/Revision_C';S=C/'review/final_sweep_20261023';B=Path((C/'review/latest_build.txt').read_text())
log=(B/'main.log').read_text(errors='replace')
for pattern in [r'^!',r'There were undefined',r'LaTeX Warning:.*undefined',r'multiply defined',r'Float too large',r'Label\(s\) may have changed']:
    assert not re.search(pattern,log,re.M),pattern
assert len(list(B.glob('pass_*.log')))==5
pdf=C/'Revision_C.pdf'
try:shutil.copy2(B/'main.pdf',pdf)
except PermissionError:
    pdf=C/'Revision_C_final_20261023.pdf';shutil.copy2(B/'main.pdf',pdf)
for name in ['main.log','main.toc']:
    shutil.copy2(B/name,C/'review'/('build.log' if name.endswith('.log') else 'build.toc'))
for p in B.glob('pass_*.log'):shutil.copy2(p,C/'review'/p.name)
def label(file,key):
    text=(B/'chapters'/file).read_text(errors='replace')
    m=re.search(r'\\newlabel\{'+re.escape(key)+r'\}\{\{([^}]+)\}\{([^}]+)\}',text)
    assert m,key;return m.group(1),int(m.group(2))
pair=[label('15-background.aux',key) for key in ['fig:rybak_cpg_circuit','fig:shevtsova_spinal_circuit']]
assert pair[0][0]=='2.8' and pair[1][0]=='2.9' and pair[0][1]==pair[1][1],pair
sadb=label('50-database.aux','ch:database')[1];sns=label('60-sns.aux','ch:sns')[1];overall=label('70-overall.aux','ch:overall')[1]
assert r'\label{fig:sadb_clusters}' in (C/'chapters/50-database.tex').read_text()
assert (C/'chapters/30-torque.tex').read_text().count(r'\begin{table}')==2
doc=fitz.open(pdf);offset=next(page for level,title,page in doc.get_toc() if level==1 and title=='Introduction')-1
out=S/'rendered';out.mkdir(exist_ok=True)
for fname,key in [('15-background.aux','fig:rybak_cpg_circuit'),('94-AppendixC.aux','tab:FlxPinPreGof'),('94-AppendixC.aux','tab:FlxPinGof'),('94-AppendixC.aux','tab:ExtPinGof'),('94-AppendixC.aux','tab:TorqueFit')]:
    n=label(fname,key)[1]+offset-1;doc[n].get_pixmap(matrix=fitz.Matrix(150/72,150/72)).save(out/f'page_{n+1:03}.png')
# Verify that relocation preserves every numeric data entry in the four table bodies.
old=(S/'before/30-torque.tex').read_text();new=(C/'chapters/94-AppendixC.tex').read_text()
def body(s,key):
    for m in re.finditer(r'\\begin\{(table|sidewaystable)\}.*?\\end\{\1\}',s,re.S):
        if r'\label{'+key+'}' in m.group():
            b=m.group().split(r'\midrule',1)[1].split(r'\bottomrule',1)[0]
            return re.findall(r'[-+]?\d+(?:\.\d+)?',b)
    raise ValueError(key)
for key in ['tab:FlxPinPreGof','tab:FlxPinGof','tab:ExtPinGof','tab:TorqueFit']:assert body(old,key)==body(new,key),key
integrity=json.loads((C/'review/source_integrity.json').read_text());assert not integrity['source_files_changed'] and not integrity['missing_labels']
assert len(json.loads((S/'comment_dispositions.json').read_text()))==124
# Prepare only compilable source, support files, and actually referenced graphics.
upload=C/'OVERLEAF_UPLOAD'
if upload.exists():
    for item in upload.rglob('*'):
        if item.is_file():
            original=C/item.relative_to(upload)
            assert original.exists() and item.read_bytes()==original.read_bytes(),f'Preserve modified upload file: {item}'
upload.mkdir(exist_ok=True)
for p in C.iterdir():
    if p.is_file() and p.suffix.lower() in ['.tex','.bib','.bst','.cls','.sty']:shutil.copy2(p,upload/p.name)
(upload/'chapters').mkdir(exist_ok=True)
for p in (C/'chapters').glob('*.tex'):shutil.copy2(p,upload/'chapters'/p.name)
assets=set()
for tex in [C/'main.tex',*sorted((C/'chapters').glob('*.tex'))]:
    text='\n'.join(re.split(r'(?<!\\)%',line,maxsplit=1)[0] for line in tex.read_text(encoding='utf-8-sig').splitlines());g=re.search(r'\\graphicspath\{((?:\{[^{}]*\})+)\}',text)
    paths=[C]
    if g:paths += [C/d for d in re.findall(r'\{([^{}]*)\}',g.group(1))]
    for name in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}',text):
        assert '\\' not in name,(tex,name)
        candidates=[]
        for base in paths:
            q=base/name
            if q.suffix:candidates.append(q)
            else:candidates += [q.with_suffix(ext) for ext in ['.pdf','.png','.jpg','.jpeg','.eps']]
        source=next((q for q in candidates if q.exists()),None);assert source,(tex,name)
        rel=source.resolve().relative_to(C.resolve());assets.add(rel)
        if source.suffix.lower()=='.eps':
            cached=source.with_name(source.stem+'-eps-converted-to.pdf')
            if cached.exists():assets.add(cached.resolve().relative_to(C.resolve()))
for rel in sorted(assets):
    target=upload/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(C/rel,target)
manifest={str(p.relative_to(upload)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in upload.rglob('*') if p.is_file()}
assert not any('report' in name.lower() or '/review/' in name.lower() for name in manifest)
assert hashlib.sha256((upload/'thesis.bib').read_bytes()).hexdigest()==hashlib.sha256((C/'thesis.bib').read_bytes()).hexdigest()
(S/'upload_manifest.json').write_text(json.dumps(manifest,indent=2))
zpath=C/'Revision_C_Overleaf.zip'
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(upload.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(upload).as_posix())
stats=json.loads((C/'review/word_counts.json').read_text())
for file in ['REPORT.md','README.md']:
    dest=S/('previous_'+file)
    if not dest.exists():shutil.copy2(C/file,dest)
report=f'''# Revision C — final local sweep, 2026-10-23

**Use `OVERLEAF_UPLOAD`, not this whole Revision_C folder, for Overleaf.** It contains only the compilable source and required graphics. `Revision_C_Overleaf.zip` contains the same clean project with `main.tex` at its root. Reports, review records, build logs, and the compiled review PDF are excluded.

## Changes made

- **SNS gap:** the opening credits Ben's Li/Deng construction and interconnected left/right RG extension, then identifies the five missing combinations in his advisor reply: human model/CPG, BPA-routed robot/CPG, Li/Deng with foot contact, graphical neural construction with maintained mechanics, and SolidWorks/BPA/neural integration. It explains the knowledge needed to compare circuit behavior across biological and artificial mechanics, names the tools built, and separates demonstrated interfaces from incomplete integrations.
- **SADb gap:** immediately after introducing the resource, the chapter states what evidence was missing for feedback design and what the database supplies. Curation notes and pathway records are defined, relevance criteria are explicit, and Ben's decision authority is clear.
- **Current app:** the embedded export was read directly. It contains **1,393 records: 943 papers and 450 non-paper records**, with **436 archived and 957 retained**. Chapter 5 and Appendix E now reflect the completed full-base export. The 561 retained-paper coverage counts (523 notes, 426 preparation tags, 86 afferent classifications, and 258 robot/simulation notes) were checked against the actual embedded data. The monthly scraper remains under development in its separate session.
- **Bubble screenshot:** the paper-only citation-neighborhood screenshot is now in the SADb main chapter, Figure {label('50-database.aux','fig:sadb_clusters')[0]}. The paper-table screenshot moves to Appendix E. The graph remains paper-only even though the table now includes all record types.
- **Chapter 4 tables:** only **two compact tables** remain in the main chapter: a four-row mean-error summary and five-row adopted-correction table. The four full fit tables move to Appendix C with their original labels and all numeric body entries unchanged. Captions and duplicate bounds detail are reduced in the main text.
- **Completed physical tests:** the committed MATLAB exports are included as Figures {label('30-torque.aux','fig:extensor_physical_final')[0]} and {label('30-torque.aux','fig:flexor_physical_final')[0]}. The scripts contain nine extensor and fifteen flexor holds, spanning -119 to +5 degrees. The figures and text distinguish measured torques, predictions at measured conditions, the 620 kPa design curves, and human targets. No Excel sanity-check torque values were used. The data support torque capacity while also showing discrepancies in pointwise prediction; no new fit scores were invented.
- **Caption/layout:** the cat caption is one short sentence retaining source and restoration credit. **Figures 2.8 and 2.9 share printed page {pair[0][1]}**, with separate captions and labels. The motor/stator centers in the existing replacement schematic were visually verified concentric. The confusing graphical abstract and deleted circuit catalogue remain absent from the Background.
- **Background/Discussion:** moment arms now have a correct signed fixed-axis definition and are distinguished from OpenSim length sensitivities. The general computational-neuroscience section no longer digresses into the dissertation's CPG. The monthly scraper and active bipulley campaign are not presented as completed results.

## Advisor-comment audit

The expanded live Overleaf Overview matches the saved ledger: **73 advisor comment threads + 51 tracked edits = 124 items** (Introduction 34, Background 49, SADb 26, SNS 15). Every item has a local disposition and source-line evidence in `COMMENT_DISPOSITIONS.md`. Ben's replies remain separately identified; his four appendix/main items were checked too. CHATGPT_REPORT.md now records this sweep. No new open advisor items were detected.

These are local editorial dispositions, not Overleaf resolutions or advisor approval. No comment was resolved, no tracked change was accepted/rejected, and no source was uploaded. The historical resolved archive remains preserved. The library comparison requested in SNS-15 is functional; an unperformed performance benchmark is not claimed. Incomplete SNS integrations remain explicit rather than being declared successful.

## Build and conservation checks

- Full fresh MiKTeX build: pdflatex, BibTeX, then three pdflatex passes; all five commands succeeded. No undefined references/citations, duplicate labels, missing figures, changed-label rerun warning, or oversized floats in the final build.
- Review PDF: **{pdf.name}, {len(doc)} physical pages**. SADb occupies printed pages **{sadb}–{sns-1} ({sns-sadb} pages)**; SNS **{sns}–{overall-1} ({overall-sns} pages)**. The additional physical-test plots and full statistical tables remain in the complete document.
- The four moved fit-table bodies were compared numerically against the pre-sweep source: every entry matches. The mean summary is transcribed from the existing mean rows, not refitted.
- Protected force chapter and figures remain byte-identical. All {integrity['protected_figure_hashes_verified']} protected-manifest entries passed; bibliography bytes and the original ProofFinal/proposal sources are unchanged by this pass.
- Selected changed pages, the paired diagrams, physical-test plots, and relocated tables were rendered at 150 dpi and inspected. The final log is in `review/build.log`; source-line and figure/number changes are recorded in the review ledgers.
- Inherited wide URLs, equations, and old appendix tables still produce overfull-box warnings. They are not declared repaired. The new compact Chapter 4 tables and paired figure have been visually checked.

## Upload

1. Keep a current Overleaf project backup using **File → Download → source ZIP** before replacing files in an actively reviewed project.
2. Open **OVERLEAF_UPLOAD** locally, select **its contents**, and drag them into the Overleaf project root. Preserve the `chapters/` and `figs/` folders and allow replacement of matching files. Do not upload the enclosing Revision_C folder or its reports.
3. Keep **main.tex** selected as the main document and use pdfLaTeX; recompile. The new Appendix D and E includes are already present in main.tex.

Overleaf file replacement may disturb source-linked review anchors; the complete current comments and local dispositions are retained outside the upload bundle. No upload or review-state change was performed by this pass.

## Files for Ben/ZCode

- `COMMENT_DISPOSITIONS.md`: all 124 current advisor IDs, request summaries, dispositions, and source-line links.
- `review/final_sweep_20261023/final_sweep.diff`: this last pass only, including moves and deletions.
- `review/changes.diff`: cumulative Revision C review diff.
- `review/final_sweep_20261023/table_moves.md`: full-table conservation and compact-summary rationale.
- `review/final_sweep_20261023/physical_figure_sources.json`: exact source paths and hashes for the script-generated plots.
- `review/final_sweep_20261023/app_payload_check.json`: actual app counts and retained-paper coverage.
- `review/final_sweep_20261023/upload_manifest.json`: {len(manifest)} upload files, including {len(assets)} referenced graphic assets/cache files, with hashes.

The title-page submission date remains the original July 2026 metadata; it was outside this content pass. No commit, push, optimization run, scraper run, data-reduction rerun, or bibliography rewrite was performed.
'''
(C/'REPORT.md').write_text(report,encoding='utf-8')
(C/'README.md').write_text(f'''# Revision C

Read **{pdf.name}**. **REPORT.md** summarizes the final sweep. **COMMENT_DISPOSITIONS.md** covers all 124 current advisor items.

For Overleaf, open **OVERLEAF_UPLOAD**, select **its contents**, and drag them into the project root while preserving **chapters/** and **figs/**. Replace matching files, keep **main.tex** as the main document, and recompile using pdfLaTeX. Keep a current source ZIP backup first because file replacement can disturb review anchors.

**Do not drag the entire Revision_C folder into Overleaf.** The clean upload folder and **Revision_C_Overleaf.zip** contain no reports, comment inventories, review logs, or compiled review PDF. The ZIP is also suitable for creating a separate Overleaf project.

The compiled copy passes the full MiKTeX/BibTeX sequence. Original ProofFinal and Overleaf were not edited. The new physical-test plots come from the committed analysis-script exports; Excel torque cells remain sanity checks. The full SADb record export is included in the prose; monthly scraper development and the active bipulley run remain separate.
''',encoding='utf-8')
(C/'review/render_manifest.json').write_text(json.dumps({'pdf':str(pdf),'pdf_pages':len(doc),'current_render_directory':str(S/'rendered'),'figures_2_8_and_2_9_printed_page':pair[0][1]},indent=2))
print(json.dumps({'pdf':str(pdf),'pages':len(doc),'SADb_pages':[sadb,sns-1],'SNS_pages':[sns,overall-1],'paired_figures_page':pair[0][1],'upload_files':len(manifest),'graphic_files':len(assets),'zip_bytes':zpath.stat().st_size},indent=2))
