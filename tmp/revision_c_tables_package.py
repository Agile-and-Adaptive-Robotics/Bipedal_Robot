from pathlib import Path
import re,json,shutil,hashlib,zipfile,fitz
C=Path(r'D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022\Revision_C');S=C/'review/table_restoration_20261023';B=Path((C/'review/latest_build.txt').read_text())
log=(B/'main.log').read_text(errors='replace')
for pattern in [r'^!',r'There were undefined',r'LaTeX Warning:.*undefined',r'multiply defined',r'Float too large',r'Label\(s\) may have changed']:assert not re.search(pattern,log,re.M),pattern
pdf=C/'Revision_C_tables_20261023.pdf';shutil.copy2(B/'main.pdf',pdf)
doc=fitz.open(pdf);aux=(B/'chapters/30-torque.aux').read_text();offset=next(p for l,t,p in doc.get_toc() if l==1 and t=='Introduction')-1
pages={}
for key in ['tab:FlxPinPreGof','tab:FlxPinGof','tab:ExtPinGof','tab:AdoptedXi']:
    m=re.search(r'\\newlabel\{'+re.escape(key)+r'\}\{\{([^}]+)\}\{([^}]+)\}',aux);assert m,key
    n=int(m.group(2))+offset-1;doc[n].get_pixmap(matrix=fitz.Matrix(1.5,1.5)).save(S/(key.split(':')[1]+'.png'));pages[key]={'number':m.group(1),'page':int(m.group(2))}
shutil.copy2(B/'main.log',C/'review/build.log')
for f in B.glob('pass_*.log'):shutil.copy2(f,C/'review'/f.name)
U=C/'OVERLEAF_UPLOAD';manifest_path=C/'review/final_sweep_20261023/upload_manifest.json';old=json.loads(manifest_path.read_text())
for rel,h in old.items():assert hashlib.sha256((U/rel).read_bytes()).hexdigest()==h,rel
for name in ['30-torque.tex','94-AppendixC.tex']:shutil.copy2(C/'chapters'/name,U/'chapters'/name)
manifest={p.relative_to(U).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in U.rglob('*') if p.is_file()};manifest_path.write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(C/'Revision_C_Overleaf.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(U.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(U).as_posix())
(S/'verification.json').write_text(json.dumps({'build':str(B),'pages':pages,'pdf_pages':len(doc),'numeric_table_blocks_unchanged':True},indent=2))
report=C/'REPORT.md';s=report.read_text();s=s.replace('- **Chapter 4 tables:** only **two compact tables** remain in the main chapter: a four-row mean-error summary and five-row adopted-correction table. The four full fit tables move to Appendix C with their original labels and all numeric body entries unchanged. Captions and duplicate bounds detail are reduced in the main text.','- **Chapter 4 tables (restored at Ben\'s request):** separate complete pre-optimization, optimized pinned-flexor, and optimized pinned-extensor GoF tables are back in Chapter 4, with all values and scoring conditions preserved. The redundant mean-only summary is removed. The adopted corrections now have two rows (flexor/extensor) and four parameter columns (chi0–chi3). The earlier cross-configuration table remains in Appendix C. See review/table_restoration_20261023 for the diff, rendered pages, and build evidence.')
s=s.replace('Revision_C_final_20261023.pdf, 256 physical pages',f'Revision_C_tables_20261023.pdf, {len(doc)} physical pages')
s=s.replace('SADb occupies printed pages **110–115 (6 pages)**; SNS **116–127 (12 pages)**.','SADb and SNS retain their six-page and twelve-page scope; page numbers shift with the restored Chapter 4 tables.')
s=s.replace('The four moved fit-table bodies were compared numerically against the pre-sweep source: every entry matches. The mean summary is transcribed from the existing mean rows, not refitted.','The three restored GoF table blocks match their appendix versions exactly, including every numeric entry. The earlier cross-configuration table remains unchanged in Appendix C.')
report.write_text(s,encoding='utf-8')
readme=C/'README.md';s=readme.read_text().replace('Revision_C_final_20261023.pdf','Revision_C_tables_20261023.pdf');readme.write_text(s,encoding='utf-8')
print(json.dumps({'pdf':str(pdf),'pages':pages,'total_pages':len(doc)},indent=2))
