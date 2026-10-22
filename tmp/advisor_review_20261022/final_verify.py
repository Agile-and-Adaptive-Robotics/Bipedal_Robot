from pathlib import Path
import json,re,hashlib,shutil,fitz
H=Path(__file__).resolve().parent;R=H.parent.parent;D=R/'Documentation/Reports and Papers/Dissertation';N=R/'Notes/advisor_review_20261022'
B=Path((H/'build_path.txt').read_text());doc=fitz.open(B/'main.pdf');log=(B/'main.log').read_text(errors='replace')
records=[]
for num in ['6.6','7.2']:
    hit=next(i for i,p in enumerate(doc) if 'Figure '+num+':' in p.get_text())
    page=doc[hit];page.get_pixmap(matrix=fitz.Matrix(150/72,150/72)).save(N/f'final_figure_{num}.png')
    spans=[s for b in page.get_text('dict')['blocks'] if 'lines' in b for l in b['lines'] for s in l['spans']]
    records.append({'figure':num,'pdf_page':hit+1,'minimum_figure_text_pt':min(s['size'] for s in spans if 'Arial' in s['font']),'lowest_text_y':max(s['bbox'][3] for s in spans),'page_height':page.rect.height})
protected=json.loads((D/'Notes/figure_polish_20261022/protected_hashes.json').read_text())
changed=[p for p,h in protected.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h]
live=[]
for n,file in [('60','60-sns.tex'),('70','70-overall.tex')]:
    expected=(H/'overleaf_expected'/file).read_text(encoding='utf-8-sig');actual=(H/f'final{n}.tex').read_text(encoding='utf-8-sig')
    live.append({'file':file,'matches_expected_targeted_edits':expected==actual})
    if expected!=actual:
        import difflib
        (H/f'unexpected_{n}.diff').write_text(''.join(difflib.unified_diff(expected.splitlines(True),actual.splitlines(True))),encoding='utf-8')
summary={'pages':len(doc),'figures':records,'undefined_refs':len(re.findall(r'undefined references|Reference .* undefined',log)),'undefined_citations':len(re.findall(r'Citation .* undefined',log)),'missing_files':len(re.findall(r'File .* not found',log)),'protected_changed':changed,'live_source_checks':live}
(N/'final_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
assert not changed
OUT=D/'ProofFinal/build_advisor_figures_20261022_final';OUT.mkdir(exist_ok=True)
for p in B.iterdir():
    if p.is_file():shutil.copy2(p,OUT/p.name)
if (B/'chapters').exists():shutil.copytree(B/'chapters',OUT/'chapters',dirs_exist_ok=True)
for n in ['review_final.json','highlights_final_additions.json','final_caption_patches.json','background_final.tex']:
    shutil.copy2(H/n,N/n)
print('SAVED',OUT)
