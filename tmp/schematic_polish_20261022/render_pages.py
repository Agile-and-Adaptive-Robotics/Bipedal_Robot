from pathlib import Path
import re,json,subprocess,fitz,hashlib
H=Path(__file__).resolve().parent;B=Path((H/'build_path.txt').read_text());pdf=B/'main.pdf'
doc=fitz.open(pdf);lof=(B/'main.lof').read_text(errors='replace')
subjects={'roadmap':'Dissertation roadmap','actuator_technology_comparison':'Representative robotic actuator','pipeline':'Simulation-to-robot pipeline','sns_library_overview':'Simulink SNS library blocks','KneeReflex_circuit':'Reduced-order knee-reflex','sns_curriculum':'Five-stage training','design_workflow':'first-pass design workflow','bpa_construction':'Construction of the braided','sensorium_hierarchy':'Control stack'}
records=[];out=H/'pages';out.mkdir(exist_ok=True)
for name,key in subjects.items():
    entry=next((l for l in lof.splitlines() if key.lower() in l.lower()),None)
    if not entry:raise RuntimeError('Missing LOF '+name)
    match=re.search(r'numberline\s*\{([^}]+)\}',entry);num=match.group(1)
    candidates=[i for i,p in enumerate(doc) if re.search(r'Figure\s+'+re.escape(num)+r'\s*:',p.get_text())]
    if len(candidates)!=1:raise RuntimeError(str((name,num,candidates)))
    page=candidates[0]+1;dst=out/(name+'.png')
    subprocess.run([r'D:\MiKTeX\miktex\bin\x64\mgs.exe','-q','-dNOPAUSE','-dBATCH','-sDEVICE=png16m','-r150','-dTextAlphaBits=4','-dGraphicsAlphaBits=4',f'-dFirstPage={page}',f'-dLastPage={page}',f'-sOutputFile={dst}',str(pdf)],check=True)
    p=doc[page-1];spans=[span for block in p.get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']]
    bad=[s for s in spans if s['bbox'][3]>p.rect.height or s['bbox'][0]<0 or s['bbox'][2]>p.rect.width]
    records.append({'name':name,'figure':num,'pdf_page':page,'out_of_page_text':len(bad)})
log=(B/'main.log').read_text(errors='replace')
protected=json.loads((H/'protected_hashes.json').read_text());changed=[p for p,h in protected.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h]
summary={'build':str(B),'pages':len(doc),'figures':records,'undefined_refs':len(re.findall(r'undefined references|Reference .* undefined',log)),'undefined_citations':len(re.findall(r'Citation .* undefined',log)),'missing_files':len(re.findall(r'File .* not found',log)),'protected_changed':changed}
(H/'validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
