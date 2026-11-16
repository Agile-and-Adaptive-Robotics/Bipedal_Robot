"""Static source gates; explicitly not a successful compilation."""
from pathlib import Path
from collections import Counter
import hashlib, json, re

ROOT=Path(__file__).resolve().parents[5]
DISS=ROOT/'Documentation/Reports and Papers/Dissertation'
FINAL=DISS/'Final'
NOTES=Path(__file__).resolve().parent

def uncomment(text):
    return '\n'.join(re.split(r'(?<!\\)%',line,maxsplit=1)[0] for line in text.splitlines())

sources={p:uncomment(p.read_text(encoding='utf-8',errors='replace')) for p in FINAL.rglob('*.tex')}
alltext='\n'.join(sources.values())
labels=Counter(re.findall(r'\\label\{([^}]+)\}',alltext))
refs=[]
for match in re.finditer(r'\\(?:ref|eqref|pageref|autoref|cref|Cref)\*?(?:\[[^]]*\])?\{([^}]+)\}',alltext):
    refs.extend(x.strip() for x in match[1].split(','))
bibfiles=list(FINAL.glob('*.bib'))
bibkeys=set()
for p in bibfiles:
    bibkeys.update(re.findall(r'@\w+\s*\{\s*([^,\s]+)\s*,',p.read_bytes().decode('latin-1')))
cites=[]
for match in re.finditer(r'\\cite\w*\*?(?:\[[^]]*\]){0,2}\{([^}]+)\}',alltext):
    cites.extend(x.strip() for x in match[1].split(','))
graphics=[]
for path,text in sources.items():
    dirs=[FINAL]
    paths=re.search(r'\\graphicspath\{((?:\{[^{}]*\})+)\}',text)
    if paths:
        dirs.extend(FINAL/p for p in re.findall(r'\{([^{}]*)\}',paths[1]))
    dirs.extend([DISS/'Figures'])
    for graphic in re.findall(r'\\includegraphics\*?(?:\[[^]]*\])?\{([^{}]+)\}',text):
        trials=[]
        for directory in dirs:
            item=directory/graphic
            trials.extend([item] if item.suffix else [item.with_suffix(x) for x in ('.pdf','.png','.jpg','.jpeg','.eps')])
        found=[str(p.resolve()) for p in trials if p.is_file()]
        graphics.append({'source':str(path),'graphic':graphic,'found':list(dict.fromkeys(found))})
missing_inputs=[]
for p,text in sources.items():
    for item in re.findall(r'\\(?:input|include)\{([^}]+)\}',text):
        trial=FINAL/item
        if not trial.suffix: trial=trial.with_suffix('.tex')
        if not trial.exists(): missing_inputs.append({'source':str(p),'input':item})
report={'status':'static_checks_only','canonical_entry':str(FINAL/'main.tex'),
        'source_count':len(sources), 'label_count':len(labels),'citation_key_count':len(set(cites)),
        'undefined_source_references':sorted(set(refs)-set(labels)-{'LastPage'}),
        'missing_citations':sorted(set(cites)-bibkeys),
        'duplicate_labels':{k:v for k,v in labels.items() if v>1},
        'missing_graphics':[x for x in graphics if not x['found']],
        'missing_input_files':missing_inputs,
        'active_graphics':graphics,
        'final_directory_files':[str(p.relative_to(FINAL)) for p in sorted(FINAL.rglob('*')) if p.is_file()],
        'source_hashes':{str(p.relative_to(FINAL)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}}
report['all_static_gates_pass']=not any(report[key] for key in ('undefined_source_references','missing_citations','duplicate_labels','missing_graphics','missing_input_files'))
(NOTES/'BUILD_STATIC_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ('active_graphics','source_hashes','final_directory_files')},indent=2))
