"""Create a source-only review package and audit, independently of TeX startup."""
from pathlib import Path
import difflib
import hashlib
import json
import re
import shutil
import zipfile
from verify_and_package import ROOT, DISS, NOTES, SCRATCH, prepare

stage = SCRATCH/'review_source_delivery'
stage.mkdir(parents=True, exist_ok=True)
prepare(stage)
sources = list(stage.glob('*.tex'))+list((stage/'chapters').glob('*.tex'))
texts = {p: re.sub(r'(?<!\\)%[^\n]*','',p.read_text(encoding='utf-8',errors='replace')) for p in sources}
alltext = '\n'.join(texts.values())
labels = re.findall(r'\\label\{([^}]+)\}', alltext)
refs = re.findall(r'\\(?:ref|pageref|eqref|autoref|cref)\{([^}]+)\}', alltext)
unknown_refs = sorted(set(refs)-set(labels))
duplicates = sorted(k for k in set(labels) if labels.count(k)>1)
bibtext = '\n'.join(p.read_text(encoding='latin-1') for p in stage.glob('*.bib'))
bibkeys = set(re.findall(r'@\w+\s*\{\s*([^,\s]+)', bibtext))
cites = set(k.strip() for group in re.findall(r'\\cite\w*\*?(?:\[[^\]]*\]){0,2}\{([^}]+)\}', alltext) for k in group.split(','))
unknown_cites = sorted(cites-bibkeys)
active = set()
missing = []
for tex,text in texts.items():
    paths = re.findall(r'\{(figs/[^{}]+/)\}',text)
    for graphic in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^{}]+)\}',text):
        found = []
        for candidate in [stage/p/graphic for p in paths]+[stage/'figs'/graphic,stage/graphic]:
            if candidate.suffix:
                if candidate.exists(): found.append(candidate)
            else:
                found.extend(candidate.with_suffix(ext) for ext in ('.pdf','.png','.jpg','.jpeg','.eps') if candidate.with_suffix(ext).exists())
        if not found: missing.append([tex.name,graphic])
        for item in found:
            active.add(item)
            if item.suffix=='.eps':
                cache=item.with_name(item.stem+'-eps-converted-to.pdf')
                if cache.exists():active.add(cache)
audit={'scope':'full source static scan; not a compilation result',
       'unknown_references':unknown_refs,'duplicate_labels':duplicates,
       'unknown_citations':unknown_cites,'missing_graphics':missing,
       'active_graphic_files':len(active)}
if any([unknown_refs,duplicates,unknown_cites,missing]):
    print(json.dumps(audit,indent=2))
    raise RuntimeError('Static package gate failed')
output=ROOT/'output/packages'
output.mkdir(parents=True,exist_ok=True)
destination=output/'Bolen_Dissertation_Chapter6_Review_20261116.zip'
with zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(stage.rglob('*')):
        if not path.is_file(): continue
        rel=path.relative_to(stage)
        if rel.parts[0]=='figs':
            if path not in active: continue
        elif path.suffix not in ('.tex','.bib','.cls','.bst'):continue
        archive.write(path,rel.as_posix())
with zipfile.ZipFile(DISS/'Bolen_Dissertation.zip') as baseline:
    patches=[]
    changes=[]
    for member in baseline.namelist():
        if Path(member).suffix not in ('.tex','.bib','.cls','.bst'):continue
        before=baseline.read(member)
        after=(stage/member).read_bytes()
        if before.replace(b'\r\n',b'\n')==after.replace(b'\r\n',b'\n'):continue
        changes.append(member)
        patches.extend(difflib.unified_diff(before.decode('utf-8').splitlines(True),after.decode('utf-8').splitlines(True),fromfile='a/'+member,tofile='b/'+member))
    if set(changes)!={'chapters/60-sns.tex','chapters/95-AppendixD.tex','chapters/70-overall.tex'}:
        raise RuntimeError(f'Unexpected source scope: {changes}')
    diff=output/'Chapter6_AppendixD_20261116.diff'
    diff.write_text(''.join(patches),encoding='utf-8')
    preserved = {name:baseline.read(name)==(stage/name).read_bytes() for name in
                 ('chapters/20-muscles.tex','chapters/10-introduction.tex','chapters/30-torque.tex','thesis.bib')}
audit.update({'modified_source_members':changes,
              'source_zip':str(destination.relative_to(ROOT)),
              'source_zip_sha256':hashlib.sha256(destination.read_bytes()).hexdigest(),
              'diff':str(diff.relative_to(ROOT)),
              'protected_muscle_source_unchanged':preserved['chapters/20-muscles.tex'],
              'latest_introduction_unchanged':preserved['chapters/10-introduction.tex'],
              'torque_source_unchanged':preserved['chapters/30-torque.tex'],
              'bibliography_unchanged':preserved['thesis.bib']})
(NOTES/'STATIC_PACKAGE_RESULT.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
for name in ('EB475WS4_FOLLOWUP.txt','EASTEREGG2_FOLLOWUP.txt'):
    shutil.copy2(NOTES/name,output/name)
print(json.dumps(audit,indent=2))
