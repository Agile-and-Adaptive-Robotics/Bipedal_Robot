from pathlib import Path
import hashlib,json,re,difflib,importlib.util
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';P=D/'Notes/advisor_revision_20261022';C=P/'Revision_C'
spec=importlib.util.spec_from_file_location('words',P/'measure_words.py');w=importlib.util.module_from_spec(spec)
helper=(P/'measure_words.py').read_text(encoding='utf-8').replace('masked2[i] = " "','\n            if i < len(masked2): masked2[i] = " "')
exec(compile(helper,str(P/'measure_words.py'),'exec'),w.__dict__)
# The campaign helper strips escaped percent signs incorrectly; fix that for both sides.
w.strip_comments=lambda t:'\n'.join(re.split(r'(?<!\\)%',line,maxsplit=1)[0] for line in t.splitlines())
counts={};diffs=[];labelchanges={};numeric=[]
for p in sorted((C/'chapters').glob('*.tex')):
    before=C/'review/before'/p.name
    if not before.exists():
        proposal=P/'proposals'/p.name
        before=proposal if proposal.exists() else D/'ProofFinal/chapters'/p.name
    old=before.read_text(encoding='utf-8-sig') if before.exists() else ''
    new=p.read_text(encoding='utf-8-sig')
    if old==new:continue
    oldname=str(before.relative_to(R)) if before.exists() else '/dev/null'
    diffs.extend(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=oldname,tofile=str(p.relative_to(R))))
    counts[p.name]={'before':w.measure(str(before))[:3] if before.exists() else [0,0,0],'after':w.measure(str(p))[:3]}
    oldlabels=set(re.findall(r'\\label\{([^}]+)\}',old));newlabels=set(re.findall(r'\\label\{([^}]+)\}',new))
    labelchanges[p.name]={'removed_or_moved':sorted(oldlabels-newlabels),'added_or_moved':sorted(newlabels-oldlabels)}
    for line in difflib.unified_diff(old.splitlines(),new.splitlines(),n=0):
        if line.startswith('-') and not line.startswith('---') and re.search(r'\d',line):numeric.append({'file':p.name,'removed_or_replaced_line':line[1:]})
oldmain=(P/'baseline/main.tex').read_text(encoding='utf-8-sig');newmain=(C/'main.tex').read_text(encoding='utf-8-sig')
diffs.extend(difflib.unified_diff(oldmain.splitlines(True),newmain.splitlines(True),fromfile=str((P/'baseline/main.tex').relative_to(R)),tofile=str((C/'main.tex').relative_to(R))))
(C/'review/changes.diff').write_text(''.join(diffs),encoding='utf-8')
(C/'review/word_counts.json').write_text(json.dumps(counts,indent=2))
(C/'review/label_moves.json').write_text(json.dumps(labelchanges,indent=2))
(C/'review/numeric_change_review.json').write_text(json.dumps(numeric,indent=2))
alltext='\n'.join(p.read_text(encoding='utf-8-sig') for p in (C/'chapters').glob('*.tex'))
labels=re.findall(r'\\label\{([^}]+)\}',alltext);refs=re.findall(r'\\(?:eqref|ref|autoref)\{([^}]+)\}',alltext)
missing=sorted(set(refs)-set(labels));assert not missing,missing
original=json.loads((C/'review/baseline_hashes.json').read_text())
hashchanged=[name for name,digest in original.items() if hashlib.sha256((R/name).read_bytes()).hexdigest()!=digest]
assert not hashchanged,hashchanged
protected=json.loads((D/'Notes/figure_polish_20261022/protected_hashes.json').read_text())
for source,digest in protected.items():
    source=Path(source);assert hashlib.sha256(source.read_bytes()).hexdigest()==digest,source
    target=(C/'chapters'/source.name) if source.suffix=='.tex' else (C/'figs'/source.relative_to(D/'Figures'))
    if target.exists():assert hashlib.sha256(target.read_bytes()).hexdigest()==digest,target
assert (C/'chapters/20-muscles.tex').read_bytes()==(D/'ProofFinal/chapters/20-muscles.tex').read_bytes()
for p in (C/'figs/20-muscles').glob('*'):
    source=D/'Figures/20-muscles'/p.name
    if source.exists():assert p.read_bytes()==source.read_bytes(),p
assert (C/'thesis.bib').read_bytes()==(D/'ProofFinal/thesis.bib').read_bytes()
(C/'review/source_integrity.json').write_text(json.dumps({'source_files_changed':hashchanged,'missing_labels':missing,'protected_chapter_20_unchanged':True,'protected_figure_hashes_verified':len(protected),'copied_chapter_20_figures_unchanged':True,'bibliography_byte_identical':True},indent=2))
for name,vals in counts.items():print(name,vals)
print('Source integrity and label checks passed')
