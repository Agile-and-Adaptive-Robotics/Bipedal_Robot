from pathlib import Path
import shutil,json,hashlib
R=Path(r'D:\Github\Bipedal_Robot');D=R/'Documentation/Reports and Papers/Dissertation';A=D/'Notes/advisor_revision_20261022';C=A/'Revision_C'
assert not C.exists(),'Revision_C already exists; do not overwrite'
(C/'chapters').mkdir(parents=True);(C/'review').mkdir()
for p in (D/'ProofFinal').iterdir():
    if p.is_file() and p.suffix.lower() in ['.tex','.bib','.bst','.cls','.sty']:
        shutil.copy2(p,C/p.name)
shutil.copy2(A/'baseline/main.tex',C/'main.tex')
hashes={}
for p in (D/'ProofFinal/chapters').glob('*.tex'):
    source=A/'proposals'/p.name
    if not source.exists():source=p
    shutil.copy2(source,C/'chapters'/p.name)
    hashes[str(source.relative_to(R))]=hashlib.sha256(source.read_bytes()).hexdigest()
for p in (D/'Figures').rglob('*'):
    if p.is_file() and not any(x.lower() in ['components','deprecated','__pycache__'] for x in p.relative_to(D/'Figures').parts[:-1]) and p.suffix.lower() in ['.png','.jpg','.jpeg','.pdf','.eps','.svg']:
        dst=C/'figs'/p.relative_to(D/'Figures');dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
(C/'review/baseline_hashes.json').write_text(json.dumps(hashes,indent=2))
for name in ['10-introduction','15-background','30-torque','50-database','60-sns']:
    target=C/'review/before'/f'{name}.tex';target.parent.mkdir(exist_ok=True);shutil.copy2(C/'chapters'/f'{name}.tex',target)
print('Created',C,'with',len(list((C/'figs').rglob('*.*'))),'figure files')
