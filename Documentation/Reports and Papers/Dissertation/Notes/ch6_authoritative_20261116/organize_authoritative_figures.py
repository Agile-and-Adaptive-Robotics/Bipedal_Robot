"""Install reviewed vector figures and preserve superseded assets with hashes."""
from pathlib import Path
import hashlib
import json
import shutil

HERE = Path(__file__).resolve().parent
DISS = HERE.parent.parent
ROOT = DISS.parents[2]
FIGS = DISS/'Figures'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
old = {
    '60-sns':['biw2lvar_air_stage2_20261116','simulation_workflow_20261116',
        'sns_knee_reflex_20261116','sns_cup_reflex_20261116'],
    '95-AppendixD':['muscle_force_compare_current','knee_transient_response','animatlab_shared_rg_original_20261116',
        'w2l_perturbation_verified','syn6_pf_timing_selected_20261116',
        'biw2l_stimulus_neural_20261116','biw2l_stimulus_mechanics_20261116']}
ledger_path = HERE/'FIGURE_MOVES.json'
ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else {'moves':[], 'installed':[]}
for chapter, stems in old.items():
    parent = FIGS/chapter
    for stem in stems:
        for p in sorted(parent.glob(stem+'*')):
            if not p.is_file():
                continue
            dest = parent/'Deprecated/20261116_authoritative'/p.name
            p.resolve().relative_to(ROOT.resolve())
            dest.resolve().relative_to(ROOT.resolve())
            if dest.exists():
                raise FileExistsError(dest)
            dest.parent.mkdir(parents=True,exist_ok=True)
            h = sha(p)
            shutil.move(str(p),str(dest))
            assert sha(dest)==h
            ledger['moves'].append({'original_path':str(p),'deprecated_path':str(dest),'sha256':h})
            ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
        for dest in sorted((parent/'Deprecated/20261116_authoritative').glob(stem+'*')):
            if dest.is_file() and not any(x['deprecated_path']==str(dest) for x in ledger['moves']):
                ledger['moves'].append({'original_path':str(parent/dest.name),
                    'deprecated_path':str(dest),'sha256':sha(dest)})
                ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')

spec = [
    ('60-sns',HERE/'air_figure','biw2lvar_air_aligned_authoritative','biw2lvar_air_aligned_authoritative'),
    ('60-sns',HERE,'BUILD_workflow_20261116','simulation_workflow_authoritative'),
    ('60-sns',HERE/'CAD_sns_compact','sns_knee_reflex_authoritative','sns_knee_reflex_authoritative'),
    ('60-sns',HERE/'CAD_sns_compact','sns_cup_reflex_authoritative','sns_cup_reflex_authoritative'),
    ('95-AppendixD',HERE,'CAD_supported_traces_20261116','CAD_supported_traces'),
    ('95-AppendixD',HERE,'PF_timing_authoritative','PF_timing_authoritative'),
    ('95-AppendixD',HERE,'AIR_stimulus_neural','AIR_stimulus_neural'),
    ('95-AppendixD',HERE,'AIR_stimulus_mechanics','AIR_stimulus_mechanics')]
for chapter,srcdir,srcstem,dststem in spec:
    for ext in ['.pdf','.png','.svg','.alt.txt']:
        src=srcdir/(srcstem+ext)
        if not src.exists():
            continue
        dest=FIGS/chapter/(dststem+ext)
        if dest.exists() and sha(dest)!=sha(src):
            prior=next((x for x in ledger['installed'] if x['active_path']==str(dest)),None)
            if prior is None or sha(dest)!=prior['sha256']:
                raise FileExistsError('Refusing overwrite of independently edited active asset: '+str(dest))
            archived=FIGS/chapter/'Deprecated/20261116_authoritative/working_draft'/(dest.stem+'_'+sha(dest)[:12]+dest.suffix)
            archived.resolve().relative_to(ROOT.resolve())
            if archived.exists():
                raise FileExistsError(archived)
            archived.parent.mkdir(parents=True,exist_ok=True)
            oldhash=sha(dest)
            shutil.move(str(dest),str(archived))
            ledger['moves'].append({'original_path':str(dest),'deprecated_path':str(archived),
                'sha256':oldhash,'reason':'Session working draft compacted before final review'})
            ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
        shutil.copy2(src,dest)
        item={'source':str(src),'active_path':str(dest),'sha256':sha(dest)}
        ledger['installed']=[x for x in ledger['installed'] if x['active_path']!=str(dest)]+[item]
ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
print(json.dumps({'preserved':len(ledger['moves']),'installed':len(ledger['installed'])}))
