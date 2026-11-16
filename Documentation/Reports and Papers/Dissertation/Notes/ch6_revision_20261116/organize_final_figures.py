"""Keep Final's referenced chapter figures active; preserve old assets in Deprecated."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil

ROOT=Path(__file__).resolve().parents[5]
DISS=ROOT/'Documentation/Reports and Papers/Dissertation'
NOTES=DISS/'Notes/ch6_revision_20261116'
ACTIVE={
 '60-sns': {'OurWalker','LI_DENG_FALL_FRAME','simulation_workflow_20261116',
            'sns_library_20261116','sns_knee_reflex_20261116','sns_cup_reflex_20261116',
            'biw2lvar_air_stage2_20261116','w2l_feedback_short_20261031','cad_supported_20261116'},
 '95-AppendixD': {'muscle_force_compare_current','w2l_crossengine_20261116',
                 'knee_transient_response','animatlab_shared_rg_original_20261116',
                 'animatlab_short_response_20261116','biw2l_stimulus_neural_20261116',
                 'biw2l_stimulus_mechanics_20261116','syn6_pf_timing_selected_20261116',
                 'gait2392_ik_overlay','w2l_perturbation_verified',
                 'spiking_ground_neuromechanics_short_20261031'}
}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def key(p):
    name=p.name
    for suffix in ('.alt.txt','_alt.txt'):
        if name.endswith(suffix):return name[:-len(suffix)]
    return p.stem

def confined(path):
    resolved=path.resolve()
    if not resolved.is_relative_to(ROOT.resolve()):raise RuntimeError(f'Outside repository: {resolved}')
    return resolved

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--apply',action='store_true')
    args=parser.parse_args()
    moves=[]
    for folder,stems in ACTIVE.items():
        directory=DISS/'Figures'/folder
        for source in sorted(directory.iterdir()):
            if source.is_file() and key(source) not in stems:
                destination=directory/'Deprecated/20261116'/source.name
                confined(source);confined(destination)
                if destination.exists():raise RuntimeError(f'Deprecation destination already exists: {destination}')
                moves.append({'source':str(source.relative_to(ROOT)),
                              'destination':str(destination.relative_to(ROOT)),
                              'sha256':digest(source)})
    recovery=DISS/'Figures/60-sns/LI_DENG_FALL_FRAME.png'
    original=ROOT/'tmp/ch6_20261116/baseline/figs/60-sns/LI_DENG_FALL_FRAME.png'
    if not original.is_file():raise RuntimeError('Missing recovery baseline')
    recovery_record={'source':str(original.relative_to(ROOT)),
                     'destination':str(recovery.relative_to(ROOT)),
                     'sha256':digest(original),'needed':not recovery.exists()}
    if recovery.exists() and digest(recovery)!=digest(original):
        raise RuntimeError('Recovery image already exists with different contents')
    if args.apply:
        if not recovery.exists():shutil.copy2(original,confined(recovery))
        for item in moves:
            source=confined(ROOT/item['source'])
            destination=confined(ROOT/item['destination'])
            destination.parent.mkdir(parents=True,exist_ok=True)
            shutil.move(str(source),str(destination))
            if digest(destination)!=item['sha256']:raise RuntimeError(f'Moved content changed: {destination}')
    # Resolve active includes in every Final chapter, including inherited search paths.
    missing=[]
    for tex in (DISS/'Final/chapters').glob('*.tex'):
        text=re.sub(r'(?<!\\)%[^\n]*','',tex.read_text(encoding='utf-8',errors='replace'))
        paths=re.findall(r'\{([^{}]+/)\}',re.search(r'\\graphicspath\{(.*)\}',text).group(1)) if '\\graphicspath' in text else []
        for figure in re.findall(r'\\includegraphics(?:\[[^\]]*\])?\{([^{}]+)\}',text):
            candidates=[DISS/'Final'/p/figure for p in paths]+[DISS/'Final'/figure]
            found=any(p.is_file() or (not p.suffix and any(p.with_suffix(e).is_file() for e in ('.pdf','.png','.jpg','.jpeg','.eps'))) for p in candidates)
            if not found:missing.append({'chapter':tex.name,'figure':figure})
    manifest_path=NOTES/'FIGURE_ORGANIZATION.json'
    if args.apply and manifest_path.exists() and not moves:
        previous=json.loads(manifest_path.read_text(encoding='utf-8'))
        for item in previous['moves']:
            destination=confined(ROOT/item['destination'])
            if not destination.is_file() or digest(destination)!=item['sha256']:
                raise RuntimeError(f'Prior deprecated asset changed: {destination}')
        previous['missing_final_graphics']=missing
        previous['repeat_verification']='Prior move ledger and hashes preserved; no additional moves required.'
        manifest_path.write_text(json.dumps(previous,indent=2),encoding='utf-8')
        print(json.dumps({'mode':'verified existing organization','moved_count':previous['moved_count'],
                          'missing_final_graphics':missing},indent=2))
        if missing:raise RuntimeError('Missing Final figures remain')
        return
    report={'mode':'applied' if args.apply else 'plan',
            'canonical_sources':['Final/chapters/60-sns.tex','Final/chapters/95-AppendixD.tex'],
            'moved_count':len(moves),'moves':moves,'recovery':recovery_record,
            'missing_final_graphics':missing,
            'historical_sources':'ProofFinal source files are preserved; their old figure paths require the Deprecated copies if reproducing historical builds.'}
    (NOTES/('FIGURE_ORGANIZATION.json' if args.apply else 'FIGURE_ORGANIZATION_PLAN.json')).write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    if args.apply and missing:raise RuntimeError('Missing Final figures remain')

if __name__=='__main__':main()
