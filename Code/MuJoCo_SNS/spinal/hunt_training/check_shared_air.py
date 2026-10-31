"""Fresh-process stages 1 and 2 with one saved configuration and weight set."""
import os,sys,json,hashlib,argparse,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

H=Path(__file__).resolve().parent
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--weights',type=Path,required=True);p.add_argument('--output-directory',type=Path,required=True)
p.add_argument('--config',type=Path);p.add_argument('--backend',choices=('sparse','native'),default='sparse')
a=p.parse_args();out=a.output_directory.resolve();weights=a.weights.resolve()
if out.exists():p.error('Choose a new output directory')
if not all((weights/name).is_file() for name in ('pf_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json')):p.error('Missing fitted weights')
out.mkdir(parents=True)
base=json.loads(a.config.read_text()) if a.config else dict(biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
base.pop('stage',None);base.update(weights_directory=str(weights),backend=a.backend,diagnose_swing=True)
if base.get('motor_teacher'):p.error('A motor teacher cannot be a neural curriculum candidate')
def hashes():return {path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in weights.glob('*.json')}
before=hashes();(out/'shared_configuration.json').write_text(json.dumps(base,indent=2))
def run(stage):
    path=out/f'stage{stage}.json';path.write_text(json.dumps({**base,'stage':stage},indent=2))
    with path.with_suffix('.log').open('w') as log:
        completed=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT)
    if completed.returncode:return dict(accepted=False,exit_code=completed.returncode,log=str(path.with_suffix('.log')))
    result=json.loads(path.with_suffix('.result.json').read_text());return dict(assessment=result['assessment'],measured=result['measured'],configuration=result['configuration'])
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(zip(('1','2'),pool.map(run,(1,2))))
unchanged=hashes()==before
report=dict(shared_configuration=base,weight_sha256=before,weights_unchanged=unchanged,stages=rows,
    accepted=bool(unchanged and all(r.get('assessment',{}).get('accepted',False) for r in rows.values())),
    scope='Same saved parameters; evaluator applies the required stage-specific sensory/interleg ablations. No later-stage promotion.')
(out/'SUMMARY.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
