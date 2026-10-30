"""Corrected stage 4 with the same passing controller weights and real ground."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'ground_initial'
assert json.loads((root/'fresh_repeat/SUMMARY.json').read_text())['accepted']
standing=json.loads((root/'standing/task_low.result.json').read_text());assert standing['assessment']['accepted']
if out.exists():raise SystemExit('Existing stage-4 probes retained')
out.mkdir();base={**standing['configuration'],'stage':4}
cases=(('spring_support',{}),('force_support70',dict(force_controlled_support=True,support_fraction=.7)))
def run(item):
    name,extra=item;path=out/(name+'.json');path.write_text(json.dumps({**base,**extra},indent=2))
    with path.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT).returncode
    return name,json.loads(path.with_suffix('.result.json').read_text()) if code==0 else dict(accepted=False,exit=code)
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,cases))
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));shutil.copyfile(Path(__file__),out/Path(__file__).name)
print('GROUND_COMPLETE',json.dumps({name:r.get('assessment',{}).get('accepted',False) for name,r in rows.items()}),flush=True)
