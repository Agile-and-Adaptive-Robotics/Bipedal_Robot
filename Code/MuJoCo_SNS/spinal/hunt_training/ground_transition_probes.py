"""Bounded stage-4 diagnostics for a left RG extension arrest; no pass assumed."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'ground_transition'
assert json.loads((root/'fresh_repeat/SUMMARY.json').read_text())['accepted']
base=json.loads((root/'ground_initial/force_support70.json').read_text())
if out.exists():raise SystemExit('Existing transition probes retained')
out.mkdir()
cases=(('onset_kick',dict(walk_kick=8.)),('reduced_timing',dict(timing_scale=.2,walk_kick=8.)),
       ('stronger_commissural',dict(commissural_gain=.6,walk_kick=8.)),('timing_off_diagnostic',dict(timing_afferents=False,walk_kick=8.)))
def run(item):
    name,extra=item;path=out/(name+'.json');path.write_text(json.dumps({**base,**extra},indent=2))
    with path.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT).returncode
    return name,json.loads(path.with_suffix('.result.json').read_text()) if code==0 else dict(accepted=False,exit=code)
with ThreadPoolExecutor(max_workers=2) as pool:
    rows={}
    for name,result in pool.map(run,cases):
        rows[name]=result;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(name,result.get('assessment',{}).get('accepted',False),flush=True)
shutil.copyfile(Path(__file__),out/Path(__file__).name)
