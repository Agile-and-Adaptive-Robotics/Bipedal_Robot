"""Stage-3 standing only after the same candidate passes both air stages twice."""
import os,sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'standing'
for folder in ('air_pair','fresh_repeat'):
    assert json.loads((root/folder/'SUMMARY.json').read_text())['accepted'],folder
if out.exists():raise SystemExit('Existing standing tests retained')
out.mkdir();base=json.loads((root/'fresh_repeat/shared_configuration.json').read_text())
base.update(stage=3,restriction=1.,height_offset=-.02,start_hip_deg=-7.,start_knee_deg=-10.,start_ankle_deg=0.)
def run(item):
    name,gains=item;path=out/(name+'.json')
    cfg={**base,'balance':dict(kp=gains[0],kd=gains[1],ka=0.,delay_s=.1,replace_inherited=False)}
    path.write_text(json.dumps(cfg,indent=2))
    with path.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT).returncode
    return name,json.loads(path.with_suffix('.result.json').read_text()) if code==0 else dict(accepted=False,exit=code)
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,(('task_off',(0.,0.)),('task_low',(.5,.2)))))
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));shutil.copyfile(Path(__file__),out/Path(__file__).name)
print('STANDING_COMPLETE',json.dumps({name:r.get('assessment',{}).get('accepted',False) for name,r in rows.items()}),flush=True)
