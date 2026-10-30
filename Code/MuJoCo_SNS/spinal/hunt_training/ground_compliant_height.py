"""Combine prescribed weight support with a damped compliant height tether."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'compliant_height'
if out.exists():raise SystemExit('Existing compliant-height probes retained')
assert json.loads((root/'stance_support/ib_0p3/air_pair/SUMMARY.json').read_text())['accepted']
out.mkdir();base=json.loads((root/'stance_support/ib_0p3.json').read_text());base['force_controlled_support']=False
plans={'k04_s70':{**base,'vertical_stiffness_scale':.04},'k04_s65':{**base,'vertical_stiffness_scale':.04,'support_fraction':.65},
       'k02_s70':{**base,'vertical_stiffness_scale':.02},'k04_s70_event':{**base,'vertical_stiffness_scale':.04,'cross_contact_gain':6.,'cross_contact_event_amplitude':2.,'cross_contact_event_decay_s':.15}}
def run(item):
    name,cfg=item;path=out/(name+'.json');path.write_text(json.dumps(cfg,indent=2))
    with path.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT).returncode
    return name,json.loads(path.with_suffix('.result.json').read_text()) if code==0 else dict(accepted=False,exit=code)
with ThreadPoolExecutor(max_workers=2) as pool:
    rows={}
    for name,result in pool.map(run,plans.items()):
        rows[name]=result;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(name,result.get('assessment',{}).get('accepted',False),flush=True)
shutil.copyfile(Path(__file__),out/Path(__file__).name)
