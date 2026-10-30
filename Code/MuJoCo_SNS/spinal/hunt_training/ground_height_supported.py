"""Hold pelvis height while assessing positive-Ib stance support and loading events."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'height_supported'
if out.exists():raise SystemExit('Existing height-supported probes retained')
out.mkdir();base=json.loads((root/'stance_support/ib_0p3.json').read_text());base.pop('support_fraction',None);base['force_controlled_support']=False
event=dict(cross_contact_gain=6.,cross_contact_event_amplitude=2.,cross_contact_event_decay_s=.15)
plans={'ib_spring':base,'ib_spring_event':{**base,**event},'ib_spring_lower':{**base,**event,'harness_vertical_offset_m':-.0002},
       'ib_lower_spring':{**base,**event,'weights_directory':str(root/'ground_contact_events/ib_lower_weights'),'harness_vertical_offset_m':-.0002}}
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
