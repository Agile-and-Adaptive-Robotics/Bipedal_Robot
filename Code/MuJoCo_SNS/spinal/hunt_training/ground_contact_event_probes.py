"""Measured crossed foot-loading events with the existing load-sharing candidate."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'ground_contact_events'
if out.exists():raise SystemExit('Existing contact-event probes retained')
out.mkdir();base=json.loads((root/'stance_support/ib_0p3.json').read_text())
plans={f'cross_{g:g}':{**base,'cross_contact_gain':g,'cross_contact_event_amplitude':2.,'cross_contact_event_decay_s':.15,'export_graph':True} for g in (2.,6.,12.)}
plans['ib_lower']={**plans['cross_6']}
weights=out/'ib_lower_weights';weights.mkdir()
for filename in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json','circuit_options.json'):shutil.copyfile(root/'weights'/filename,weights/filename)
z=json.loads((weights/'reflex_weights.json').read_text())
for key in z['edges']:
    if key.startswith(('IBEXC_hip_ext_','IBEXC_knee_ext_','IBEXC_ankle_pf_')):z['edges'][key]+=.25
(weights/'reflex_weights.json').write_text(json.dumps(z,indent=2));plans['ib_lower']['weights_directory']=str(weights)
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
