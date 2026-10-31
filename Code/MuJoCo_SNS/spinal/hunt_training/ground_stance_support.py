"""Separate inadequate stance recruitment from the crossed swing trigger.

Diagnostic changes use existing PF/IBEXC neurons and isolated weight copies.
Any ground finalist must return through both air-stage checks before promotion.
"""
import os,sys,json,shutil,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/tonic_abduction_20261030_s1_v1/0p02';out=root/'stance_support'
if out.exists():raise SystemExit('Existing stance-support probes retained')
out.mkdir();base=json.loads((root/'ground_initial/force_support70.json').read_text());plans={}
for family,increment in (('ib',.1),('ib',.3),('knee_pf',.1),('knee_pf',.3)):
    name=family+'_'+str(increment).replace('.','p');folder=out/name;weights=folder/'weights';weights.mkdir(parents=True)
    for filename in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json','circuit_options.json'):shutil.copyfile(root/'weights'/filename,weights/filename)
    filename='reflex_weights.json' if family=='ib' else 'mn_weights.json';path=weights/filename;z=json.loads(path.read_text());changes=[]
    for key,gain in z['edges'].items():
        pre,post=key.split(' -> ')
        match=pre.startswith(('IBEXC_hip_ext_','IBEXC_knee_ext_','IBEXC_ankle_pf_')) if family=='ib' else pre.startswith('PF_KNEE-E_') and post[3:].startswith(('vas_med_','vas_int_','vas_lat_'))
        if match:z['edges'][key]=gain+increment;changes.append(dict(edge=key,old=gain,new=gain+increment))
    assert changes,(family,increment)
    z['experimental_stance_support']=dict(kind=family,increment=increment,changes=changes);path.write_text(json.dumps(z,indent=2))
    cfg={**base,'weights_directory':str(weights)};plans[name]=cfg
for gain in (40.,100.):plans[f'contra_heel_{gain:g}']={**base,'cross_contact_gain':gain,'export_graph':True}
(out/'PLAN.json').write_text(json.dumps(plans,indent=2))
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
