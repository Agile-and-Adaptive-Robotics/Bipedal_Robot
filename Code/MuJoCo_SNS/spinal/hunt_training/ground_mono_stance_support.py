"""Validate selective monoarticular stance recruitment before ground probes.

Retain the passing early-stage candidate. Add conductance only on existing
IBEXC-group outputs to gluteus maximus, vasti, and soleus, avoiding extra
rectus-femoris hip flexion and gastrocnemius knee flexion. These group cells
also receive RG-E drive; this is not a pure Ib-force intervention.
"""
import sys,json,shutil,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent
root=H/'runs/tonic_abduction_20261030_s1_v1/0p02'
out=root/'mono_stance_support'
if out.exists():raise SystemExit('Existing mono-stance probes retained')
out.mkdir();shutil.copyfile(Path(__file__),out/Path(__file__).name)
base=json.loads((root/'ground_initial/force_support70.json').read_text())
standing=json.loads((root/'standing/task_low.json').read_text())
variants={};rows={}
for increment in (.3,.6):
    name='mono_ib_'+str(increment).replace('.','p');folder=out/name
    weights=folder/'weights';weights.mkdir(parents=True)
    for filename in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json','circuit_options.json'):
        shutil.copyfile(root/'weights'/filename,weights/filename)
    path=weights/'reflex_weights.json';z=json.loads(path.read_text());changes=[]
    for edge,value in z['edges'].items():
        pre,post=edge.split(' -> ')
        match=(pre.startswith('IBEXC_hip_ext_') and post.startswith('MN_glut_max') or
               pre.startswith('IBEXC_knee_ext_') and post.startswith(('MN_vas_med_','MN_vas_int_','MN_vas_lat_')) or
               pre.startswith('IBEXC_ankle_pf_') and post.startswith('MN_soleus_'))
        if match:z['edges'][edge]=value+increment;changes.append(dict(edge=edge,old=value,new=value+increment))
    assert len(changes)==14,len(changes)
    z['selective_stance_support']=dict(increment=increment,changes=changes,
        scope='Existing phase-and-load group cells; not pure force-only feedback.')
    path.write_text(json.dumps(z,indent=2));variants[name]=weights
    with (folder/'air_pair.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(weights),'--output-directory',str(folder/'air_pair')],stdout=log,stderr=subprocess.STDOUT).returncode
    pair=json.loads((folder/'air_pair/SUMMARY.json').read_text()) if code==0 else dict(accepted=False,exit_code=code)
    rows[name]=dict(air_pair=pair);print(name,'air',pair['accepted'],flush=True)
    (out/'SUMMARY.json').write_text(json.dumps(rows,indent=2))
def physical(item):
    name,cfg=item;path=out/(name+'.json');path.write_text(json.dumps(cfg,indent=2))
    with path.with_suffix('.log').open('w') as log:
        code=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(path)],stdout=log,stderr=subprocess.STDOUT).returncode
    result=json.loads(path.with_suffix('.result.json').read_text()) if code==0 else dict(assessment=dict(accepted=False),exit_code=code)
    return name,result
plans={name:{**standing,'weights_directory':str(weights)} for name,weights in variants.items() if rows[name]['air_pair']['accepted']}
with ThreadPoolExecutor(max_workers=2) as pool:
    for name,result in pool.map(physical,plans.items()):
        rows[name]['standing']=result;print(name,'standing',result['assessment']['accepted'],flush=True)
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2))
plans={}
for name,weights in variants.items():
    if not rows[name].get('standing',{}).get('assessment',{}).get('accepted',False):continue
    plans[name+'_ground']={**base,'weights_directory':str(weights)}
    plans[name+'_events']={**base,'weights_directory':str(weights),'cross_contact_gain':6.,'cross_contact_event_amplitude':2.,'cross_contact_event_decay_s':.15}
with ThreadPoolExecutor(max_workers=2) as pool:
    for name,result in pool.map(physical,plans.items()):
        rows[name]=result;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2))
        print(name,'ground',result['assessment']['accepted'],flush=True)
