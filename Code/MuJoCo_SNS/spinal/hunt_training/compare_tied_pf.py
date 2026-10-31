"""Test symmetric PF inhibitory conductances on isolated copies of the same fit."""
import os,sys,json,hashlib,shutil,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np

H=Path(__file__).resolve().parent;source=H/'runs/native_joint_pf_mono_20261030_v1/weights'
out=H/'runs/tied_pf_20261030_s1_v1'
if out.exists():raise SystemExit('Existing experiment retained; choose a new name in a reviewed copy')
out.mkdir();plans={}
for method in ('arithmetic','geometric'):
    root=out/method;weights=root/'weights';weights.mkdir(parents=True)
    for name in ('pf_weights.json','bi_weights.json','mn_weights.json','reflex_weights.json','timing_weights.json','circuit_options.json'):
        shutil.copyfile(source/name,weights/name)
    value=json.loads((weights/'pf_weights.json').read_text());changes=[]
    for side in ('r','l'):
        for joint in ('HIP','KNEE'):
            keys=[f'PF_IN_{joint}-E_{side} -> PF_{joint}-F_{side}',f'PF_IN_{joint}-F_{side} -> PF_{joint}-E_{side}']
            old=np.array([value['edges'][k] for k in keys]);shared=float(old.mean() if method=='arithmetic' else np.sqrt(old.prod()))
            for key in keys:value['edges'][key]=shared
            changes.append(dict(edges=keys,original=old.tolist(),shared=shared))
    value['experimental_constraint']=dict(type='Equal opposite inhibitory output conductances within each joint PF pair',method=method,changes=changes)
    (weights/'pf_weights.json').write_text(json.dumps(value,indent=2));plans[method]=dict(source=str(source),changes=changes)
(out/'PLAN.json').write_text(json.dumps(plans,indent=2));shutil.copyfile(Path(__file__),out/Path(__file__).name)
def run(method):
    root=out/method
    with (root/'pair.log').open('w') as log:
        result=subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(root/'weights'),'--output-directory',str(root/'air_pair')],stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:return dict(accepted=False,exit=result.returncode)
    return json.loads((root/'air_pair/SUMMARY.json').read_text())
# Each pair launcher runs two fresh processes; run methods sequentially to cap load.
rows={method:run(method) for method in plans}
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(json.dumps({k:v['accepted'] for k,v in rows.items()}),flush=True)
