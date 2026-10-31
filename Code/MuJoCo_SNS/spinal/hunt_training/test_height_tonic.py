"""Combine regular lowered-height stepping with shared neural cadence drive."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'height_tonic_v1';out.mkdir()
cases={'k1000_i08':('neural_k1000',.8),'k1000_i12':('neural_k1000',1.2),'k5000_i08':('neural_k5000',.8),'k5000_i12':('neural_k5000',1.2)}
def run(item):
 name,(source,current)=item;cfg=json.loads((root/'lowered_height_guide_v1'/(source+'.json')).read_text());cfg['rg_tonic_current']=current;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}).get('accepted'),flush=True)
