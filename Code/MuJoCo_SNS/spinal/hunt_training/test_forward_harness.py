"""Explicit forward harness guidance; results are assisted, not autonomous gait."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'forward_harness_v1';out.mkdir();base=json.loads((root/'shared_force_v3_ground/corrected_init.json').read_text())
cases={'v008':dict(harness_forward_speed=.08),'v016':dict(harness_forward_speed=.16),'v016_tonic08':dict(harness_forward_speed=.16,rg_tonic_current=.8)}
def run(item):
 name,extra=item;p=out/(name+'.json');p.write_text(json.dumps({**base,**extra},indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}).get('accepted'),flush=True)
