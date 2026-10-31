"""Lower weight support after fixing initial contacts and knee extension stops."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'aligned_support_ladder_v1';out.mkdir()
base=json.loads((root/'contact_feasible_start_v1/ground_guide.json').read_text());cases={'repeat_0p5':base}
for support in (.3,.1):cases['support_'+str(support).replace('.','p')]={**base,'support_fraction':support}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment'),r.get('foot_clearance'),flush=True)
