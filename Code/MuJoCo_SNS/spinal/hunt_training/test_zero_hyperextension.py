"""Revalidate knee stops at every saved physics step, with no pose clipping."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'zero_hyperextension_v1';out.mkdir()
cases={}
for name,source in [('ground_original','height_tonic_v1/k5000_i12.json'),('ground_no_comm','ground_final_checks_v1/no_comm_support_0.7.json'),('ground_half_support','ground_final_checks_v1/no_comm_support_0.5.json'),('standing','ground_final_checks_v1/standing_lumbar_0.0.json')]:cases[name]=json.loads((root/source).read_text())
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment'),r.get('foot_clearance'),flush=True)
