"""Test whether unilateral lumbar phase drive produces unequal swing clearance."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'bilateral_midline_drive_v1';out.mkdir()
base=json.loads((root/'ground_final_checks_v1/no_comm_support_0.5.json').read_text());cases={}
for support in (.7,.5,.4):cases['support_'+str(support)]={**base,'bilateral_midline_drive':True,'support_fraction':support}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment'),r.get('foot_clearance'),flush=True)
