"""Ground diagnostics with explicit shared-static-law and coupling variants."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'shared_force_coupling_v1';out.mkdir()
base=json.loads((root/'swing_diagnostics_v1/swing_0p15.json').read_text());force=H.parent/'physics_fidelity/shared_force_law_functional_v1.json'
plans={}
for gain in (.6,1.2,2.4):plans['shared_c'+str(gain).replace('.','p')]={**base,'commissural_gain':gain,'force_law_parameters_file':str(force)}
plans['old_force_c2p4']={**base,'commissural_gain':2.4}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,plans.items()):
  rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}).get('accepted'),flush=True)
