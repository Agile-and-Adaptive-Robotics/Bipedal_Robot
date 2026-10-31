"""Standing pose sensitivity, coupled-circuit ablation on ground, and support ladder."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'ground_final_checks_v1';out.mkdir()
base=json.loads((root/'height_tonic_v1/k5000_i12.json').read_text());base['balance']={**base['balance'],'bilateral_recruitment':True};cases={}
for angle in (0.,-5.):cases['standing_lumbar_'+str(angle)]={**base,'stage':3,'support_fraction':0.,'vertical_stiffness_scale':1.,'vertical_damping_multiplier':1.,'harness_lowering_m':0.,'start_lumbar_extension_deg':angle}
for support in (.4,.3):cases['support_'+str(support)]={**base,'reflex_scale':.5,'support_fraction':support}
for support in (.7,.5):cases['no_comm_support_'+str(support)]={**base,'commissural_gain':0.,'commissural_e_gain':0.,'timing_scale':.01,'support_fraction':support}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):
  rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}),flush=True)
