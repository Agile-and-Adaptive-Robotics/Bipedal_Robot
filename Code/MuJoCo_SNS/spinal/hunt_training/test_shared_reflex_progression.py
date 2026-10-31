"""Check one common reflex attenuation with mirrored balance recruitment."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'shared_reflex_progression_v1';out.mkdir()
base=json.loads((root/'height_tonic_v1/k5000_i12.json').read_text());cases={}
for gain in (.5,.25):
 for stage in (2,3,4):
  c={**base,'stage':stage,'reflex_scale':gain,'balance':{**base['balance'],'bilateral_recruitment':True}}
  if stage<4:c.update(support_fraction=0.,force_controlled_support=False,vertical_stiffness_scale=1.,vertical_damping_multiplier=1.,harness_lowering_m=0.)
  if stage==2:
   for key in ('restriction','height_offset','start_hip_deg','start_knee_deg','start_ankle_deg'):c.pop(key,None)
  cases[f'reflex_{gain}_s{stage}']=c
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):
  rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}),flush=True)
