"""Controlled timing/coupling ablation, standing sensory weighting, and support ladder."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'ground_regression_causes_v1';out.mkdir()
base=json.loads((root/'height_tonic_v1/k5000_i12.json').read_text());base['balance']={**base['balance'],'bilateral_recruitment':True};cases={}
for name,timing,comm in [('no_timing',0.,True),('low_timing',.01,True),('no_comm',.1,False),('low_timing_no_comm',.01,False)]:
 c={**base,'stage':2,'support_fraction':0.,'vertical_stiffness_scale':1.,'vertical_damping_multiplier':1.,'harness_lowering_m':0.,'timing_scale':timing}
 for key in ('restriction','height_offset','start_hip_deg','start_knee_deg','start_ankle_deg'):c.pop(key,None)
 if not comm:c.update(commissural_gain=0.,commissural_e_gain=0.)
 cases['air_'+name]=c
for name,kp,kd,weight in [('strong_somatic',4.,1.5,1.),('mixed',4.,1.5,.5),('orientation',4.,1.5,0.)]:
 c={**base,'stage':3,'support_fraction':0.,'vertical_stiffness_scale':1.,'vertical_damping_multiplier':1.,'harness_lowering_m':0.,'balance':{**base['balance'],'kp':kp,'kd':kd,'somatic_weight':weight}}
 cases['standing_'+name]=c
for support in (.6,.5):cases[f'support_{support}']={**base,'reflex_scale':.5,'support_fraction':support}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):
  rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}),flush=True)
