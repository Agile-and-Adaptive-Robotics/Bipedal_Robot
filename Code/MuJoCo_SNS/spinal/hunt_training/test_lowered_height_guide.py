"""Declared, damped harness-height guidance with measured net assistance."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'lowered_height_guide_v1';out.mkdir();base=json.loads((root/'shared_force_v3_ground/corrected_init.json').read_text())
common=dict(force_controlled_support=False,vertical_damping_multiplier=4.2,harness_lowering_m=-.04,harness_lowering_start_s=2.,harness_lowering_duration_s=2.)
cases={'neural_k1000':dict(vertical_stiffness_scale=.002),'neural_k5000':dict(vertical_stiffness_scale=.01),'teacher_k5000':dict(vertical_stiffness_scale=.01,motor_teacher=dict(mode='feedback',period_s=2.466,kp=400.,kd=20.,stride=5,sagittal_only=True))}
def run(item):
 name,extra=item;p=out/(name+'.json');p.write_text(json.dumps({**base,**common,**extra},indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}).get('accepted'),flush=True)
