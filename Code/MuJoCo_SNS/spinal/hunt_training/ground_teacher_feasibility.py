"""Ground feasibility oracle only; clamped teacher results are never neural passes."""
import sys,json,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'ground_teacher_v1';out.mkdir();base=json.loads((root/'corrected_initializer_v1/same_controller.json').read_text())
cases={'slow_pd':dict(period_s=2.466,kp=180.,kd=12.),'slow_stiff_pd':dict(period_s=2.466,kp=400.,kd=20.),'normal_pd':dict(period_s=1.233,kp=400.,kd=20.)}
def run(item):
 name,teacher=item;p=out/(name+'.json');cfg={**base,'motor_teacher':dict(mode='feedback',stride=5,**teacher)};p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment',{}).get('teacher_motion_gate'),flush=True)
