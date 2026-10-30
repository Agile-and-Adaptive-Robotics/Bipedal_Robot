"""Determine whether feasible, low-gain task torque can produce the requested air motion."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
H=Path(__file__).resolve().parent;B.O=H/'runs/swing_release_20261030_v1/sagittal_teachers';B.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/fourcell_dynamic_20261029_v1/curriculum/stage2.json').read_text());base['restriction']=.05;base['diagnose_swing']=True
cases=[('slow_low',2.466,40.,8.),('fast_low',1.233,40.,8.),('fast_medium',1.233,80.,12.)]
def run(item):
    name,T,kp,kd=item;r=B.run(name,{**base,'motor_teacher':dict(mode='feedback',period_s=T,kp=kp,kd=kd,sagittal_only=True,stride=5)})
    print(name,json.dumps(r['assessment']) if r else None,flush=True)
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(run,cases))
