"""Capture successful air teachers as targets; do not qualify them as neural motion."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
H=Path(__file__).resolve().parent;B.O=H/'runs/swing_release_20261030_v1/feasible_targets';B.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/fourcell_dynamic_20261029_v1/curriculum/stage2.json').read_text());base['restriction']=.05
def run(item):
    name,T=item;r=B.run(name,{**base,'motor_teacher':dict(mode='feedback',period_s=T,kp=40.,kd=8.,sagittal_only=True,stride=5)})
    assert r and r['assessment'].get('teacher_motion_gate'),(name,r)
    print(name,'TEACHER_TARGET_ACCEPTED_ONLY',flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[('slow',2.466),('fast',1.233)]))
