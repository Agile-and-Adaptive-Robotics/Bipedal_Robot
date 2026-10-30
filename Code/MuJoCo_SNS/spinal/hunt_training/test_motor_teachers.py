"""Ablate neural-fit error versus inadequate motor targets in the failed air stage."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
H=Path(__file__).resolve().parent;O=H/'runs/swing_release_20261030_v1/teachers_cached';O.mkdir(exist_ok=True);B.O=O
base=json.loads((H/'runs/fourcell_dynamic_20261029_v1/curriculum/stage2.json').read_text());base['diagnose_swing']=True
cases=[('reference_replay',{'mode':'reference'}),('feedback_teacher',{'mode':'feedback'})]
def run(item):
    name,teacher=item;r=B.run(name,{**base,'motor_teacher':teacher});print(name,json.dumps(r['assessment']) if r else None,flush=True)
with ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,cases))
