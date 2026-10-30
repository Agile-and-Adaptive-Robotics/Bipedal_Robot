"""Diagnostic gradual lowering / constant-speed harness assistance, no stride clock."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
from ground_contact_check import check
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1';B.O=root/'ground_assistance';B.O.mkdir(exist_ok=True)
base=json.loads((root/'curriculum/stage3.json').read_text());base.update(stage=4,reflex_scale=.5,commissural_gain=.3,restriction=.2)
cases=[('lower_'+str(lower)+'_speed_'+str(speed),{**base,'harness_lowering_m':lower,'harness_forward_speed':speed}) for lower in (-.04,-.08,-.12) for speed in (0.,.04,.08)]
def work(item):
    name,cfg=item;r=B.run(name,cfg)
    if r:
        strict=check(str(B.O/name)+'.assistance.npz',r['measured']['body_weight_N']);r['strict_ground']=strict
        (B.O/(name+'.strict.json')).write_text(json.dumps(strict,indent=2));print(name,'strict',json.dumps(strict),flush=True)
    return name,r
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(work,cases))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured'],strict_ground=v['strict_ground']) for k,v in rows.items() if v},indent=2))
