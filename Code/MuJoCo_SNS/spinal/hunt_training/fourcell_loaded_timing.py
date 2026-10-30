"""Loaded timing ablation with intact stance output; diagnostics, not a refit."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
from ground_contact_check import check
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1';B.O=root/'loaded_timing';B.O.mkdir(exist_ok=True)
base=json.loads((root/'curriculum/stage3.json').read_text());base.update(stage=4,restriction=.2,harness_lowering_m=-.06,motor_scales={},reflex_scale=.5)
cases=[('timing_'+str(timing)+'_lateral_'+str(lateral),{**base,'timing_scale':timing,'lateral_stiffness_scale':lateral,'harness_forward_speed':.04}) for timing in (0.,.25,1.) for lateral in (.001,1.)]
def work(item):
    name,cfg=item;r=B.run(name,cfg)
    if r:
        strict=check(str(B.O/name)+'.assistance.npz',r['measured']['body_weight_N']);r['strict_ground']=strict
        (B.O/(name+'.strict.json')).write_text(json.dumps(strict,indent=2));print(name,json.dumps(strict),flush=True)
    return name,r
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(work,cases))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured'],strict_ground=v['strict_ground']) for k,v in rows.items() if v},indent=2))
