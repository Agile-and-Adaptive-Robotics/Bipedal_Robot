"""Check whether air-oriented output scaling suppresses ground load support."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
from ground_contact_check import check
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1';B.O=root/'stance_load';B.O.mkdir(exist_ok=True)
base=json.loads((root/'curriculum/stage3.json').read_text());base.update(stage=4,restriction=.2,harness_lowering_m=-.06)
cases=[]
for ext,flex in ((1.,1.),(1.,2.),(.6,1.)):
    for reflex in (.25,.5,1.):
        name=f'ext_{ext}_flex_{flex}_reflex_{reflex}'
        cfg={**base,'reflex_scale':reflex,'motor_scales':{'knee_ext_r':ext,'knee_ext_l':ext,'knee_flex_r':flex,'knee_flex_l':flex}}
        cases.append((name,cfg))
def work(item):
    name,cfg=item;r=B.run(name,cfg)
    if r:
        strict=check(str(B.O/name)+'.assistance.npz',r['measured']['body_weight_N']);r['strict_ground']=strict
        (B.O/(name+'.strict.json')).write_text(json.dumps(strict,indent=2));print(name,json.dumps(strict),flush=True)
    return name,r
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(work,cases))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured'],strict_ground=v['strict_ground']) for k,v in rows.items() if v},indent=2))
