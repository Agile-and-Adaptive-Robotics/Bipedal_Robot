"""Record paired PF and selected MN activity while isolating ground feedback."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
from ground_contact_check import check
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1';B.O=root/'ground_isolation';B.O.mkdir(exist_ok=True)
base=json.loads((root/'stance_load/ext_1.0_flex_1.0_reflex_0.5.json').read_text())
base['trace_cells']=[prefix+side for side in ('r','l') for prefix in ('PF_BI_JOINT_HAM_','PF_BI_JOINT_RECT_','PF_IN_BI_PAIR_HAM_','PF_IN_BI_PAIR_RECT_','MN_rect_fem_','MN_bifemlh_','MN_bifemsh_','MN_vas_lat_')]
cases=[('full_feedback',{}),('no_motor_reflex',dict(reflex_scale=0.)),('no_timing_or_motor_reflex',dict(reflex_scale=0.,timing_afferents=False))]
def work(item):
    name,change=item;r=B.run(name,{**base,**change})
    if r:
        r['strict_ground']=check(str(B.O/name)+'.assistance.npz',r['measured']['body_weight_N']);(B.O/(name+'.strict.json')).write_text(json.dumps(r['strict_ground'],indent=2));print(name,json.dumps(r['strict_ground']),flush=True)
    return name,r
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(work,cases))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured'],strict_ground=v['strict_ground']) for k,v in rows.items() if v},indent=2))
