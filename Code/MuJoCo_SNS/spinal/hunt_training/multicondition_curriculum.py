"""Recheck the corrected curriculum after fitting both cadence conditions."""
import json,os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
from ground_contact_check import check
H=Path(__file__).resolve().parent;campaign=os.environ.get('HUNT_CAMPAIGN','fourcell_multicondition_20261029_v1')
if Path(campaign).name!=campaign:raise ValueError('Campaign must be a directory name under hunt_training/runs')
root=H/'runs'/campaign;assert (root/'TRAINING_COMPLETE.json').exists()
B.O=root/'curriculum';B.O.mkdir(exist_ok=True)
base=dict(weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
cases=[('stage1',dict(stage=1)),('stage1_restrict05',dict(stage=1,restriction=.05)),('stage1_restrict20',dict(stage=1,restriction=.2)),
    ('stage2',dict(stage=2)),('stage2_axis',dict(stage=2,axis_feedback=1.5,hip_abduction_goal=.08)),
    ('stage2_axis_halfreflex',dict(stage=2,axis_feedback=1.5,hip_abduction_goal=.08,reflex_scale=.5)),
    ('stage2_restrict05',dict(stage=2,restriction=.05)),('stage2_restrict20',dict(stage=2,restriction=.2))]
def work(item):
    name,change=item;r=B.run(name,{**base,**change})
    if r:print(name,json.dumps(r['assessment']),flush=True)
    return name,r
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(work,cases))
air=[r for name,r in rows.items() if name.startswith('stage2') and r]
valid=[r for r in air if r['assessment']['accepted']]
best=max(valid or air,key=lambda r:r['assessment']['score'])
stand={**best['configuration'],'stage':3,'restriction':.2,'height_offset':-.02,'start_hip_deg':-7.,'start_knee_deg':-10.,'start_ankle_deg':0.,'balance':dict(kp=.5,kd=.2,ka=0.,delay_s=.1,replace_inherited=False)}
rows['stage3']=B.run('stage3',stand)
for name,change in [('stage4_static',{}),('stage4_lowering',dict(harness_lowering_m=-.06)),('stage4_towed',dict(harness_lowering_m=-.06,harness_forward_speed=.04,lateral_stiffness_scale=.001))]:
    r=B.run(name,{**stand,'stage':4,**change});rows[name]=r
    if r:
        r['strict_ground']=check(str(B.O/name)+'.assistance.npz',r['measured']['body_weight_N']);(B.O/(name+'.strict.json')).write_text(json.dumps(r['strict_ground'],indent=2));print(name,json.dumps(r['strict_ground']),flush=True)
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured'],strict_ground=v.get('strict_ground')) for k,v in rows.items() if v},indent=2))
