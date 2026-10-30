"""Check physical cadence transfer and fresh-process robustness, stage 2 only."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,copy,numpy as np
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/stage2_cadence_v1';campaign.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/stage2_coordination_v1/ext0.3_comm0.3_axis1.5.json').read_text())
cases=[('fresh_replay',base)]
rng=np.random.default_rng(291029)
for i in range(4):
    cfg=copy.deepcopy(base)
    for key in ('reciprocal_gain','commissural_gain','axis_feedback','hip_abduction_goal'):cfg[key]*=rng.uniform(.99,1.01)
    for key in cfg['motor_scales']:cfg['motor_scales'][key]*=rng.uniform(.99,1.01)
    cases.append(('basin_'+str(i),cfg))
for comm in (.3,.6):
    for timing in (1.5,2.,3.):cases.append((f'comm{comm}_timing{timing}',{**base,'commissural_gain':comm,'timing_scale':timing}))
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda case:campaign.run(*case),cases))
