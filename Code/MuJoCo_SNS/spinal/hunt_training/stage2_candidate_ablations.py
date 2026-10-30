"""Identify which changes are necessary in the passing slow stage-2 candidate."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import copy,json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/stage2_candidate_ablations_v1';campaign.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/stage2_cadence_v1/fresh_replay.json').read_text())
cases=[]
for name,key,value in [('no_reciprocal','reciprocal_gain',0.),('no_commissural','commissural_gain',0.),('no_axis_feedback','axis_feedback',0.)]:cases.append((name,{**base,key:value}))
cfg=copy.deepcopy(base);cfg.pop('motor_scales');cases.append(('original_motor_gains',cfg))
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda c:campaign.run(*c),cases))
