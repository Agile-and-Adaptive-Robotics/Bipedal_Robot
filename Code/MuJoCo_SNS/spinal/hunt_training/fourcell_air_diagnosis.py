"""Separate frontal mechanics, sensory refinement, and nominal timing failures."""
import json
from pathlib import Path
import balance_campaign as B
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1';B.O=root/'air_diagnosis';B.O.mkdir(exist_ok=True)
base=dict(weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True,restriction=1.)
rows={}
cases=[('stage1_sagittal',dict(stage=1)),('stage2_sagittal',dict(stage=2)),('stage2_reflex_quarter',dict(stage=2,reflex_scale=.25)),('stage2_reflex_off',dict(stage=2,reflex_scale=0.)),('stage2_fixed_reciprocal',dict(stage=2,learn_reciprocal=False,reciprocal_gain=.1))]
for name,change in cases:
    r=B.run(name,{**base,**change})
    if r:rows[name]=r;print(name,json.dumps(r['assessment']),flush=True)
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured']) for k,v in rows.items()},indent=2))
