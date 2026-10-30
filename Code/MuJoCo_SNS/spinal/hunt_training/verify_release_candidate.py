"""Use one controller configuration for both early curriculum stages."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json,balance_campaign as B
H=Path(__file__).resolve().parent;root=H/'runs/native_joint_pf_mono_20261030_v1'
B.O=root/'release_candidate';B.O.mkdir(exist_ok=True)
base=dict(weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True,motor_scales={'hip_abd_r':8.,'hip_abd_l':8.})
def run(stage):return str(stage),B.run('stage'+str(stage),{**base,'stage':stage,'diagnose_swing':True})
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,(1,2)))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured']) for k,v in rows.items() if v},indent=2))
print('MATCHED_CANDIDATE_COMPLETE',json.dumps({k:v['assessment'] for k,v in rows.items() if v}),flush=True)
