"""Rerun failed neural air stages; never advance on teacher or failed trials."""
import argparse,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
H=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--campaign',type=Path,default=H/'runs/net_torque_pf_20261030_v1')
root=parser.parse_args().campaign.resolve();assert (root/'TRAINING_COMPLETE.json').exists();B.O=root/'air_retry';B.O.mkdir(exist_ok=True)
base=dict(weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
cases=[('stage1',{'stage':1}),('stage2',{'stage':2,'diagnose_swing':True}),('stage2_axis',{'stage':2,'axis_feedback':1.5,'hip_abduction_goal':.08})]
def run(item):
    name,change=item;return name,B.run(name,{**base,**change})
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(run,cases))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=r['configuration'],assessment=r['assessment'],measured=r['measured']) for k,r in rows.items() if r},indent=2))
print('FAILED_STAGE_RECHECK_COMPLETE',json.dumps({k:r['assessment'] for k,r in rows.items() if r}),flush=True)
