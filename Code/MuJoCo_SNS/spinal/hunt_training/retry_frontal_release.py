"""Check abductor feedforward strength against the unchanged stage-1 gates."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,json,balance_campaign as B
H=Path(__file__).resolve().parent;root=H/'runs/native_joint_pf_mono_20261030_v1'
B.O=root/'frontal_retry';B.O.mkdir(exist_ok=True)
base=dict(stage=1,weights_directory=str(root/'weights'),biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True)
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--scales',type=float,nargs='+',default=[1.1,1.25,1.5]);args=parser.parse_args()
def run(scale):
    return str(scale),B.run('abductors_'+str(scale).replace('.','p'),{**base,'motor_scales':{'hip_abd_r':scale,'hip_abd_l':scale}})
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(run,args.scales))
summary=json.loads((B.O/'summary.json').read_text()) if (B.O/'summary.json').exists() else {}
summary.update({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured']) for k,v in rows.items() if v})
(B.O/'summary.json').write_text(json.dumps(summary,indent=2))
print('FRONTAL_RECHECK_COMPLETE',json.dumps({k:v['assessment'] for k,v in rows.items() if v}),flush=True)
