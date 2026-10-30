"""Ground-transfer probes after stage-2 replay/perturbations and stage 3 pass."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;source=H/'runs/ground_after_stage2_v1'
if not json.loads((source/'stage3.result.json').read_text())['assessment']['accepted']:raise RuntimeError('Standing gate not satisfied')
base=json.loads((source/'stage4_height-0.02.json').read_text());campaign.O=H/'runs/ground_contact_rules_v1';campaign.O.mkdir(exist_ok=True)
cases=[]
for cross in (.5,2.,4.):
    for lateral in (1.,.01):cases.append((f'cross{cross}_lateral{lateral}',{**base,'cross_contact_gain':cross,'lateral_stiffness_scale':lateral}))
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda c:campaign.run(*c),cases))
