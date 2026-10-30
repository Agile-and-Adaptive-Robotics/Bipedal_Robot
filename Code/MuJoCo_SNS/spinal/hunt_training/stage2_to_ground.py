"""Advance only after fresh, perturbed stage-2 movement/frontal gates pass."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;prior=H/'runs/stage2_cadence_v1'
required=[prior/'fresh_replay.result.json',*[prior/f'basin_{i}.result.json' for i in range(4)]]
if not all(p.exists() and json.loads(p.read_text())['assessment']['accepted'] for p in required):raise RuntimeError('Stage-2 robustness screen not satisfied')
base=json.loads((prior/'fresh_replay.json').read_text())
campaign.O=H/'runs/ground_after_stage2_v1';campaign.O.mkdir(exist_ok=True)
balance=json.loads((H/'runs/balance_20261029_v3/winner.json').read_text())['configuration']['balance']
standing={**base,'stage':3,'restriction':1.,'height_offset':-.02,'start_hip_deg':-7.,'start_knee_deg':-10.,'start_ankle_deg':0.,'balance':balance}
r=campaign.run('stage3',standing)
if not r or not r['assessment']['accepted']:raise SystemExit('Stage 3 fails with revised controller: no ground stepping progression')
cases=[(f'stage4_height{h}',{**standing,'stage':4,'height_offset':h}) for h in (-.02,-.04,-.06)]
with ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(lambda c:campaign.run(*c),cases))
passed=[(c,r) for c,r in zip(cases,results) if r and r['assessment']['accepted']]
if passed:
    case,r=max(passed,key=lambda x:x[1]['assessment']['score'])
    replay=campaign.run('stage4_fresh',case[1])
    if replay and replay['assessment']['accepted']:
        for gain in (.05,.15):campaign.run('stage5_contact'+str(gain),{**case[1],'stage':5,'contact_pf_scale':gain})
else:print('No stage-4 candidate passed; no stage-5 progression',flush=True)
