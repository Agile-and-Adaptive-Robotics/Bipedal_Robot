"""Explicit physical-calibration probes; changing timing gains requires refitting."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
H=Path(__file__).resolve().parent;root=H/'runs/fourcell_hunt_20261029_v1';B.O=root/'phase_probes';B.O.mkdir(exist_ok=True)
base=json.loads((root/'curriculum/stage2_scaled_probe.json').read_text())
cases=[('comm_'+str(g)+'_reflex_'+str(f),{**base,'commissural_gain':g,'reflex_scale':f}) for g in (.3,.6,1.2,2.4) for f in (.5,1.)]
cases += [('stage1_restriction_'+str(r),{**base,'stage':1,'restriction':r}) for r in (0.,.05,.2)]
def work(item):
    name,cfg=item;r=B.run(name,cfg)
    if r:print(name,json.dumps(r['assessment']),flush=True)
    return name,r
with ThreadPoolExecutor(max_workers=3) as pool:rows=dict(pool.map(work,cases))
(B.O/'summary.json').write_text(json.dumps({k:dict(configuration=v['configuration'],assessment=v['assessment'],measured=v['measured']) for k,v in rows.items() if v},indent=2))
