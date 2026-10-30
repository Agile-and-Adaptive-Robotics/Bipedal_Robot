"""Compare currents and forces in previously failed air and ground cases."""
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as B
H=Path(__file__).resolve().parent;prior=H/'runs/fourcell_dynamic_20261029_v1';B.O=H/'runs/swing_release_20261030_v1';B.O.mkdir(exist_ok=True)
def work(item):
    name,filename=item;cfg=json.loads((prior/'curriculum'/filename).read_text());cfg['diagnose_swing']=True
    return name,B.run(name,cfg)
with ThreadPoolExecutor(max_workers=2) as pool:results=dict(pool.map(work,[('air','stage2.json'),('ground','stage4_static.json')]))
print('DIAGNOSTICS_DONE',flush=True)
