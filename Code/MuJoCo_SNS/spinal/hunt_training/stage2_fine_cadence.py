from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/stage2_fine_cadence_v1';campaign.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/stage2_coordination_v1/ext0.3_comm0.3_axis1.5.json').read_text())
cases=[(f'comm{comm}_timing{timing}',{**base,'commissural_gain':comm,'timing_scale':timing}) for comm in (.3,.6,1.) for timing in (2.3,2.6)]
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda c:campaign.run(*c),cases))
