"""Optional four-cell HAM/RECT PF pairing; original candidate remains intact."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/biarticular_pair_v1';campaign.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/stage2_cadence_v1/fresh_replay.json').read_text())
cases=[(f'pair_gain{gain}',{**base,'biarticular_pair_gain':gain}) for gain in (.1,.3,.6)]
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda c:campaign.run(*c),cases))
