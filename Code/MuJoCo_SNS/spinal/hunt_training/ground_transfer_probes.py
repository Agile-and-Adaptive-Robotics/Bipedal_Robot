"""Test finite walk-start initialization and approach ground from above."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/ground_transfer_v1';campaign.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/ground_after_stage2_v1/stage4_height-0.02.json').read_text())
cases=[(f'height{height}_kick{kick}',{**base,'height_offset':height,'walk_kick':kick}) for height in (0.,.02,.04) for kick in (5.,12.)]
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda c:campaign.run(*c),cases))
