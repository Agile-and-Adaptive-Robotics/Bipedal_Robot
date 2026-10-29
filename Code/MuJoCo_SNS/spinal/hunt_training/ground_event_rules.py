"""Physical contact-onset feedback, not a periodic external phase command."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import json
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/ground_event_rules_v1';campaign.O.mkdir(exist_ok=True)
base=json.loads((H/'runs/ground_after_stage2_v1/stage4_height-0.02.json').read_text())
cases=[(f'cross{cross}_onset{onset}',{**base,'cross_contact_gain':cross,'contact_onset':onset}) for cross in (.5,2.) for onset in (1.,3.,6.)]
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda c:campaign.run(*c),cases))
