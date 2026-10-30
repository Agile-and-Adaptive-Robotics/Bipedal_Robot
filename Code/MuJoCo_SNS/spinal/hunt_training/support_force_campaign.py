"""Actual prescribed vertical support, distinct from stiffness scaling."""
import json
from pathlib import Path
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/balance_20261029_force_support';campaign.O.mkdir(exist_ok=True)
source=H/'runs/balance_20261029_v3/winner.json'
if not source.exists():raise RuntimeError('Finish additive gain fitting first')
base=json.loads(source.read_text())['configuration']
for stage in (3,4):
    for fraction in (.6,.4):
        cfg={**base,'stage':stage,'force_controlled_support':True,'support_fraction':fraction}
        if stage==4:cfg['motor_scales']={'knee_ext_l':.5,'knee_ext_r':.5,'knee_flex_l':2.,'knee_flex_r':2.}
        campaign.run(f's{stage}_support{fraction}',cfg)
