"""Small isolated wiring/gain probes; no downstream stage is run here."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/stage2_rules_v1';campaign.O.mkdir(exist_ok=True)
weights=str(H/'runs/encoder_fixed_v1')
cases=[]
for reciprocal,comm in ((.1,0.),(.5,0.),(.1,.05),(.5,.15)):
    for timing in (1.,.25):
        cfg=dict(stage=2,weights_directory=weights,reciprocal_gain=reciprocal,commissural_gain=comm,timing_scale=timing)
        cases.append((f'recip{reciprocal}_comm{comm}_timing{timing}',cfg))
for comm in (0.,.1):
    for ext in (.3,.6):cases.append((f'motor_ext{ext}_comm{comm}',dict(stage=2,weights_directory=weights,commissural_gain=comm,reciprocal_gain=.1,motor_scales={'knee_ext_l':ext,'knee_ext_r':ext,'knee_flex_l':2.,'knee_flex_r':2.})))
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda case:campaign.run(*case),cases))
