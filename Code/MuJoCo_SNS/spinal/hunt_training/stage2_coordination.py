"""Coordination and explicit proprioceptive axis-control probes, stage 2 only."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/stage2_coordination_v1';campaign.O.mkdir(exist_ok=True)
cases=[]
for motor in (.3,.6):
    for comm in (.3,.6):
        for axis in (.5,1.5):
            cfg=dict(stage=2,weights_directory=str(H/'runs/encoder_fixed_v1'),reciprocal_gain=.1,commissural_gain=comm,axis_feedback=axis,hip_abduction_goal=.08,
                motor_scales={'knee_ext_l':motor,'knee_ext_r':motor,'knee_flex_l':2.,'knee_flex_r':2.})
            cases.append((f'ext{motor}_comm{comm}_axis{axis}',cfg))
with ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(lambda case:campaign.run(*case),cases))
