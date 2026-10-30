"""Preserve original weights; correct encoder parity and refit only Hunt D."""
import shutil,json
from pathlib import Path
import motor_fit as M
import controller as C
H=Path(__file__).resolve().parent;O=H/'runs/encoder_fixed_v1';O.mkdir(exist_ok=True)
for name in ('pf_weights.json','mn_weights.json','reflex_weights.json'):
    if not (O/name).exists():shutil.copyfile(H/'runs/fitted'/name,O/name)
real=C.build
def build(*args,**kwargs):
    options=kwargs.pop('options',{}) or {};options['weights_directory']=str(O)
    return real(*args,options=options,**kwargs)
C.build=build;M.O=O
M.replay('REFLEX');M.validate()
import balance_campaign as campaign
campaign.O=H/'runs/stage2_encoder_fixed_v1';campaign.O.mkdir(exist_ok=True)
for name,timing,scale in [('full',True,1.),('timing_only',True,0.),('reflex_only',False,1.),('half_reflex',True,.5)]:
    campaign.run(name,dict(stage=2,timing_afferents=timing,reflex_scale=scale,weights_directory=str(O)))
