"""Isolate RG timing feedback from MN reflex feedback before changing weights."""
from pathlib import Path
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/stage2_diagnosis_v1';campaign.O.mkdir(exist_ok=True)
for name,timing,reflex in [('both_off',False,0.),('timing_only',True,0.),('reflex_only',False,1.),('both_on',True,1.)]:
    campaign.run(name,dict(stage=2,timing_afferents=timing,reflex_scale=reflex))
