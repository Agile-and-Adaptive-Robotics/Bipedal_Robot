"""Air stages use the calibrated air start, not the neutral standing pose."""
from pathlib import Path
import balance_campaign as campaign
H=Path(__file__).resolve().parent;campaign.O=H/'runs/balance_20261029_v3'
for stage in (1,2):campaign.run('curriculum_air_stage_'+str(stage),dict(stage=stage))
