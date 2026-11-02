"""Same declared standing search/gates/budget for the two remaining families."""
import json
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
source=json.loads((HERE/'GAIT_SYN6_STANDING_MANIFEST.json').read_text())
for family in ('s3k','w2lvar'):
 cfg=dict(source);cfg['family']=family;cfg['output_dir']=str(ROOT/('rustuna_'+family+'_v1'))
 cfg['scope']='Same shared bilateral biological-Hill standing mechanics, rig, corrected pose, objective gates, search bounds, seed, twelve-trial/600s budget, and initial candidates as syn6. Existing family baseline/PF0 are reused. Stage-specific recruitment calibration; no walking or unchanged literature-replication claim.'
 path=ROOT/('GAIT_'+family.upper()+'_STANDING_MANIFEST.json')
 if path.exists() and json.loads(path.read_text())!=cfg:raise FileExistsError(path)
 path.write_text(json.dumps(cfg,indent=2),encoding='utf-8');print(path)
