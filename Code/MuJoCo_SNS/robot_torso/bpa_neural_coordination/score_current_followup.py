"""Keep the one-dimensional current-based diagnostic separate from equal-budget TPE."""
import json
from pathlib import Path
from gait_standing_objective import score
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');p=ROOT/'w2lvar_current_followup/w2lvar_current_balance0p25_v1';r=json.loads(p.with_suffix('.result.json').read_text());value,a=score(r,p)
report=dict(machine='easteregg2',shared_pf_mn_gain=.25,ib_group_gain=0.,standing_bias_scale=.8333333333333334,result=str(p.with_suffix('.result.json')),value=value,assessment=a,upward_rig_bw=r['mean_upward_rig_support_bw'],downward_rig_bw=r['mean_downward_rig_force_bw'],scope='One separate current-balance diagnostic after the equal-budget12trial comparison, not part of that budget. Predicted gain0.25 tests residual IIIN inhibition versus lost PF excitation with other gains fixed. No source rewiring. Fresh repeat required before adoption.')
(p.parent/'CURRENT_BALANCE_DIAGNOSTIC.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))
