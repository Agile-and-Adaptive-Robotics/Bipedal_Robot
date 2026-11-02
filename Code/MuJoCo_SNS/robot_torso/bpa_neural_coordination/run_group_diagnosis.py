"""One-path/group interventions; frozen family standing gains, actual input capture."""
import argparse,json
from pathlib import Path
import gait_standing_family as F
from gait_standing_objective import score
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis'

def main(case,measured=False):
 family=case if case in ('s3k','syn6','w2lvar') else 'w2lvar';best=json.loads((ROOT/('rustuna_'+family+'_v1/BEST.json')).read_text());params=best['trial']['params'];scales={}
 if case not in ('s3k','syn6','w2lvar'):scales={'IIIN:'+case:0.}
 tag=case+('_measured_v2' if measured else '_v1')
 r=F.run(family,OUT,tag=tag,gain=params['pf_mn_gain'],ib_group_gain=params['ib_group_gain'],standing_bias_scale=params['standing_bias_scale'],pathway_group_scales=scales,capture_neural_inputs=True)
 prefix=OUT/(family+'_'+tag);report=dict(case=case,assessment=score(r,prefix)[1],result=str(prefix.with_suffix('.result.json')),scope='Matched mechanics/rig, frozen selected family standing gains; intervention affects only the named presynaptic pathway onto a primary muscle group, bilaterally. No permanent source edits.')
 prefix.with_suffix('.assessment.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case',required=True);p.add_argument('--measured',action='store_true');args=p.parse_args();main(args.case,args.measured)
