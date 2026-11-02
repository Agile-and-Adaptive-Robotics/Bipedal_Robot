"""Two distinct ±1% corners: both fitted connection multipliers are exactly zero."""
import json,subprocess,sys
from pathlib import Path
from gait_standing_objective import score
HERE=Path(__file__).resolve().parent;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');BEST=ROOT/'rustuna_s3k_v1/BEST.json'

def main():
 params=json.loads(BEST.read_text())['trial']['params'];assert params['pf_mn_gain']==0 and params['ib_group_gain']==0
 out=ROOT/'rustuna_s3k_v1/basin';out.mkdir(exist_ok=True);rows=[]
 for sign in (-1,1):
  tag='minus1pct_v1' if sign<0 else 'plus1pct_v1';p=dict(params);p['standing_bias_scale']*=1+sign*.01;config=out/(tag+'.json')
  if config.exists() and json.loads(config.read_text())!=p:raise FileExistsError(config)
  config.write_text(json.dumps(p,indent=2),encoding='utf-8');path=out/('s3k_'+tag+'.result.json')
  if not path.exists():
   with (out/(tag+'.log')).open('w',encoding='utf-8') as log:r=subprocess.run([sys.executable,str(HERE/'gait_standing_family.py'),'--family','s3k','--out',str(out),'--tag',tag,'--parameters-json',str(config)],stdout=log,stderr=subprocess.STDOUT)
   if r.returncode:raise RuntimeError(tag)
  result=json.loads(path.read_text());a=score(result,path.with_suffix('').with_suffix(''))[1];rows.append(dict(sign=sign,parameters=p,assessment=a,result=str(path)));print(json.dumps(rows[-1]),flush=True)
 report=dict(machine='easteregg2',unique_corners=2,nominal_three_parameter_corner_count=8,both_zero_parameters_unchanged_under_relative_perturbation=True,all_pass=all(r['assessment']['accepted'] for r in rows),cases=rows,scope='Closed-loop8s standing at the two DISTINCT ±1% fitted-parameter corners. PF and aggregate-Ib fitted multipliers are zero, so all eight nominal three-dimensional corners collapse to these two. This does not test absolute positive activation of disconnected branches, other sensor gains, pose/load uncertainty, or BPA actuation.')
 (out/'BASIN_CHECK.json').write_text(json.dumps(report,indent=2),encoding='utf-8')

if __name__=='__main__':main()
