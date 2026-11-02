"""Fresh/long/sensitivity checks of the explicit all-path-preserving mapping."""
import argparse,json,sys,subprocess
from pathlib import Path
import numpy as np
import gait_standing_family as F
from gait_standing_objective import score
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'group_diagnosis';CHECK=OUT/'normalized_checks'

def run(case):
 CHECK.mkdir(exist_ok=True);source=json.loads((OUT/'w2lvar_fanin_mean_v1.assessment.json').read_text());params=json.loads((ROOT/'rustuna_w2lvar_v1/BEST.json').read_text())['trial']['params'];scales=dict(source['scales']);duration=20. if case=='hold20' else 8.
 if case.startswith('basin'):
  rng=np.random.default_rng(100040+int(case[5:]));scales={k:v*(1+rng.choice([-1.,1.])*.01) for k,v in scales.items()};params['standing_bias_scale']*=1+rng.choice([-1.,1.])*.01
 r=F.run('w2lvar',CHECK,duration=duration,tag=case+'_v1',gain=params['pf_mn_gain'],ib_group_gain=params['ib_group_gain'],standing_bias_scale=params['standing_bias_scale'],pathway_group_scales=scales,capture_neural_inputs=True);prefix=CHECK/('w2lvar_'+case+'_v1');report=dict(case=case,parameters=params,scales=scales,assessment=score(r,prefix)[1],result=str(prefix.with_suffix('.result.json')))
 prefix.with_suffix('.assessment.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report),flush=True)

def sequence():
 CHECK.mkdir(exist_ok=True)
 for case in ('fresh','hold20')+tuple('basin'+str(i) for i in range(8)):
  target=CHECK/('w2lvar_'+case+'_v1.assessment.json')
  if not target.exists():
   with (CHECK/(case+'.log')).open('w',encoding='utf-8') as log:finished=subprocess.run([sys.executable,str(Path(__file__)),'--case',case],stdout=log,stderr=subprocess.STDOUT)
   if finished.returncode:raise RuntimeError(case)
  row=json.loads(target.read_text());print(json.dumps(dict(case=case,assessment=row['assessment'])),flush=True)
  if not row['assessment']['accepted']:break

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case');a=p.parse_args();run(a.case) if a.case else sequence()
