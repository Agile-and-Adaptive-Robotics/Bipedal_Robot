"""Bounded prepared-state stress protocol; one fresh physics process at a time."""
import argparse,json,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'pinned_preparation'

def launch(case):
 ledger=json.loads((OUT/'SUMMARY.json').read_text())
 if case in ledger['cases']:return ledger['cases'][case]
 with (ROOT/('prepared_'+case+'.log')).open('w',encoding='utf-8') as log:subprocess.run([sys.executable,str(HERE/'run_pinned_preparation.py'),'--case',case],stdout=log,stderr=subprocess.STDOUT,check=True)
 return json.loads((OUT/'SUMMARY.json').read_text())['cases'][case]

def main(loosen,startup_controls):
 gate=json.loads((OUT/'STATE_MATCH_AUDIT.json').read_text())
 if not gate['nominal_exact_fresh_gate']:raise ValueError('Fresh exact nominal gate not met')
 for case in ('nominal_replay','harsh_live','harsh_replay','harsh_late_live','harsh_late_replay'):
  row=launch(case);print(json.dumps(dict(case=case,accepted=row['assessment']['accepted'],failure=row['failure'])),flush=True)
 subprocess.run([sys.executable,str(HERE/'audit_pinned_preparation.py')],check=True)
 audit=json.loads((OUT/'STATE_MATCH_AUDIT.json').read_text())
 if not audit['checks']['nominal_replay'].get('whole_trace_bit_exact'):raise ValueError('Nominal sensory replay parity gate failed')
 for case in ('harsh_live','harsh_replay','harsh_late_live','harsh_late_replay'):
  if not audit['checks'][case]['release_matches_common_nominal']:raise ValueError('Unmatched initialized state: '+case)
 if startup_controls:
  for case in ('harsh_zero_live','harsh_zero_replay','harsh_cold_live','harsh_cold_replay'):
   row=launch(case);print(json.dumps(dict(case=case,accepted=row['assessment']['accepted'],failure=row['failure'])),flush=True)
 if loosen:
  for case in ('rot0p8','rot0p5','lat0p8','vertical0p8'):
   row=launch(case);print(json.dumps(dict(case=case,accepted=row['assessment']['accepted'])),flush=True)
   if not row['assessment']['accepted']:break
   (OUT/'SELECTED_RUNG.json').write_text(json.dumps(dict(case=case,result=row['result'],rig_scales=json.loads(Path(row['result']).read_text())['rig_scales']),indent=2),encoding='utf-8')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--loosen',action='store_true');p.add_argument('--startup-controls',action='store_true');a=p.parse_args();main(a.loosen,a.startup_controls)
