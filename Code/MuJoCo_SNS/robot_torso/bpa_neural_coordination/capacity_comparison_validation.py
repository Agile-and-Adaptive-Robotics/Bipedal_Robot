"""Restartable frozen selection, fresh30s parity, then untouched validation."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
import json,sys,subprocess,argparse,hashlib
from pathlib import Path
import numpy as np
from capacity_comparison_adapter import registration,evaluate_case
from register_capacity_comparison import DEST

def select(mode):
 p=DEST/mode/'SELECTED.json'
 if p.exists():return json.loads(p.read_text())
 lines=DEST/mode/'rustuna/TRIALS.jsonl';rows=[json.loads(line) for line in lines.read_text().splitlines()];good=[row for row in rows if row['diagnostics']['all_training_accepted']]
 if len(rows)!=4:raise ValueError('Registered trial budget must finish before selection')
 if not good:return None
 winner=min(good,key=lambda row:row['value']);result=dict(params=winner['params'],training_trial=winner['number'],training_value=winner['value'],selection='Minimum registered training objective among alltraining passers only. Frozen before fresh30s stability and before new validation.',training_file_sha256=hashlib.sha256(lines.read_bytes()).hexdigest())
 p.write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def child(mode,name):
 folder=DEST/mode/name;summary=folder/'CASE_SUMMARY.json'
 if summary.exists():return json.loads(summary.read_text())
 folder.mkdir(parents=True,exist_ok=True)
 with (folder/'RUN.log').open('w',encoding='utf-8') as log:proc=subprocess.run([sys.executable,str(Path(__file__)),'--mode',mode,'--case',name],stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
 if proc.returncode:raise RuntimeError('Failed case; inspect '+str(folder/'RUN.log'))
 return json.loads(summary.read_text())

def run(mode):
 cfg=registration();selection=select(mode);status=dict(mode=mode,alltraining_candidate=selection is not None,new_validation_started=False)
 if selection is None:
  status['scope']='No candidate passes alltraining. No freshfinalist or untouched validation run.'
 else:
  a=child(mode,'fresh_nominal_a');b=child(mode,'fresh_nominal_b');status['nominal30s']=[a,b];checks={}
  if a['accepted'] and b['accepted']:
   for p in sorted((DEST/mode/'fresh_nominal_a').glob('w2lvar_case*.npz')):
    q=DEST/mode/'fresh_nominal_b'/p.name
    with np.load(p) as x,np.load(q) as y:
     assert set(x.files)==set(y.files);checks[p.name]={key:bool(np.array_equal(x[key],y[key])) for key in x.files}
   exact=all(v for group in checks.values() for v in group.values());status['fresh_complete_array_parity']=checks;status['fresh_exact']=exact
   if exact:
    status['new_validation_started']=True;status['validation']=[dict(condition=f['name'],**child(mode,f['name'])) for f in cfg['validation']]
  status['scope']='Parameters remain frozen regardless of validation outcome. No postvalidation mode or parameter reselection; compare modes with all failures shown.'
 (DEST/mode/'VALIDATION_STATUS.json').write_text(json.dumps(status,indent=2),encoding='utf-8');print(json.dumps(status,indent=2))

def main():
 p=argparse.ArgumentParser();p.add_argument('--mode',choices=['pressure','force'],required=True);p.add_argument('--case');a=p.parse_args()
 if a.case:
  cfg=registration();selection=json.loads((DEST/a.mode/'SELECTED.json').read_text());fresh=a.case.startswith('fresh_nominal_');spec=cfg['training'][0] if fresh else next(f for f in cfg['validation'] if f['name']==a.case);r=evaluate_case(selection['params'],spec,DEST/a.mode/a.case,duration=30. if fresh else 20.,full=fresh);print(json.dumps(r,indent=2))
 else:run(a.mode)

if __name__=='__main__':main()
