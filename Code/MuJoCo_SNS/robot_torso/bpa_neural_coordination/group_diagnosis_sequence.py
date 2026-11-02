import sys,subprocess,json,argparse
from pathlib import Path
from run_group_diagnosis import ROOT,OUT
HERE=Path(__file__).resolve().parent;OUT.mkdir(exist_ok=True)
p=argparse.ArgumentParser();p.add_argument('--ablations',action='store_true');p.add_argument('--measured',action='store_true');args=p.parse_args()
for case in (('ankle_df','knee_ext','hip_ext','ALL') if args.ablations else ('s3k','syn6','w2lvar')):
 family=case if case in ('s3k','syn6','w2lvar') else 'w2lvar';p=OUT/(family+'_'+case+('_measured_v2' if args.measured else '_v1')+'.assessment.json')
 if p.exists():continue
 with (OUT/(case+('_measured' if args.measured else '')+'.log')).open('w',encoding='utf-8') as log:r=subprocess.run([sys.executable,str(HERE/'run_group_diagnosis.py'),'--case',case]+(['--measured'] if args.measured else []),stdout=log,stderr=subprocess.STDOUT)
 if r.returncode:raise RuntimeError(case)
 print(p,flush=True)
