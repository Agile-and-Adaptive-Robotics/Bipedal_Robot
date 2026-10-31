"""Test shared RF recruitment reduction after restoring its source force curve."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'rectus_swing_recruitment_v1';out.mkdir()
base=json.loads((root/'ground_final_checks_v1/no_comm_support_0.5.json').read_text());cases={}
for scale in (.5,.25):
 folder=out/('rf_'+str(scale));weights=folder/'weights';shutil.copytree(base['weights_directory'],weights)
 p=weights/'mn_weights.json';spec=json.loads(p.read_text());changed=[]
 for key in spec['edges']:
  if key.split(' -> ')[1] in ('MN_rect_fem_r','MN_rect_fem_l'):spec['edges'][key]*=scale;changed.append(key)
 assert changed and len(changed)%2==0
 from bilateral import assert_tied
 assert_tied(spec['edges']);p.write_text(json.dumps(spec,indent=2));cfg={**base,'weights_directory':str(weights)};cases['rf_'+str(scale)]=cfg
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
rows={}
with ThreadPoolExecutor(max_workers=2) as pool:
 for n,r in pool.map(run,cases.items()):rows[n]=r;(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2));print(n,r.get('assessment'),r.get('foot_clearance'),flush=True)
