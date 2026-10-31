"""Stage-4 diagnostics, not curriculum promotion past a failed air gate."""
import sys,json,shutil,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from bilateral import assert_tied
H=Path(__file__).resolve().parent;source=H/'runs/symmetric_loadedstress_20261030_s1_v1/tonic_0p02_exact/weights';out=source.parent/'ground_diagnostic'
out.mkdir();base=json.loads((H/'runs/tonic_abduction_20261030_s1_v1/0p02/proximal_stance_v1/heel_event.json').read_text())
plans={}
for variant in ('shared_only','shared_proximal'):
 folder=out/variant;W=folder/'weights';W.mkdir(parents=True)
 for p in source.glob('*.json'):shutil.copyfile(p,W/p.name)
 if variant=='shared_proximal':
  p=W/'reflex_weights.json';z=json.loads(p.read_text())
  for key in z['edges']:
   if key.startswith(('IBEXC_hip_ext_','IBEXC_knee_ext_')):z['edges'][key]+=.3
  assert_tied(z['edges']);p.write_text(json.dumps(z,indent=2))
 plans[variant]={**base,'weights_directory':str(W),'diagnostic_scope':'Not curriculum promotion; slow air stage remains failed.'}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,plans.items()))
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2))
for n,r in rows.items():print(n,r.get('assessment',{}).get('accepted'),flush=True)
