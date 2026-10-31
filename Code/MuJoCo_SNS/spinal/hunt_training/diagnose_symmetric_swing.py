"""Restore canonical ankle PF rhythm and test symmetric swing recruitment.

Diagnostic candidates must revalidate early stages before any promotion.
No external joint trajectory, motor clamp, or stride clock is used.
"""
import sys,json,shutil,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from bilateral import assert_tied
H=Path(__file__).resolve().parent;source=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/weights';out=source.parent/'swing_diagnostics_v1';out.mkdir()
base=json.loads((H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact/ground.json').read_text());plans={}
for increment in (0.,.15,.3):
 name='swing_'+str(increment).replace('.','p');W=out/name/'weights';W.mkdir(parents=True)
 for p in source.glob('*.json'):shutil.copyfile(p,W/p.name)
 p=W/'circuit_options.json';z=json.loads(p.read_text());z['direct_ankle_pf_gain']=2.4;p.write_text(json.dumps(z,indent=2))
 p=W/'pf_weights.json';z=json.loads(p.read_text())
 for edge in z['edges']:
  pre,post=edge.split(' -> ')
  if pre.startswith('HUNT_BANK') and post.startswith('PF_FUNC_ANKLE_Z_'):z['edges'][edge]=0.
 assert_tied(z['edges']);p.write_text(json.dumps(z,indent=2))
 p=W/'mn_weights.json';z=json.loads(p.read_text())
 for s in ('r','l'):
  z['edges']['PF_FUNC_ANKLE_Z_N_'+s+' -> MN_tib_ant_'+s]=.3
  for pre,post in [('PF_KNEE-F_','MN_bifemsh_'),('PF_HIP-F_','MN_iliacus_'),('PF_HIP-F_','MN_psoas_')]:
   key=pre+s+' -> '+post+s;assert key in z['edges'];z['edges'][key]+=increment
 assert_tied(z['edges']);p.write_text(json.dumps(z,indent=2))
 plans[name]={**base,'weights_directory':str(W),'diagnostic_scope':'Symmetric swing-expression test; early stages require revalidation.'}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,plans.items()))
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2))
for n,r in rows.items():print(n,r.get('assessment',{}).get('accepted'),flush=True)
