"""Add the established symmetric tonic support and require early gates."""
import sys,json,shutil,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent;source=Path(sys.argv[1]).resolve();out=source.parent/'tonic_0p02_exact';W=out/'weights'
if out.exists():raise SystemExit('Existing evidence retained')
W.mkdir(parents=True)
for p in source.glob('*.json'):shutil.copyfile(p,W/p.name)
from bilateral import assert_tied
p=W/'mn_weights.json';z=json.loads(p.read_text())
reference=json.loads((H/'runs/tonic_abduction_20261030_s1_v1/0p02/weights/mn_weights.json').read_text())
abductor_edges={r['edge'] for r in reference['experimental_tonic_abductor_support']['changes']}
for key in z['edges']:
 pre,post=key.split(' -> ')
 if key in abductor_edges:z['edges'][key]+=.02
assert_tied(z['edges']);z['tonic_abductor_increment']=.02;p.write_text(json.dumps(z,indent=2))
subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(W),'--output-directory',str(out/'air_pair')],stdout=(out/'air_pair.log').open('w'),stderr=subprocess.STDOUT,check=True)
pair=json.loads((out/'air_pair/SUMMARY.json').read_text());print('AIR_PAIR',pair['accepted'],flush=True)
if not pair['accepted']:raise SystemExit(0)
base=json.loads((H/'runs/tonic_abduction_20261030_s1_v1/0p02/standing/task_low.json').read_text());base['weights_directory']=str(W)
for name,cfg in [('standing',base),('ground',{**base,'stage':4,'force_controlled_support':True,'support_fraction':.7,'regional_contacts':True,'cross_contact_gain':6.,'cross_contact_event_amplitude':2.})]:
 p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT,check=True)
 r=json.loads(p.with_suffix('.result.json').read_text());print(name,r['assessment']['accepted'],flush=True)
 if not r['assessment']['accepted']:break
