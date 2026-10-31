"""Fresh repeated ground winner plus matched-controller early-stage checks."""
import sys,json,subprocess,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
H=Path(__file__).resolve().parent;root=H/'runs/symmetric_bothcadence_20261030_s1_v2/tonic_0p02_exact';out=root/'ground_winner_validation_v1';out.mkdir()
base=json.loads((root/'height_tonic_v1/k5000_i12.json').read_text());(out/'winning_configuration.json').write_text(json.dumps(base,indent=2))
air={k:v for k,v in base.items() if k not in ('stage','support_fraction','force_controlled_support','vertical_stiffness_scale','vertical_damping_multiplier','harness_lowering_m','harness_lowering_start_s','harness_lowering_duration_s','height_offset','restriction','start_hip_deg','start_knee_deg','start_ankle_deg','diagnostic_scope')}
air.update(force_controlled_support=False,vertical_stiffness_scale=1.,vertical_damping_multiplier=1.)
(out/'shared_controller.json').write_text(json.dumps(air,indent=2))
with (out/'air_pair.log').open('w') as log:subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',base['weights_directory'],'--output-directory',str(out/'air_pair'),'--config',str(out/'shared_controller.json')],stdout=log,stderr=subprocess.STDOUT,check=True)
print('AIR_PAIR',json.loads((out/'air_pair/SUMMARY.json').read_text())['accepted'],flush=True)
standing={**base,'stage':3,'support_fraction':0.,'force_controlled_support':False,'vertical_stiffness_scale':1.,'vertical_damping_multiplier':1.,'harness_lowering_m':0.}
def run(item):
 name,cfg=item;p=out/(name+'.json');p.write_text(json.dumps(cfg,indent=2))
 with p.with_suffix('.log').open('w') as log:c=subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT).returncode
 return name,json.loads(p.with_suffix('.result.json').read_text()) if c==0 else dict(exit_code=c)
with ThreadPoolExecutor(max_workers=2) as pool:rows=dict(pool.map(run,[('standing',standing),('ground_repeat',base)]))
(out/'SUMMARY.json').write_text(json.dumps(rows,indent=2))
for n,r in rows.items():print(n,r.get('assessment',{}).get('accepted'),flush=True)
for name in ('controller.py','physical_eval.py','contact_events.py','regional_contacts.py','bilateral.py'):
 shutil.copyfile(H/name,out/name)
