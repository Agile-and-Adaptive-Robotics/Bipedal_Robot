"""Bounded common motor-expression gains; no left/right compensation."""
import sys,json,subprocess
from pathlib import Path
H=Path(__file__).resolve().parent;source=H/'runs/symmetric_loadedstress_20261030_s1_v1/tonic_0p02_exact/weights';out=source.parent/'motor_scale'
out.mkdir()
for scale in (1.25,1.5,2.):
 folder=out/str(scale).replace('.','p');folder.mkdir();cfg=dict(biarticular_pair_gain=.3,commissural_gain=.3,learn_reciprocal=True,motor_scale=scale)
 p=folder/'config.json';p.write_text(json.dumps(cfg,indent=2))
 with (folder/'air.log').open('w') as log:subprocess.run([sys.executable,str(H/'check_shared_air.py'),'--weights',str(source),'--output-directory',str(folder/'air_pair'),'--config',str(p)],stdout=log,stderr=subprocess.STDOUT,check=True)
 r=json.loads((folder/'air_pair/SUMMARY.json').read_text());print('AIR',scale,r['accepted'],flush=True)
 if not r['accepted']:continue
 base=json.loads((H/'runs/tonic_abduction_20261030_s1_v1/0p02/standing/task_low.json').read_text());base.update(cfg,weights_directory=str(source))
 for name,config in [('standing',base),('ground',{**base,'stage':4,'force_controlled_support':True,'support_fraction':.7,'regional_contacts':True,'cross_contact_gain':6.,'cross_contact_event_amplitude':2.})]:
  p=folder/(name+'.json');p.write_text(json.dumps(config,indent=2))
  with p.with_suffix('.log').open('w') as log:subprocess.run([sys.executable,str(H/'physical_eval.py'),str(p)],stdout=log,stderr=subprocess.STDOUT,check=True)
  r=json.loads(p.with_suffix('.result.json').read_text());print(name,scale,r['assessment']['accepted'],flush=True)
  if not r['assessment']['accepted']:break
