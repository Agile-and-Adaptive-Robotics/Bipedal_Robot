"""Frozen successful standing family transfer; w2lvar remains gated out."""
import sys,json,subprocess
from pathlib import Path
import numpy as np
from gait_standing_objective import score
HERE=Path(__file__).resolve().parent;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main():
 accepted=json.loads((ROOT/'OTHER_FAMILY_STANDING_CHECKPOINT.json').read_text())['families'][0]
 if accepted['family']!='s3k' or not accepted['assessment']['accepted'] or not accepted['bit_exact']:raise ValueError('Fresh standing gate not met')
 out=ROOT/'other_family_pressure';out.mkdir(exist_ok=True);cases=[]
 for tag in ('bio_prepared_v1','bpa_prepared_v1','bpa_prepared_fresh_v1'):
  path=out/('s3k_'+tag+'.result.json')
  if not path.exists():
   cmd=[sys.executable,str(HERE/'gait_standing_family.py'),'--family','s3k','--out',str(out),'--tag',tag,'--parameters-json',str(ROOT/'rustuna_s3k_v1/BEST.json'),'--pinned-actuator-settling']
   if tag.startswith('bpa'):cmd+=['--bpa-configuration',str(ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json'),'--settling-reference-configuration',str(ROOT/'pinned_preparation/PREPARATION_REFERENCE.json')]
   with (out/(tag+'.log')).open('w',encoding='utf-8') as log:done=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
   if done.returncode:raise RuntimeError(tag+' failed; inspect log')
  r=json.loads(path.read_text());a=score(r,path.with_suffix('').with_suffix(''))[1];cases.append(dict(tag=tag,path=str(path),assessment=a,failure=r['failure'],upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m')));print(json.dumps(cases[-1]),flush=True)
  if not a['accepted']:break
 exact=None
 if len(cases)==3:
  exact=True
  for suffix in ('.npz','.apparatus.npz','.pressure.npz','.release_state.npz'):
   with np.load(out/('s3k_bpa_prepared_v1'+suffix)) as a,np.load(out/('s3k_bpa_prepared_fresh_v1'+suffix)) as b:exact=exact and set(a.files)==set(b.files) and all(np.array_equal(a[k],b[k]) for k in a.files)
 report=dict(machine='easteregg2',cases=cases,fresh_bpa_bit_exact=exact,scope='Same frozen s3k standing gains, matched0.6s preparation policy and same baseline rig for biological versus pressure actuation. Virtual BPA sleeve design was previously sized from syn6 standing only; failures are retained and no resizing is performed. w2lvar pressure transfer is not run because its biological standing load/sway gate fails. No gait/curriculum or mounted CAD claim.')
 (out/'S3K_PREPARED_TRANSFER.json').write_text(json.dumps(report,indent=2),encoding='utf-8')

if __name__=='__main__':main()
