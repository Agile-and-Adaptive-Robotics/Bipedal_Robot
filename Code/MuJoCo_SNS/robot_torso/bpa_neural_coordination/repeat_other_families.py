"""Fresh-process finalist repeats, including an informative nonpassing candidate."""
import subprocess,sys,json
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
from gait_standing_objective import score

def main():
 rows=[]
 for family in ('s3k','w2lvar'):
  folder=ROOT/('rustuna_'+family+'_v1');best=json.loads((folder/'BEST.json').read_text());original=Path(best['trial']['diagnostics']['result']);out=folder/'fresh';target=out/(family+'_exact_v1.result.json')
  if not target.exists():
   out.mkdir(exist_ok=True)
   with (out/'run.log').open('w',encoding='utf-8') as log:
    done=subprocess.run([sys.executable,str(HERE/'gait_standing_family.py'),'--family',family,'--out',str(out),'--tag','exact_v1','--parameters-json',str(folder/'BEST.json')],stdout=log,stderr=subprocess.STDOUT)
   if done.returncode:raise RuntimeError(family+' repeat failed')
  r=json.loads(target.read_text());prefix=target.with_suffix('').with_suffix('');checks={}
  for suffix in ('.npz','.apparatus.npz','.apparatus_extended.npz'):
   a=original.with_suffix('').with_suffix(suffix);b=prefix.with_suffix(suffix)
   if not a.exists():continue
   with np.load(a) as za,np.load(b) as zb:
    if set(za.files)!=set(zb.files):raise ValueError('Field mismatch')
    checks[suffix]={k:bool(np.array_equal(za[k],zb[k])) for k in za.files}
  exact=all(v for group in checks.values() for v in group.values());row=dict(family=family,params=best['trial']['params'],fresh_result=str(target),assessment=score(r,prefix)[1],bit_exact=exact,comparisons=checks,upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),knee_max_deg=r.get('max_knee_angle_deg'));rows.append(row)
 report=dict(machine='easteregg2',families=rows,scope='Same twelve-trial/600s biological standing calibration as syn6. Frozen exact-JSON fresh repeats. A best numerical objective is not a pass; failed ground-load/stance gates remain failed. PF/aggregate-Ib gains are stage-specific overrides, not a walking controller or faithful unchanged literature replication.')
 (ROOT/'OTHER_FAMILY_STANDING_CHECKPOINT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
