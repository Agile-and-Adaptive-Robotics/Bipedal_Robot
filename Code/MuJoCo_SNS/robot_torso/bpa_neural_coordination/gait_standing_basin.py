"""Restartable eight-corner ±1% standing-parameter check, one physics worker."""
import itertools,json
from pathlib import Path
from gait_standing_family import run
from gait_standing_objective import score

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main():
 best=json.loads((ROOT/'rustuna_syn6_v1/BEST.json').read_text());params=best['trial']['params']
 keys=('pf_mn_gain','ib_group_gain','standing_bias_scale');out=ROOT/'standing_basin';out.mkdir(exist_ok=True)
 results=[]
 for signs in itertools.product((-1,1),repeat=3):
  label=''.join('p' if s>0 else 'm' for s in signs);tag='basin_'+label
  candidate={k:params[k]*(1.+.01*s) for k,s in zip(keys,signs)}
  path=out/('syn6_'+tag+'.result.json')
  if path.exists():
   r=json.loads(path.read_text());reused=True
   actual=dict(pf_mn_gain=r['shared_pf_mn_gain'],ib_group_gain=r['ib_group_gain'],standing_bias_scale=r['standing_bias_scale'])
   if actual!=candidate:raise ValueError('Existing corner has a different parameter assignment')
  else:
   r=run('syn6',out,tag=tag,gain=candidate['pf_mn_gain'],ib_group_gain=candidate['ib_group_gain'],standing_bias_scale=candidate['standing_bias_scale']);reused=False
  value,assessment=score(r,path.with_suffix('').with_suffix(''))
  results.append(dict(corner=label,parameters=candidate,reused=reused,result=str(path),assessment=assessment,value=value,upward_rig_bw=r['mean_upward_rig_support_bw'],visual_mesh_height_range_m=r['actual_visual_foot_mesh_clearance_range_m'],max_knee_deg=r['max_knee_angle_deg']))
  summary=dict(machine='easteregg2',family='syn6',condition='Biological standing; all mirrored weights constrained together',relative_perturbation=.01,fit_parameter_count=3,requested_corners=8,completed=len(results),accepted=sum(x['assessment']['accepted'] for x in results),all_corners_complete=len(results)==8,baseline_parameters=params,results=results,scope='Eight-second standing-mode ±1% fitted-parameter corners only; not a walking, harness-free, or force-model uncertainty gate.')
  (out/'BASIN_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
  print(json.dumps(dict(corner=label,accepted=assessment['accepted'],completed=len(results))),flush=True)

if __name__=='__main__':main()
