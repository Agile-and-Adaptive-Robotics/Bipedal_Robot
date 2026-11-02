"""Small shared-bilateral Gait standing calibration, not a walking objective."""
import os,sys,json,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;SPINAL=HERE.parents[1]/'spinal'
sys.path[:0]=[str(HERE),str(SPINAL/'functional_gait2392')]
import gait_standing_family as F
from curriculum_metrics import assess
from frontal_metrics import assess_frontal
KEYS={'pf_mn_gain','ib_group_gain','standing_bias_scale'}

def score(result,path):
 if result['failure'] is not None or result['completed_samples']<3600:return 10000.,dict(accepted=False,safety_failure=result['failure'])
 a=assess(str(path.with_suffix('.npz')),result['metrics'],3);frontal=assess_frontal(str(path.with_suffix('.npz')));a.update(frontal)
 knees=result.get('max_knee_angle_deg') or [99.,99.];a['knee_extension_gate']=max(knees)<=0.;a['measured_ground_load_gate']=result['mean_ground_load_bw']>.5
 a['accepted']=bool(a['accepted'] and a['frontal_gate'] and a['knee_extension_gate'] and a['measured_ground_load_gate'])
 # Existing stage3 score with stronger distance-to-sway criterion; failure never a pass.
 penalty=max(0.,a['sway_m']-.06)*1000+max(0.,a['tilt_deg']-15.)*3
 return float(-a['score']+penalty+(0. if a['accepted'] else 500.)),a

def handle_request(request):
 cfg=request['manifest'];family=cfg['family'];root=Path(cfg['campaign_root']);baseroot=root/'standing'
 if family not in ('s3k','w2lvar','syn6'):raise ValueError('Unknown Gait neural family')
 if request['op']=='preflight':
  gate=json.loads((root/'family_parity'/f'{family}_parity.json').read_text());baseline=json.loads((baseroot/f'{family}_baseline_v2.result.json').read_text())
  ready=bool(gate['sparse_pairwise_bit_exact'] and baseline['failure'] is None and max(baseline['max_knee_angle_deg'])<=0)
  return dict(ok=True,ready_for_optimization=ready,family=family,actuation='biological Hill',baseline_evaluated_not_passed=True,baseline_result=str(baseroot/f'{family}_baseline_v2.result.json'),scope='Shared standing calibration only; no stepping/curriculum promotion or BPA result')
 if request['op']!='evaluate':raise ValueError('Unknown operation')
 params=request['params']
 if set(params)!=KEYS:raise ValueError('Unused/missing optimization parameters')
 cached=None
 if params==dict(pf_mn_gain=1.,ib_group_gain=1.,standing_bias_scale=1.):cached=baseroot/f'{family}_baseline_v2'
 if params==dict(pf_mn_gain=0.,ib_group_gain=1.,standing_bias_scale=1.):cached=baseroot/f'{family}_pf_motor0_v1'
 if cached is not None:
  r=json.loads(cached.with_suffix('.result.json').read_text());path=cached;reused=True
 else:
  out=root/'tuning_trials';tag=request['evaluation_id'];r=F.run(family,out,duration=8.,tag=tag,gain=params['pf_mn_gain'],standing_bias_scale=params['standing_bias_scale'],ib_group_gain=params['ib_group_gain']);path=out/(family+'_'+tag);reused=False
 value,a=score(r,path);return dict(ok=True,value=value,diagnostics=dict(assessment=a,result=str(path.with_suffix('.result.json')),reused_completed_trial=reused,mean_upward_rig_support_bw=r.get('mean_upward_rig_support_bw'),knee_max_deg=r.get('max_knee_angle_deg')))
