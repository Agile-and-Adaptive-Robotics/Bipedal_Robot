"""Two active causal group gains, frozen global standing configuration."""
import json
from pathlib import Path
import gait_standing_family as F
from gait_standing_objective import score
KEYS={'knee_ext_iiin','hip_ext_iiin'}

def handle_request(request):
 cfg=request['manifest'];root=Path(cfg['campaign_root']);out=root/'group_diagnosis';audit=json.loads((out/'PER_MOTOR_CURRENT_AUDIT.json').read_text())
 if request['op']=='preflight':return dict(ok=True,ready_for_optimization=audit['shared_static_targets_bit_exact'] and max(r['membrane_balance_max_error_mV'] for r in audit['families'])<1e-9,scope='Only group-II antagonist conductance onto primary knee_ext/hip_ext, bilaterally; other global winner gains frozen.')
 if request['op']!='evaluate' or set(request['params'])!=KEYS:raise ValueError('Unknown/unconsumed operation or keys')
 p=request['params'];reuse={(1.,1.):'w2lvar',(0.,1.):'knee_ext',(1.,0.):'hip_ext'};case=reuse.get((p['knee_ext_iiin'],p['hip_ext_iiin']))
 if case:
  prefix=out/('w2lvar_'+case+'_v1');r=json.loads(prefix.with_suffix('.result.json').read_text())
 else:
  global_params=json.loads((root/'rustuna_w2lvar_v1/BEST.json').read_text())['trial']['params'];tag=request['evaluation_id'];r=F.run('w2lvar',out/'group_trials',tag=tag,gain=global_params['pf_mn_gain'],ib_group_gain=global_params['ib_group_gain'],standing_bias_scale=global_params['standing_bias_scale'],pathway_group_scales={'IIIN:knee_ext':p['knee_ext_iiin'],'IIIN:hip_ext':p['hip_ext_iiin']},capture_neural_inputs=True);prefix=out/'group_trials'/('w2lvar_'+tag)
 value,a=score(r,prefix);return dict(ok=True,value=value,diagnostics=dict(assessment=a,result=str(prefix.with_suffix('.result.json')),reused_completed_trial=bool(case),upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw')))
