"""Fixed held-out pressure/map errors and torso pulse; no tuning on these cases."""
import json,copy,argparse
from pathlib import Path
import gait_standing_family as F
from gait_standing_objective import score
import numpy as np

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main(series_domain_recheck=False,consistent_forward_recheck=False):
 best=json.loads((ROOT/'rustuna_syn6_v1/BEST.json').read_text())['trial']['params']
 nominal=ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json'
 config=copy.deepcopy(json.loads(nominal.read_text()))
 for row in config['actuators'].values():row.update(tau_fill_s=.08,tau_vent_s=.12,delay_s=.016)
 config['uncertainty']=dict(force_gain=.8,directional_gain=.12,direction_tau_s=.08,velocity_scale_m_s=.03)
 config['uncertainty_provenance']='Fixed hypothetical held-out stressor, not fitted or measured BPA dynamics: doubled pressure time constants,16msdelay,20%static surface loss,12%loading-direction force modulation with80msmemory. Inverse map remains nominal.'
 heldout=ROOT/'virtual_bpa_design/HELDOUT_SLOW_WEAK_MEMORY.json'
 if heldout.exists() and json.loads(heldout.read_text())!=config:raise ValueError('Preserve existing held-out specification')
 heldout.write_text(json.dumps(config,indent=2),encoding='utf-8')
 replay=ROOT/'pressure_standing/syn6_nominal_fresh_v1.inputs.npz'
 cases=[('slowweak_live',heldout,None,0.),('slowweak_replay',heldout,replay,0.),('torso_push_live',nominal,None,.02),('torso_push_replay',nominal,replay,.02)]
 for label,override in (('maperror',dict(force_gain=.8)),('delay',None)):
  single=copy.deepcopy(json.loads(nominal.read_text()))
  if override:single['uncertainty']=override
  else:
   for row in single['actuators'].values():row['delay_s']=.016
  single['uncertainty_provenance']='Single held-out causal isolation after combined stressor left the declared contraction domain; no neural retuning.'
  spec=ROOT/('virtual_bpa_design/HELDOUT_'+label.upper()+'.json')
  if spec.exists() and json.loads(spec.read_text())!=single:raise ValueError('Preserve existing isolated uncertainty specification')
  spec.write_text(json.dumps(single,indent=2),encoding='utf-8')
  cases.extend([(label+'_live',spec,None,0.),(label+'_replay',spec,replay,0.)])
 if series_domain_recheck:cases=[(tag+'_series_v2',spec,sense,push) for tag,spec,sense,push in cases if 'torso_push' not in tag]
 if consistent_forward_recheck:cases=[(tag+'_equilibrium_v3',spec,sense,push) for tag,spec,sense,push in cases if tag.startswith(('slowweak','maperror'))]
 out=ROOT/('pressure_forward_recheck' if consistent_forward_recheck else 'pressure_series_recheck' if series_domain_recheck else 'pressure_heldouts');out.mkdir(exist_ok=True);results=[]
 for tag,spec,sense,push in cases:
  path=out/('syn6_'+tag+'.result.json')
  if path.exists():r=json.loads(path.read_text());reused=True
  else:r=F.run('syn6',out,tag=tag,gain=best['pf_mn_gain'],ib_group_gain=best['ib_group_gain'],standing_bias_scale=best['standing_bias_scale'],bpa_configuration=spec,sensory_replay=sense,mechanical_push_bw=push);reused=False
  value,a=score(r,path.with_suffix('').with_suffix(''))
  neural_error=None
  trace=path.with_suffix('').with_suffix('.apparatus.npz')
  if trace.exists():
   z=np.load(trace);ref=np.load(ROOT/'pressure_standing/syn6_nominal_fresh_v1.apparatus.npz')
   if z['neuro'].shape==ref['neuro'].shape:neural_error=float(np.max(abs(z['neuro']-ref['neuro'])))
  results.append(dict(case=tag,result=str(path),reused=reused,assessment=a,neural_change_from_nominal_max_mV=neural_error,upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),mesh_height_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),failure=r['failure']))
  summary=dict(machine='easteregg2',completed=len(results),requested_cases=len(cases),controller_parameters=best,results=results,scope='Virtual Gait pressure-actuated standing, original apparatus retained. No retuning on held-outs. Nominal afferent replay disables SNS state corrections while preserving command biases and local length-aware pressure inversion. Declared uncertainty is not hardware identification.')
  (out/'HELDOUT_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
  print(json.dumps(dict(case=tag,accepted=a['accepted'],neural_change_mV=neural_error)),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();group=p.add_mutually_exclusive_group();group.add_argument('--series-domain-recheck',action='store_true');group.add_argument('--consistent-forward-recheck',action='store_true');a=p.parse_args();main(a.series_domain_recheck,a.consistent_forward_recheck)
