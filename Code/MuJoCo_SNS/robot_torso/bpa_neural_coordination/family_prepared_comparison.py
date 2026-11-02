"""Durable scoped comparison; legacy act is never treated as activation."""
import json,csv
from pathlib import Path
import numpy as np
from gait_standing_objective import score
from native_activation_replay import replay
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def row(family,encoding,prefix,recipe):
 r=json.loads(prefix.with_suffix('.result.json').read_text());a=score(r,prefix)[1];data=dict(family=family,encoding=encoding,force_model='virtual BPA pressure' if r.get('bpa_configuration') else 'biological native muscle',recipe=recipe,result=str(prefix.with_suffix('.result.json')),accepted=a.get('accepted',False),failure=r['failure'],study_duration_s=r['duration_requested_s'],completed_physics_samples=r['completed_samples'],free_physics_s=r['completed_samples']*.002,tilt_deg=a.get('tilt_deg'),sway_mm=a.get('sway_m',0)*1000 if a.get('sway_m') is not None else None,ground_load_bw=r.get('mean_ground_load_bw'),upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),rig_scales=r.get('rig_scales',dict(rotational=1.,lateral=1.,vertical=1.)),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),knee_max_deg=r.get('max_knee_angle_deg'),actual_release_activation_mean=None,actual_postrelease_activation_mean=None,activation_provenance='Not inferred from legacy act field',postrelease_pressure_mean_kpa=None,postrelease_pressure_peak_kpa=None,release_pressure_mean_kpa=None,release_pressure_peak_kpa=None,force_units='N',frame=r.get('frame'))
 if prefix.with_suffix('.release_state.npz').exists():
  with np.load(prefix.with_suffix('.release_state.npz')) as z:
   if not r.get('bpa_configuration'):data['actual_release_activation_mean']=float(z['act'].mean());data['activation_provenance']='Actual native activation state at release'
   if 'pressure'in z:data['release_pressure_mean_kpa']=float(z['pressure'].mean());data['release_pressure_peak_kpa']=float(z['pressure'].max())
 if prefix.with_suffix('.motor_state.npz').exists():
  with np.load(prefix.with_suffix('.motor_state.npz')) as z:
   data['mean_neural_excitation']=float(z['excitation'].mean())
   if not r.get('bpa_configuration'):data['actual_postrelease_activation_mean']=float(z['activation'].mean());data['activation_provenance']='Actual full native activation telemetry; release state separately listed'
   else:data['activation_provenance']='No biological activation drives BPA force; pressure/tension are causal. Native act bookkeeping is deliberately not presented as BPA activation.'
 elif prefix.with_suffix('.npz').exists():
  with np.load(prefix.with_suffix('.npz')) as z:data['mean_neural_excitation']=float(z['ctrl'][-r['completed_samples']:].mean()) if r['completed_samples'] else None
 else:data['mean_neural_excitation']=None;data['activation_provenance']+='; legacy run stopped before native log finalization; no synthetic replacement'
 if prefix.with_suffix('.pressure.npz').exists():
  with np.load(prefix.with_suffix('.pressure.npz')) as z:
   pressure=z['samples'][:,1,:];data['postrelease_pressure_mean_kpa']=float(pressure.mean());data['postrelease_pressure_peak_kpa']=float(pressure.max());data['pressure_active_routes']=len(z['names'])
 if not r.get('bpa_configuration') and prefix.with_suffix('.release_state.npz').exists() and prefix.with_suffix('.npz').exists():
  values,provenance=replay(prefix);data['native_activation_replay_validation']=provenance
  if data['actual_postrelease_activation_mean'] is None:data['actual_postrelease_activation_mean']=float(values.mean());data['activation_provenance']='Exact native-kernel state reconstruction from actual measured release activation and full excitation; method verified against all92 new prepared native activation channels'
 return data

def main():
 specs=[('s3k','graded',ROOT/'other_family_pressure/s3k_bio_prepared_v1','shared standing gains'),('s3k','graded',ROOT/'other_family_pressure/s3k_bpa_prepared_v1','frozen transfer; failed domain preserved'),('syn6','graded',ROOT/'pinned_preparation/syn6_prepared_bio_nominal_v1','shared standing gains'),('syn6','graded',ROOT/'pinned_preparation/syn6_prepared_nominal_v1','frozen transfer'),('syn6','sensory LIF hybrid',ROOT/'prepared_spiking/syn6_bio_spike_v2','frozen sensory substitution'),('syn6','sensory LIF hybrid',ROOT/'prepared_spiking/syn6_bpa_spike_v2','frozen sensory substitution')]
 for case in ('bio_graded','bio_spike','bpa_graded','bpa_spike'):
  specs.append(('w2lvar mean-IIIN','sensory LIF hybrid' if 'spike'in case else 'graded',ROOT/('normalized_prepared_transfer/w2lvar_'+case+'_v1'),'all IIIN connections retained; standing-specific normalized mapping'))
 rows=[row(*spec) for spec in specs if spec[2].with_suffix('.result.json').exists()]
 report=dict(machine='easteregg2',rows=rows,scope=['Descriptive family comparison, not an equal-training-compute ranking. Global standing fits used12trials each; w2lvar additionally required causal mapping diagnosis and an explicit mean-pool variant.','Same corrected-pose body, force calibration,0.6s pinned preparation policy,8s study, and baseline apparatus for the listed nominal cases. Each actuator/encoder treatment naturally has different internal release states; equality is required within nominal/replay/fresh pairs, not between force models.','Virtual BPA dimensions were sized from syn6 standing and remain unchanged for all families, so this tests transfer into that declared design, not intrinsic controller superiority. Static Festo surface is characterized; dimensions, pressure lags, delay, and directional memory are engineering assumptions. No mounted CAD identification.','Biological activation is from actual release/native telemetry, or exactly replayed through the native activation kernel using measured release activation and full excitation. The replay is independently verified against all92 newly logged prepared activation channels. Each row labels provenance. Legacy act is excitation. BPA act bookkeeping does not drive force and is not interpreted as muscle activation.','Pressure values are actual model states. Postrelease peaks do not claim the maximum over unlogged preparation; final release pressure is listed separately. Domain failures may have only partial postrelease records.','Six pruned routes are disabled in the86-route pressure design; biological model retains92native routes. Strong lateral/rotational guidance and logged mesh penetration remain. No walking, independent standing, or complete five-stage curriculum claim.'])
 out=ROOT/'family_comparison';out.mkdir(exist_ok=True);(out/'CURRENT_PREPARED_FAMILY_COMPARISON.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 fields=['family','encoding','force_model','accepted','study_duration_s','completed_physics_samples','free_physics_s','tilt_deg','sway_mm','ground_load_bw','upward_rig_bw','downward_rig_bw','actual_release_activation_mean','actual_postrelease_activation_mean','mean_neural_excitation','release_pressure_mean_kpa','release_pressure_peak_kpa','postrelease_pressure_mean_kpa','postrelease_pressure_peak_kpa','activation_provenance','failure','result']
 with (out/'CURRENT_PREPARED_FAMILY_COMPARISON.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
 common={};reference=specs[0][2]
 with np.load(reference.with_suffix('.release_state.npz')) as z:reference_pose={k:z[k].copy() for k in ('qpos','qvel')}
 ref=json.loads(reference.with_suffix('.result.json').read_text())
 for spec in specs:
  prefix=spec[2]
  if not prefix.with_suffix('.release_state.npz').exists():continue
  r=json.loads(prefix.with_suffix('.result.json').read_text())
  with np.load(prefix.with_suffix('.release_state.npz')) as z:checks={k:bool(np.array_equal(reference_pose[k],z[k])) for k in reference_pose}
  checks['static_targets_exact']=bool(np.array_equal(ref['standing_activation_targets'],r['standing_activation_targets']));checks['force_law_same']=ref['force_law']==r['force_law'];common[str(prefix)]=checks
 report['common_mechanical_reference_audit']=common
 followup=[]
 for case in ('rot0p8_graded','rot0p8_spike','hold20_graded','hold20_spike'):
  path=ROOT/('normalized_prepared_transfer/w2lvar_'+case+'_v1.comparison.json')
  if path.exists():
   c=json.loads(path.read_text());r=json.loads(Path(c['result']).read_text());followup.append(dict(case=case,accepted=c['assessment']['accepted'],free_physics_s=r['completed_samples']*.002,failure=r['failure'],result=c['result']))
 report['normalized_w2lvar_bpa_followups']=followup
 valve_summary=ROOT/'pressure_reserve_valve/SUMMARY.json'
 valve=[]
 if valve_summary.exists():
  payload=json.loads(valve_summary.read_text())
  tags={'w2lvar_winner_fresh','w2lvar_zero_current','w2lvar_heldout_delay','w2lvar_off_delay','w2lvar_heldout_map','w2lvar_off_map','w2lvar_heldout_pulse','w2lvar_off_pulse','w2lvar_heldout_harsh','w2lvar_off_harsh'}
  valve=[{k:r[k] for k in ('case','condition','accepted','moving_s','failure','result')} for r in payload['rows'] if r['case'] in tags]
  report['distinct_neural_valve_followups']=dict(scope='Additional two-neuron engineering architecture, not unchanged original controller.20s study includes0.6s preparation. Nominal/pulse pass; delay/map/combined generalization rejected. No post-heldout reselection.',rows=valve,record=str(valve_summary.with_name('EXPERIMENT_RECORD.md')))
 (out/'CURRENT_PREPARED_FAMILY_COMPARISON.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 text='# Prepared family comparison — easteregg2\n\n| Family | Encoding | Force | Pass | Up/down rig (%BW) | Mean moving activation | Actual postrelease P mean/peak (kPa) |\n| --- | --- | --- | --- | ---: | ---: | ---: |\n'
 for v in rows:
  act='—' if v['actual_postrelease_activation_mean'] is None else f"{v['actual_postrelease_activation_mean']:.4f}";pressure='—' if v['postrelease_pressure_mean_kpa'] is None else f"{v['postrelease_pressure_mean_kpa']:.1f}/{v['postrelease_pressure_peak_kpa']:.1f}";up='—' if v['upward_rig_bw'] is None else f"{100*v['upward_rig_bw']:.3f}/{100*v['downward_rig_bw']:.3f}";text+=f"| {v['family']} | {v['encoding']} | {v['force_model']} | {v['accepted']} | {up} | {act} | {pressure} |\n"
 text+='\nThe nominal table covers 8 s studies (7.4 s moving after pinned preparation). Longer and looser-guidance tests are separate:\n\n'
 for v in followup:text+='- '+v['case']+': '+('pass' if v['accepted'] else 'FAIL')+f", {v['free_physics_s']:.3f} s moving; "+str(v['failure'])+'\n'
 if valve:
  text+='\nDistinct two-neuron valve-reference overlay (20 s study, 19.4 s moving when complete; not the unchanged nominal controller):\n\n| Case | Pass | Moving (s) |\n| --- | --- | ---: |\n'
  for v in valve:text+=f"| {v['case']} | {v['accepted']} | {v['moving_s']:.3f} |\n"
  text+='\nNominal/fresh and pulse pass. Longer-delay, force-map-error, and combined tests fail, so general robustness is rejected. The frozen selection was not replaced after held-out results. Full circuit, actual pressure checks, signed apparatus forces, mesh penetration, and failure details are in `'+str(valve_summary.with_name('EXPERIMENT_RECORD.md'))+'`. This overlay has not been matched across s3k/syn6 and is not a family ranking.\n'
 text+='\n'+''.join('- '+s+'\n' for s in report['scope']);text+='\nFull physical durations, signed apparatus forces, actual activation provenance, force/pressure values, failure details, and absolute result paths are in the JSON/CSV. All cases retain their individual raw records.\n'
 (out/'CURRENT_PREPARED_FAMILY_COMPARISON.md').write_text(text,encoding='utf-8');print(json.dumps(dict(rows=len(rows),folder=str(out)),indent=2))

if __name__=='__main__':main()
