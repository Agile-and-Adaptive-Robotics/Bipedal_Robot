"""Check prepared-state equality and preserve explicit apparatus qualifications."""
import json
from pathlib import Path
import numpy as np
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def arrays(a,b):
 with np.load(a) as x,np.load(b) as y:
  if set(x.files)!=set(y.files):raise ValueError('State keys differ')
  return {k:dict(shape_equal=x[k].shape==y[k].shape,bit_exact=x[k].shape==y[k].shape and np.array_equal(x[k],y[k]),max_abs_difference=float(np.max(abs(x[k]-y[k]))) if x[k].shape==y[k].shape and x[k].size else 0.) for k in x.files if x[k].dtype.kind in 'fiub'}

def main():
 folder=ROOT/'pinned_preparation';ledger=json.loads((folder/'SUMMARY.json').read_text());base=folder/'syn6_prepared_nominal_v1';results={}
 if 'nominal_fresh'in ledger['cases']:
  other=folder/'syn6_prepared_nominal_fresh_v1';checks={}
  for suffix in ('.npz','.apparatus.npz','.pressure.npz','.release_state.npz'):
   checks.update({suffix+':'+k:v for k,v in arrays(base.with_suffix(suffix),other.with_suffix(suffix)).items()})
  results['nominal_fresh']=dict(accepted=ledger['cases']['nominal_fresh']['assessment']['accepted'],bit_exact=all(v['bit_exact']for v in checks.values()),comparisons=checks)
 for name,row in ledger['cases'].items():
  if name.startswith('bio_nominal'):
   with np.load(row['release_state']) as biological:results[name]=dict(actuator_class='native biological muscle',activation_mean=float(biological['act'].mean()),activation_max=float(biological['act'].max()),scope='Different actuator-state class; not compared to pneumatic stored-state arrays')
   continue
  path=Path(row['release_state']);checks=arrays(base.with_suffix('.release_state.npz'),path)
  results.setdefault(name,{})['release_matches_common_nominal']=all(v['bit_exact']for v in checks.values())
  results[name]['release_comparisons']=checks
  prepared=path.with_name(path.name.replace('.release_state.npz','.prepared_state.npz'))
  if prepared.exists():results[name]['before_reset_matches_common_nominal']=all(v['bit_exact'] for v in arrays(base.with_suffix('.release_state.npz'),prepared).values())
 if 'nominal_replay'in ledger['cases']:
  other=folder/'syn6_prepared_nominal_replay_v1';checks={}
  for suffix in ('.npz','.apparatus.npz','.pressure.npz'):
   checks.update({suffix+':'+k:v for k,v in arrays(base.with_suffix(suffix),other.with_suffix(suffix)).items()})
  results['nominal_replay']['whole_trace_bit_exact']=all(v['bit_exact']for v in checks.values())
 with np.load(base.with_suffix('.release_state.npz')) as initial:
  release=dict(pressure_mean_kpa=float(initial['pressure'].mean()),pressure_max_kpa=float(initial['pressure'].max()),series_spring_energy_J=float(initial['series_spring_energy_J'].sum()),activation_max=float(initial['act'].max()),physics_time_s=float(initial['physics_time']))
 pairs={}
 for left,right in (('harsh_live','harsh_replay'),('harsh_late_live','harsh_late_replay'),('harsh_zero_live','harsh_zero_replay'),('harsh_cold_live','harsh_cold_replay')):
  if left in ledger['cases'] and right in ledger['cases']:pairs[left+'__'+right]=all(v['bit_exact'] for v in arrays(Path(ledger['cases'][left]['release_state']),Path(ledger['cases'][right]['release_state'])).values())
 selected={}
 if 'selected_fresh'in ledger['cases']:
  chosen=json.loads((folder/'SELECTED_RUNG.json').read_text());a=Path(chosen['result']).with_suffix('').with_suffix('');b=folder/'syn6_prepared_selected_fresh_v1';c={}
  for suffix in ('.npz','.apparatus.npz','.pressure.npz','.release_state.npz'):c.update({suffix+':'+k:v for k,v in arrays(a.with_suffix(suffix),b.with_suffix(suffix)).items()})
  selected=dict(accepted=ledger['cases']['selected_fresh']['assessment']['accepted'],bit_exact=all(v['bit_exact']for v in c.values()),comparisons=c)
 biological_fresh={}
 if 'bio_nominal_fresh'in ledger['cases']:
  a=folder/'syn6_prepared_bio_nominal_v1';b=folder/'syn6_prepared_bio_nominal_fresh_v1';checks={}
  for suffix in ('.npz','.apparatus.npz','.release_state.npz'):checks.update({suffix+':'+k:v for k,v in arrays(a.with_suffix(suffix),b.with_suffix(suffix)).items()})
  biological_fresh=dict(accepted=ledger['cases']['bio_nominal_fresh']['assessment']['accepted'],bit_exact=all(v['bit_exact']for v in checks.values()),comparisons=checks)
 report=dict(machine='easteregg2',release_summary=release,checks=results,matched_live_replay_release_pairs=pairs,selected_rung_fresh=selected,biological_fresh=biological_fresh,nominal_exact_fresh_gate=bool(results.get('nominal_fresh',{}).get('accepted') and results.get('nominal_fresh',{}).get('bit_exact')),scope='Arrays include neural persistent-channel gates, qpos/qvel/activation, pressure and full command history, modeled series-spring energy, and engine warm-start state. Pneumatic energy is not quantified without a volume model.')
 (folder/'STATE_MATCH_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in report.items() if k!='checks'},indent=2))

if __name__=='__main__':main()
