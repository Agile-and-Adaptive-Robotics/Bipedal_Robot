"""New registered internal-capacity / measured-tension neural comparison."""
import json,copy,hashlib,shutil
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');DEST=ROOT/'pressure_reserve_capacity'

def main():
 DEST.mkdir(exist_ok=True);folder=DEST/'fixtures';folder.mkdir(exist_ok=True);old=json.loads((ROOT/'pressure_reserve_uncertainty/PREREGISTRATION.json').read_text());base=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text())
 train=[]
 for f in old['fixtures']:
  if f['split']=='train':train.append(dict(f,sensor_noise_seed=811000+len(train)))
 validation=[dict(name='new_valid_interior',seed=911031,fill=[.058,.073],delay=[.008,.014],gain=.85,pulse_bw=0.,sensor_noise_seed=921031),dict(name='new_valid_slow',seed=911032,fill=[.085,.095],delay=[.018,.022],gain=.9,pulse_bw=0.,sensor_noise_seed=921032),dict(name='new_valid_weak',seed=911033,fill=[.047,.067],delay=[.006,.012],gain=.78,pulse_bw=0.,sensor_noise_seed=921033),dict(name='new_valid_pulse',seed=911034,fill=[.04,.04],delay=[.004,.004],gain=1.,pulse_bw=.10,sensor_noise_seed=921034)]
 for spec in validation:
  rng=np.random.default_rng(spec['seed']);cfg=copy.deepcopy(base);pairs={}
  for name,row in cfg['actuators'].items():
   pair=name[:-2] if name.endswith(('_r','_l')) else name
   if pair not in pairs:pairs[pair]=(float(rng.uniform(*spec['fill'])) if spec['fill'][0]!=spec['fill'][1] else spec['fill'][0],.002*round((float(rng.uniform(*spec['delay'])) if spec['delay'][0]!=spec['delay'][1] else spec['delay'][0])/.002))
   fill,delay=pairs[pair];row.update(tau_fill_s=fill,tau_vent_s=1.5*fill,delay_s=delay)
  cfg['uncertainty']={'force_gain':spec['gain']};cfg['fixture_provenance']=dict(seed=spec['seed'],scope='New untouched registered validation; no actual casegain supplied to neural estimators')
  path=folder/(spec['name']+'.json');data=json.dumps(cfg,indent=2)
  if path.exists() and path.read_text()!=data:raise ValueError('Validation fixture changed')
  path.write_text(data,encoding='utf-8');spec.update(split='validation',path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
 common=dict(machine='easteregg2',training=train,validation=validation,prior_data='All earlier nominal and stress/uncertainty-training results are known training evidence. Previous unused validation fixtures remain unused; this study registers a new validation list before fitting.',common_preparation=old['common_preparation'],gates=old['gates'],measurement=dict(additional_transport_delay_s=.006,engine_sensor_latency_s=.002,pressure_noise_std_kpa=.02,tension_noise_std_N=.05,noise_clip_sigma=3.,noise_policy='SeededPCG64 precomputed table; raw pressure/tension/routeL/v delayed together; noise begins after release; no runtime RNG draws',scope='Declared model sensor assumptions, not hardware-identified latency/noise. Existing coreIb inputs remain unchanged; only new reserve branch has these sensor models.'),architecture=dict(new_cells=2,leak_uS=1.,capacitance_nF=6.,relay_tau_s=.006,input_current_cap_nA=10.,derivative_filter_s=.01,reserve_target_N=1.,symmetric_all_parameters=True,actual_Fmax_edited=False,actual_uncertainty_gain_supplied=False),modes=dict(pressure=dict(search_space=dict(capacity_gain=dict(low=.75,high=1.),pressure_cap_kpa=dict(low=1.,high=5.),preview_s=dict(low=.10,high=.18)),fixed_parameters=dict(estimator_mode='pressure'),initial_candidates=[dict(capacity_gain=.8,pressure_cap_kpa=5.,preview_s=.18),dict(capacity_gain=.85,pressure_cap_kpa=4.,preview_s=.16),dict(capacity_gain=.75,pressure_cap_kpa=5.,preview_s=.18)]),force=dict(search_space=dict(pressure_cap_kpa=dict(low=1.,high=5.),preview_s=dict(low=.10,high=.18)),fixed_parameters=dict(estimator_mode='force'),initial_candidates=[dict(pressure_cap_kpa=1.,preview_s=.10),dict(pressure_cap_kpa=3.,preview_s=.16),dict(pressure_cap_kpa=5.,preview_s=.18)])),trials_per_mode=4,wall_budget_per_mode_s=1400.,trial_timeout_s=450.,objective='Same hard standing/knee/apparatus/mesh/domain gates; prioritize alltraining success, then worst sway/tilt/rig burden plus actual addedreference/output saturation effort; no metabolicenergy claim.',validation_gate='For each mode independently: alltraining acceptance,30s nominal stability and fresh complete-state parity before new validation; no postvalidation parameter/mode reselection.',telemetry='Full traces RAM for exact gates; compact windows/scalars/hashes for search. Full default parity/preparation and at most2unique causal failures; full selectedfinalist/fresh traces only.')
 path=DEST/'PREREGISTRATION.json';data=json.dumps(common,indent=2)
 if path.exists() and path.read_text()!=data:raise ValueError('Preregistration changed')
 path.write_text(data,encoding='utf-8')
 for i,mode in enumerate(('pressure','force')):
  plan=common['modes'][mode];manifest=dict(adapter_path=str(HERE/'capacity_comparison_adapter.py'),output_dir=str(DEST/mode/'rustuna'),registration_path=str(path),registration_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),estimator_mode=mode,search_space=plan['search_space'],fixed_parameters=plan['fixed_parameters'],initial_candidates=plan['initial_candidates'],seed=1011031+i,max_trials=4,n_startup_trials=3,wall_budget_s=1400.,trial_timeout_s=450.)
  (DEST/(mode+'_MANIFEST.json')).write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 snapshots=DEST/'source_snapshot';snapshots.mkdir(exist_ok=True);sources={}
 for name in ('register_capacity_comparison.py','capacity_comparison_adapter.py','capacity_estimator_neural.py','compact_trial_traces.py','pressure_valve_neural.py','pressure_reserve_neural.py','reserve_valve_runner.py','inverse_feasibility_policy.py'):
  p=HERE/name
  if p.exists():
   digest=hashlib.sha256(p.read_bytes()).hexdigest();target=snapshots/(p.stem+'_'+digest[:8]+p.suffix)
   if not target.exists():shutil.copyfile(p,target)
   sources[str(p)]=dict(sha256=digest,snapshot=str(target))
 (snapshots/'MANIFEST.json').write_text(json.dumps(sources,indent=2),encoding='utf-8');print(json.dumps(dict(registration=str(path),modes=['pressure','force'],trials_each=4)))

if __name__=='__main__':main()
