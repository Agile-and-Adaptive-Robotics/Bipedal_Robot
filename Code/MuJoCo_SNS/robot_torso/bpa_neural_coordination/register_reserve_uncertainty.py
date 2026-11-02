"""Preregister fixtures/search/independent validation before uncertainty fitting."""
import json,copy,hashlib,shutil
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
DEST=ROOT/'pressure_reserve_uncertainty'

def main():
 DEST.mkdir(exist_ok=True);fixtures=DEST/'fixtures';fixtures.mkdir(exist_ok=True)
 base=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text())
 prep=json.loads((ROOT/'pinned_preparation/PREPARATION_REFERENCE.json').read_text());prep['pressure_history_capacity_s']=.024
 prep_path=fixtures/'COMMON_PREPARATION_24MS_HISTORY.json';prep_text=json.dumps(prep,indent=2)
 if prep_path.exists() and prep_path.read_text()!=prep_text:raise ValueError('Common preparation changed')
 prep_path.write_text(prep_text,encoding='utf-8')
 specs=[dict(name='train_nominal',split='train',seed=601031,fill=[.04,.04],delay=[.004,.004],gain=1.,pulse_bw=0.),dict(name='train_mid',split='train',seed=601032,fill=[.055,.065],delay=[.010,.014],gain=.9,pulse_bw=0.),dict(name='train_strong',split='train',seed=601033,fill=[.072,.080],delay=[.014,.016],gain=.8,pulse_bw=0.),dict(name='valid_interior',split='validation',seed=701031,fill=[.062,.070],delay=[.008,.010],gain=.85,pulse_bw=0.),dict(name='valid_delay_extrapolation',split='validation',seed=701032,fill=[.082,.090],delay=[.018,.020],gain=.9,pulse_bw=0.),dict(name='valid_map_extrapolation',split='validation',seed=701033,fill=[.045,.055],delay=[.006,.008],gain=.75,pulse_bw=0.),dict(name='valid_new_pulse',split='validation',seed=701034,fill=[.050,.060],delay=[.010,.010],gain=.95,pulse_bw=.075)]
 for spec in specs:
  cfg=copy.deepcopy(base);rng=np.random.default_rng(spec['seed']);pairs={}
  for name,row in cfg['actuators'].items():
   pair=name[:-2] if name.endswith(('_r','_l')) else name
   if pair not in pairs:pairs[pair]=(float(rng.uniform(*spec['fill'])) if spec['fill'][0]!=spec['fill'][1] else spec['fill'][0],.002*round((float(rng.uniform(*spec['delay'])) if spec['delay'][0]!=spec['delay'][1] else spec['delay'][0])/.002))
   fill,delay=pairs[pair];row.update(tau_fill_s=fill,tau_vent_s=1.5*fill,delay_s=delay)
  cfg['uncertainty']={'force_gain':spec['gain']};cfg['fixture_provenance']=dict(seed=spec['seed'],scope='New preregistered paired-route pressure-response fixture; geometry unchanged. Actual force gain never supplied to neural reserve estimator.')
  p=fixtures/(spec['name']+'.json');data=json.dumps(cfg,indent=2)
  if p.exists() and p.read_text()!=data:raise ValueError('Registered fixture changed')
  p.write_text(data,encoding='utf-8');spec['path']=str(p);spec['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
 registration=dict(machine='easteregg2',architecture='Existing two-cell reserve-to-valve reference architecture; fixed1N equivalent reserve target,6ms neural leak,10ms derivative filter,10nA input cap',earlier_data='All prior nominal/delay/map/pulse/harsh observations are exposed prior information, not new validation. New seeded profiles and combinations below are fixed before fitting.',fixtures=specs,search_space=dict(pressure_cap_kpa=dict(low=1.,high=5.),preview_s=dict(low=.08,high=.18),sensor_gain=dict(low=.5,high=2.)),fixed_parameters={},max_trials=8,wall_budget_s=2400.,trial_timeout_s=450.,seed=810031,n_startup_trials=4,initial_candidates=[dict(pressure_cap_kpa=1.,preview_s=.1,sensor_gain=1.),dict(pressure_cap_kpa=3.,preview_s=.15,sensor_gain=1.),dict(pressure_cap_kpa=5.,preview_s=.18,sensor_gain=1.5),dict(pressure_cap_kpa=4.,preview_s=.12,sensor_gain=2.)],gates=dict(original_standing_accepted=True,knee_extension_max_deg=0.,minimum_ground_load_bw=.5,maximum_mean_upward_rig_bw=.05,maximum_mean_downward_rig_bw=.05,minimum_visual_mesh_clearance_m=-.008),objective='First minimize failed training conditions and physical incompleteness; all-pass cases minimize worst sway/tilt/apparatus burden plus actual added-reference and output-saturation effort. Pressure integral is not metabolic or pneumatic energy.',validation_gate='All training cases accepted, nominal stability, and fresh complete-state parity before any new validation. No validation-driven reselection.',telemetry='Full arrays in RAM for gates. Sparse/windowed exploratory telemetry only; full traces for baselines/unique causal failures/finalists/fresh checks. One experiment source snapshot, none pertrial.')
 registration['common_preparation']=dict(path=str(prep_path),sha256=hashlib.sha256(prep_path.read_bytes()).hexdigest(),history_capacity_s=.024,scope='Same nominal pressure preparation for every case. Longer ring capacity accommodates new20ms validation delays; nominal dynamics/parity must match prior winner. No lag/pressure-state change.')
 p=DEST/'PREREGISTRATION.json';data=json.dumps(registration,indent=2)
 if p.exists() and p.read_text()!=data:raise ValueError('Preregistration changed')
 p.write_text(data,encoding='utf-8')
 manifest=dict(registration);manifest.update(adapter_path=str(HERE/'reserve_uncertainty_adapter.py'),output_dir=str(DEST/'rustuna'),registration_path=str(p),registration_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 (DEST/'RUSTUNA_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
 snapshot=DEST/'source_snapshot';snapshot.mkdir(exist_ok=True);sources={}
 for name in ('register_reserve_uncertainty.py','reserve_uncertainty_adapter.py','compact_trial_traces.py','pressure_valve_neural.py','pressure_reserve_neural.py','reserve_valve_runner.py','inverse_feasibility_policy.py'):
  src=HERE/name
  if src.exists():
   digest=hashlib.sha256(src.read_bytes()).hexdigest();target=snapshot/(src.stem+'_'+digest[:8]+src.suffix)
   if not target.exists():shutil.copyfile(src,target)
   sources[str(src)]=dict(sha256=digest,snapshot=str(target))
 (snapshot/'MANIFEST.json').write_text(json.dumps(sources,indent=2),encoding='utf-8');print(json.dumps(dict(registration=str(p),manifest=str(DEST/'RUSTUNA_MANIFEST.json'))))

if __name__=='__main__':main()
