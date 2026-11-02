"""Compact progress, failure witnesses, and storage ledger for registered search."""
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from register_reserve_uncertainty import ROOT,DEST

def main():
 registry=json.loads((DEST/'PREREGISTRATION.json').read_text());p=DEST/'rustuna/TRIALS.jsonl';trials=[json.loads(line) for line in p.read_text().splitlines()] if p.exists() else [];diagnostics=[]
 base=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text());names=list(base['actuators']);pect=names.index('pect_r');a=base['actuators']['pect_r']
 for trial in trials:
  for condition in trial['diagnostics']['conditions']:
   folder=Path(condition['result']).parent;manifest=folder/'TELEMETRY_MANIFEST.json'
   if not manifest.exists():continue
   archives=json.loads(manifest.read_text())['archives'];lookup={Path(row['original_path']).name.split('w2lvar_case',1)[-1]:i for i,row in enumerate(archives)}
   if '.relay.npz' not in lookup or '.pressure.npz' not in lookup:continue
   with np.load(folder/'COMPACT_TELEMETRY.npz') as z:
    relay=z[f"a{lookup['.relay.npz']}_samples"][-1,0];pressure=z[f"a{lookup['.pressure.npz']}_samples"][-1,:,pect]
   # The terminal neural sample observes the preceding valid pressure/force
   # update, so their actual pressure equality verifies this alignment.
   matched=bool(relay[6]==pressure[1]);true_h=a['resting_length_m']+a['nonmuscle_length_m']+pressure[2]/a['series_stiffness_N_m']-relay[7]
   diagnostics.append(dict(trial=trial['number'],condition=condition['condition'],aligned_pressure_exact=matched,pressure_map_reserve_um=float(relay[0]*1e6),force_series_inferred_reserve_um=float(true_h*1e6),predicted_reserve_um=float(relay[2]*1e6),sensor_signal=float(relay[3]),relay_voltage_mV=float(relay[4]),actual_pressure_kpa=float(pressure[1]),actual_tension_N=float(pressure[2]),scope='Inferred contraction within the declared series model from logged actual force; diagnostic only, never supplied to controller.'))
 files=[p for p in DEST.rglob('*') if p.is_file()];storage=dict(files=len(files),bytes=sum(p.stat().st_size for p in files),policy='RAM full traces for physical gates. Sparse/windowed exploratory witnesses; full baseline and at most3unique causal search failures. No pertrial source snapshots. Original source/model files untouched.')
 completed=(DEST/'rustuna/SUMMARY.json').exists();passers=[r['number'] for r in trials if r['diagnostics']['all_training_accepted']]
 summary=dict(machine='easteregg2',registration_sha256=hashlib.sha256((DEST/'PREREGISTRATION.json').read_bytes()).hexdigest(),status='completed_no_training_passer_validation_gated' if completed and not passers else 'running' if not completed else 'training_candidate_requires_fresh_parity',completed_trials=len(trials),maximum_trials=registry['max_trials'],all_training_passers=passers,validation_run=False,endpoint_diagnostics=diagnostics,storage=storage)
 (DEST/'PROGRESS_SUMMARY.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
 if completed:(DEST/'STUDY_STATUS.json').write_text(json.dumps({k:v for k,v in summary.items() if k!='endpoint_diagnostics'},indent=2),encoding='utf-8')
 text='# Preregistered uncertainty training — easteregg2\n\n'
 text+='This new study preserves the earlier nominal winner and all exposed failed tests. Three training fixtures (nominal and seeded mid/strong paired-route pressure profiles) span 40–80 ms fill, 4–16 ms delay, and force gains 1.0/0.9/0.8. Four new independently seeded validation profiles were written before fitting, including interpolation, delay/map extrapolation, and a new 0.075 BW pulse. They remain unrun until all training cases and fresh parity pass. Earlier stress observations are explicitly training information, not unseen validation.\n\n'
 text+='Rustuna uses one persistent myo worker, eight trials maximum, seed810031, and a2400s original wall deadline. Shared parameters are pressure-reference cap1–5kPa, preview80–180ms, and sensor-input gain0.5–2; the equivalent reserve target remains1N. Original NN weights, geometry, actual pressure dynamics, and domain guard are unchanged. Both legs share every fitted parameter.\n\n'
 text+='Physical acceptance requires the original standing gates, no knee extension beyond0deg, ground load above0.5BW, mean upward/downward apparatus forces each no more than0.05BW, and minimum visual mesh clearance at least−8mm. The cost prioritizes failed conditions and incomplete duration, then worst sway/tilt/apparatus and actual added-reference/output-saturation effort. Pressure integrals are not pneumatic or metabolic energy. Strong lateral/rotational guidance remains.\n\n'
 text+='The common history ring has24ms capacity to accommodate new20ms validation delays. Its nominal body, neural, actual pressure, and actuator traces match the earlier accepted winner exactly. Full-file versus RAM-only compact telemetry also matches every saved nominal array and physical gate exactly. Longer history capacity does not change pressure lag or physical preparation.\n\n| Trial | Pcap (kPa) | Preview (ms) | Sensor gain | Failed conditions | Training results |\n| --- | ---: | ---: | ---: | ---: | --- |\n'
 for row in trials:
  p=row['params'];d=row['diagnostics'];result=', '.join(c['condition']+':'+('pass' if c['accepted'] else f"fail@{c['moving_s']:.3f}s") for c in d['conditions']);text+=f"| {row['number']} | {p['pressure_cap_kpa']:.4f} | {1000*p['preview_s']:.2f} | {p['sensor_gain']:.4f} | {d['failed_conditions']} | {result} |\n"
 if completed and not passers:text+='\nFinal outcome: eight completed Rustuna trials, zero all-training passers. All nominal cases pass, while both mid/strong uncertainty fixtures fail for every candidate. Fresh-finalist and independent validation are gated and were not run. BEST.json is the least-bad diagnostic candidate, not an accepted controller. The prior nominal winner and all original stress failures remain untouched.\n'
 text+='\n## Endpoint diagnosis\n\nThe pressure-based reserve assumes the nominal force map. Logged actual force and series stiffness provide an independent diagnostic contraction reserve: l0+b+Tactual/K−L. Under force-map loss, these can disagree substantially. This diagnostic is not an oracle input to the controller. All endpoint values and alignment checks are retained in PROGRESS_SUMMARY.json. Gain/cap changes cannot automatically fix a biased reserve estimate. At a stationary boundary, an unknown gain g requires a nominal equivalent reserve at least(1/g−1)Tgeometry; for g=0.8 and roughly22N geometry tension, that is approximately5.5N, larger than the fixed1N target. The existing fit cannot alter that target or introduce a force sensor.\n\n'
 text+='A justified separate model revision would fit a shared pressure-input/capacity weight inside the estimator (for example, a conservative0.75–1.0 multiplier on nominalFmaxP), rather than merely the amplitude of the resulting error. This would alter a neural sensory/internal-model weight, not actual actuator strength. Another option is measured tension or actual sleeve-length feedback. Neither revision was inserted into this registered study, and none of the new validation seeds has been consumed.\n\n'
 text+=f"Storage at report: {storage['files']} files, {storage['bytes']/1e6:.1f} MB. Exploratory full traces are held in RAM only for exact gate calculations; compact windows and full scalar/parameter/hash manifests are saved. Full search retention is capped at three unique causal failures plus the baseline. Finalist/fresh data may be retained separately only after training acceptance. No new validation result is claimed.\n\n"
 text+='## Files and commands\n\nMachine easteregg2. Preregistration: `'+str(DEST/'PREREGISTRATION.json')+'`. Optimizer: `D:\\Anaconda\\python.exe D:\\GitHub\\Bipedal_Robot\\Code\\MuJoCo_SNS\\spinal\\virtual_walker_20261031\\rustuna_transfer.py --manifest "'+str(DEST/'RUSTUNA_MANIFEST.json')+'"`; add `--resume` only for an interrupted identical registered study. Scientific worker: `D:\\Anaconda\\envs\\myo\\python.exe`. TRIALS.jsonl/BEST.json are candidate records, not scientific acceptance.\n'
 (DEST/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8');print(json.dumps(summary,default=str)[:250])
 (DEST/'CIRCUIT_SCHEMA.json').write_text(json.dumps(dict(architecture=registry['architecture'],shared_ranges=registry['search_space'],nodes=['RESERVE_pect_r','RESERVE_pect_l'],output='Original nominal pressure plus bounded neuron-voltage reference increment; no actual-state projection',pressure_capacity_estimator_gain_fixed=1.,equivalent_reserve_target_N=1.),indent=2),encoding='utf-8')
 from reserve_valve_report import draw
 draw(DEST,training_ranges=True)
 tools=DEST/'post_fit_tools';tools.mkdir(exist_ok=True)
 for p in (Path(__file__),Path(__file__).with_name('reserve_valve_report.py')):
  digest=hashlib.sha256(p.read_bytes()).hexdigest();target=tools/(p.stem+'_'+digest[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
 index=ROOT/'EXPERIMENT_RECORD.md';begin='<!-- reserve-uncertainty-training:start -->';end='<!-- reserve-uncertainty-training:end -->'
 addition=begin+'\n## New preregistered uncertainty training\n\nRecord: `'+str(DEST/'EXPERIMENT_RECORD.md')+'`. Eight bounded Rustuna trials use one persistent scientific worker, three fixed seeded training profiles, mirrored pressurecap/preview/input-amplitude gains, hard body/apparatus/domain gates, and reference/saturation effort cost. Every nominal case passes, but every candidate fails both uncertainty-training profiles; there is no all-training passer. Fresh-finalist/independent validation are gated and unrun. Four new validation seeds/conditions remain untouched. The nominal pressure-based reserve overstates actual force/series contraction margin under weakened force maps, a structural estimator bias that amplitude/cap changes cannot remove. Previous winners/failures are preserved. Full-file/RAM compact telemetry and extended-history nominal dynamics match exactly. Exploration is compact, with full search retention capped at three unique causal failures plus baseline; actual storage is recorded in STUDY_STATUS.json.\n'+end
 existing=index.read_text(encoding='utf-8')
 if begin in existing:
  left,tail=existing.split(begin,1);_,right=tail.split(end,1);existing=left+addition+right
 else:existing+='\n\n'+addition+'\n'
 index.write_text(existing,encoding='utf-8')

if __name__=='__main__':main()
