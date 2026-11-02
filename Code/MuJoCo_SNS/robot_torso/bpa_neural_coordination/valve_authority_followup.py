"""Sequential fresh/validation/causal controls, with compact durable reporting."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
import sys,json,time,subprocess,csv,hashlib,shutil
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;DEST=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\valve_authority')

def report():
 p=DEST/'pressure/rustuna/TRIALS.jsonl';trials=[json.loads(x) for x in p.read_text().splitlines()] if p.exists() else [];rows=[]
 for t in trials:
  for c in t['diagnostics']['conditions']:rows.append(dict(trial=t['number'],cap_kpa=t['params']['pressure_cap_kpa'],condition=c['condition'],accepted=c['accepted'],moving_s=c['moving_s'],peak_actual_pressure_kpa=c['actual_pressure_peak_kpa'],upward_bw=c['upward_bw'],downward_bw=c['downward_bw'],mesh_min_m=c['mesh_min_m'],mean_added_reference_kpa=c['mean_added_reference_kpa'],saturation_fraction=c['saturation_fraction'],failure=c['failure']))
 if rows:
  with (DEST/'TRAINING_COMPARISON.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 counts={s:sum(r['accepted'] for r in rows if r['condition']==s) for s in ('train_nominal','train_mid','train_strong')};status=dict(machine='easteregg2',completed_trials=len(trials),alltraining_passers=sum(t['diagnostics']['all_training_accepted'] for t in trials),condition_pass_counts=counts)
 summary_path=DEST/'pressure/rustuna/SUMMARY.json'
 if summary_path.exists():status['first_study_process_wall_s']=json.loads(summary_path.read_text())['this_process_wall_s']
 v=DEST/'pressure/VALIDATION_STATUS.json'
 if v.exists():status['validation']=json.loads(v.read_text())
 c=DEST/'CAUSAL_CONTROL_AUDIT.json'
 if c.exists():status['causal_controls']=json.loads(c.read_text())
 (DEST/'PROGRESS_SUMMARY.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
 text=['# Frozen-estimator neural valve-authority experiment','','Machine: **easteregg2**. New experimental variant, not a replacement for prior checkpoints. Gait2392 virtual-BPA guided standing only; no CAD hardware or gait/curriculum claim.','','The learned pressure-capacity estimator is frozen at sharedgamma0.75 and180ms preview, with the previous2 mirrored6ms leak neurons,1N-equivalent target,6ms additional sensor delay, and seeded clipped pressure/tension noise. Only the shared voltage-to-pressure-reference increment cap is fitted over6–10kPa. PhysicalFmax,620kPa ceiling,40/60ms nominal pressure response, declared uncertainty profiles, rate/lag laws, geometry, and contraction-domain guard remain unchanged. No actual-pressure projection or hidden force floor.','','Common nominal pinned preparation is identical before faults switch at the first released physics step. Noise/filter/transport states and original neural/mechanical/pressure histories are included in common-state and exact restoration audits. Previous models/source files remain frozen.','','## Registered gates and evidence','','Four Rustuna trials, same three known training profiles, alltraining acceptance before two fresh30s nominal runs and all numeric/state equality; only then the carried-forward four untouched validation profiles. Selection is frozen before validation. Domain/knee/actual footload/apparatus/visual mesh gates unchanged; actual pressure ceiling and native force parity added explicitly. Reference-increment and saturation effort penalized; this is not metabolic or pneumatic energy.','','Baseline gate: inherited cap5 path must match prior physical full-array hashes. Source snapshot and registration hashes are saved before fitting. Search keeps compact windowed telemetry/full-array hashes; full traces only parity/finalists and at most two unique causal failures.','','## Current results','',f"Completed trials: {len(trials)}/4. Alltraining passers: {status['alltraining_passers']}.",f"Nominal/moderate/strong pass counts: {counts['train_nominal']}/{counts['train_mid']}/{counts['train_strong']}.",'','The model retains strong lateral/rotational guidance. Mean upward/downward apparatus forces, actual visual-foot penetration, and actual pressure are shown in TRAINING_COMPARISON.csv. Native activation arrays are unused bookkeeping for pressure actuators; excitation is never mislabeled activation.']
 if v.exists():text+=['',f"Independent validation started: {status['validation']['new_validation_started']}. {status['validation']['scope']}"]
 attribution=DEST/'FAILURE_ATTRIBUTION.json'
 if attribution.exists():
  d=json.loads(attribution.read_text());text+=['','## Failure attribution','','At10kPa, the last exact neural excitation requests17.39N pectineus tension, while the declared series geometry already requires67.94N to keep nonnegative sleeve contraction. The additional reserve neuron is saturated at10mV/full10kPa. Offline static weak-map correction consumes approximately8.81kPa of that increment, leaving1.19kPa for lag/delay. This explains a growing posture/recruitment versus virtual actuator-domain mismatch; increasing a fixed offset delays it but does not solve it. True mapgain is used only for offline attribution, never controller input. These are model-domain aborts while supported, not observed falls or hardware infeasibility claims. Virtual sleeve dimensions and400kN/m pair stiffness remain assumed. See FAILURE_ATTRIBUTION.json for allfour cases.']
 if v.exists() and not status['validation'].get('alltraining_candidate'):text+=['','Fresh30s, independent validation, nominal-replay, and correction-off controls were not executed because no candidate passed alltraining. No new necessity or general robustness claim is made. All carried-forward validation profiles remain unused.']
 text+=['','## Causal-control interpretation','','After fresh stability only: replay the nominal reserve-neuron input time course while keeping original neural feedback live; compare with live reserve feedback and correction-off under identical complete preparation. This isolates necessity of the added reserve pathway, not the entire SNS. Nominal replay must reproduce the nominal physical prefix exactly. Both passing does not demonstrate necessity. No postvalidation parameter or mode reselection.','','Main source: `'+str(HERE/'valve_authority_experiment.py')+'`. Registration: `'+str(DEST/'PREREGISTRATION.json')+'`.']
 (DEST/'EXPERIMENT_RECORD.md').write_text('\n'.join(text)+'\n',encoding='utf-8');print(json.dumps({k:v for k,v in status.items() if k not in ('validation','causal_controls')}),flush=True)
 protocol={};snap=DEST/'source_snapshot'
 for name in ('valve_authority_followup.py','valve_authority_diagnostics.py'):
  source=HERE/name;digest=hashlib.sha256(source.read_bytes()).hexdigest();target=snap/(source.stem+'_'+digest[:8]+source.suffix)
  if not target.exists():shutil.copyfile(source,target)
  protocol[str(source)]=dict(sha256=digest,snapshot=str(target))
 (snap/'PROTOCOL_MANIFEST.json').write_text(json.dumps(protocol,indent=2),encoding='utf-8')

def launch(args,name):
 with (DEST/(name+'.log')).open('a',encoding='utf-8') as log:p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
 if p.returncode:raise RuntimeError(str(DEST/(name+'.log')))

def controls():
 status=json.loads((DEST/'pressure/VALIDATION_STATUS.json').read_text())
 if not status.get('fresh_exact'):return
 from valve_authority_experiment import child,registration
 cfg=registration();rows=[]
 for s in cfg['training']:
  for mode in ('replay','zero'):rows.append(dict(condition=s['name'],mode=mode,**child(mode+'_'+s['name'])))
 checks={};a=DEST/'pressure/fresh_nominal_a/w2lvar_case';b=DEST/'pressure/replay_train_nominal/w2lvar_case'
 # Both include the same300 preparation neural ticks; compare time-series prefix.
 for suffix in ('.npz','.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.relay.npz','.valve.npz'):
  with np.load(a.with_suffix(suffix)) as x,np.load(b.with_suffix(suffix)) as y:
   checks[suffix]={}
   for key in y.files:
    u=x[key];v=y[key]
    if u.shape!=v.shape and u.ndim and v.ndim and u.shape[1:]==v.shape[1:] and len(u)>=len(v):u=u[:len(v)]
    checks[suffix][key]=bool(np.array_equal(u,v))
 nominal_exact=all(v for group in checks.values() for v in group.values());audit=dict(controls=rows,nominal_replay_physical_prefix_exact=nominal_exact,nominal_prefix_checks=checks,necessity_by_condition={s['name']:dict(live_accepted=True,replay_accepted=next(r['accepted'] for r in rows if r['condition']==s['name'] and r['mode']=='replay'),zero_accepted=next(r['accepted'] for r in rows if r['condition']==s['name'] and r['mode']=='zero')) for s in cfg['training']},scope='No necessity claim if both live and control pass; only additional reserve pathway tested. Nominal replay physical prefix parity required before interpretation.')
 (DEST/'CAUSAL_CONTROL_AUDIT.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
 if not nominal_exact:raise ValueError('Nominal reserve-input replay does not reproduce physical prefix')

def main():
 report();deadline=time.monotonic()+1900.
 while not (DEST/'pressure/rustuna/SUMMARY.json').exists():
  if time.monotonic()>deadline:raise TimeoutError('No second physics worker started; training did not finish')
  time.sleep(4.)
 report()
 # Explicit launch gate: absence of a status file never authorizes validation.
 trials=[json.loads(line) for line in (DEST/'pressure/rustuna/TRIALS.jsonl').read_text().splitlines()]
 eligible=[t for t in trials if t['diagnostics']['all_training_accepted'] and len(t['diagnostics']['conditions'])==3 and all(c['accepted'] and all(c['gates'].values()) for c in t['diagnostics']['conditions'])]
 if len(trials)!=4:raise ValueError('Registered four-trial budget incomplete; no downstream scientific process launched')
 if not eligible:
  (DEST/'pressure/VALIDATION_STATUS.json').write_text(json.dumps(dict(mode='pressure',alltraining_candidate=False,new_validation_started=False,scope='Explicit prelaunch gate: no alltraining passer, so no fresh30s or independentvalidation process is launched.'),indent=2),encoding='utf-8')
 elif not (DEST/'pressure/VALIDATION_STATUS.json').exists():launch([sys.executable,str(HERE/'valve_authority_experiment.py'),'--validate'],'validation')
 report()
 downstream=json.loads((DEST/'pressure/VALIDATION_STATUS.json').read_text())
 fresh_rows=downstream.get('nominal30s',[])
 controls_authorized=bool(eligible and downstream.get('alltraining_candidate') and downstream.get('fresh_exact') and len(fresh_rows)==2 and all(row['accepted'] and all(row['gates'].values()) for row in fresh_rows))
 if controls_authorized and not (DEST/'CAUSAL_CONTROL_AUDIT.json').exists():controls()
 report();(DEST/'COMPLETE.json').write_text(json.dumps(dict(completed=True,scope='Every downstream stage obeys registered prerequisite gates'),indent=2),encoding='utf-8')

if __name__=='__main__':report() if '--report' in sys.argv else main()
