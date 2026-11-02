"""Explicitly gated sequential finalist/validation/valve-only comparison."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
import sys,json,time,subprocess,csv,hashlib,shutil
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;DEST=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\source_domain_recruitment')

def report():
 p=DEST/'pressure/rustuna/TRIALS.jsonl';trials=[json.loads(x) for x in p.read_text().splitlines()] if p.exists() else [];rows=[]
 for t in trials:
  for c in t['diagnostics']['conditions']:rows.append(dict(trial=t['number'],MN_current_cap_nA=t['params']['mn_current_cap_nA'],condition=c['condition'],accepted=c['accepted'],moving_s=c['moving_s'],mean_MN_current_nA=c['mean_added_MN_current_nA'],max_MN_current_nA=c['max_added_MN_current_nA'],pressure_peak_kpa=c['actual_pressure_peak_kpa'],lateral_rig_rms_bw=c['lateral_rig_rms_bw'],lateral_rig_peak_bw=c['lateral_rig_peak_bw'],upward_bw=c['upward_bw'],downward_bw=c['downward_bw'],mesh_min_m=c['mesh_min_m'],failure=c['failure']))
 if rows:
  with (DEST/'TRAINING_COMPARISON.csv').open('w',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 counts={s:sum(r['accepted'] for r in rows if r['condition']==s) for s in ('train_nominal','train_mid','train_strong')};status=dict(machine='easteregg2',completed_trials=len(trials),alltraining_passers=sum(t['diagnostics']['all_training_accepted'] for t in trials),condition_pass_counts=counts)
 for name in ('VALIDATION_STATUS','VALVE_ONLY_CONTROLS'):
  path=DEST/'pressure'/(name+'.json')
  if path.exists():status[name]=json.loads(path.read_text())
 (DEST/'PROGRESS_SUMMARY.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
 text=['# Source-domain recruitment audit and bounded neural comparison','','Machine: **easteregg2**. Virtual-BPA guided human-model standing only. No calibrated CAD hardware, harness-free balance, gait/curriculum, or literature-controller reproduction claim.','','## Offline source/domain result','','All86 configured routes checked at12 observed endpoints. Only pect_r/l have neural targets below geometric minimum. Model route-length reconstruction matches saved failing geometry<1e−10m; equality-aware pulling moment arms agree with local tangent finite differences<1e−6m. Read OFFLINE_RECRUITMENT_AUDIT.json. A latent neural-target gap is NOT automatically extra physical torque: nearest-feasible nominal reference inversion already supplies higher actual tension when requests are infeasible.','','Explicit measured-force probes (not live force overrides) separate hypothetical small5N physical headroom from that latent gap. Strong-case extra force restores overall posture and dissipates mechanical power; lateral COM acceleration worsens relative initialworldreference by~0.0036m/s². Recomputed loaded support margins104–109mm mean this is not proof of support loss. The audit qualifies only a bounded stronglyguided trial, not globalbalance stability. Current-guide and lateral-force constraints remain. Read MEASURED_TENSION_LOAD_PROBE.json.','','## One registered architecture','','Two mirrored6ms leak relays receive predicted geometric minimum plus5N minus measured actualtension, normalized over5N. Prediction uses measuredL+180ms positivevelocity plus6ms transport compensation. Only known assumedrest/fixedlength/K and measuredL,v,T are used; never actual uncertainforcegain. Input current is0..10nA; output is a positive capped conductance current into existing dedicated pectMN input ports, Erev40mV. No original synapse/current cancellation, direct force/voltage projection, or enforced physicalminimum. Source MN fullactivation scale is5mV, not20mV.','','The existing pressure-capacity/valve observer stays frozen atgamma0.75,180ms,10kPa cap. Thus the comparison isolates addedMN recruitment. Only shared MN current cap1–8nA is fitted; no independentleft/rightweights. Scope is pectpair identified by offline data, not allpossible futurepostures. All original model sources and geometry remain unchanged.','','Actual tension drives the error, rather than blindly forcing latent NN requests to a minimum when physical tension is already sufficient. Requested force is logged only as efference-copy diagnostic. Five-newton headroom is a neural error target; actual tension is never assigned or floored.','','Independent source-SNS neuron error0mV, currentconverter error0nA, and32-step fullrecruitment-core state repeat exact. Whole transformed native state restoration is automatically checked upon real recruitment. Zero MN current reproduces inherited valve10 physical arrays/hashes exactly, including actual pressures, forces, neural inputs, mechanics, and apparatus. Source telemetry logging was corrected BEFORE fitting to use the real5mV MN scale; this changed only an inactive zero-gain diagnostic expression, not the verified baseline dynamics.','','## Registered protocol and progress','','Four Rustuna trials on nominal and the same two known uncertainty profiles. Existing620kPa/domain/knee/footload/mesh/vertical apparatus gates remain; added lateralrig RMS≤5%BW and peak≤15%BW, positive cappedMNcurrent gates. Soft current/reference/lateral effort costs, not metabolicenergy. Alltraining acceptance is required before two fresh30s accepted exact-state repeats, then the four untouched carried-forward validation profiles. No heldout reselection. Only after that freshgate are matched gain0 valve-only controls allowed.','','Completed: '+str(len(trials))+'/4; alltraining passers: '+str(status['alltraining_passers'])+'. Nominal/moderate/strong counts: '+str(counts)+'.','','## Files','','PREREGISTRATION.json and RUSTUNA_MANIFEST.json freeze the experiment; BASELINE_GATE.json and INDEPENDENT_RECRUITMENT_CHECKS.json hold exact unit/default evidence. TRAINING_COMPARISON.csv records percondition physics/current/pressure/gates. Search keeps full arrays in RAM and compact windowed witnesses on disk; full evidence only selected/fresh or limited unique failures. PREPARATION_REFERENCE.npz retains the complete shared prepared state. Legacy NPZ act remains selectedexcitation, never nativeactivation.','']
 selection_path=DEST/'pressure/SELECTED.json'
 if selection_path.exists():text+=['','## Frozen result','','Selected sharedMN cap: '+str(json.loads(selection_path.read_text())['params']['mn_current_cap_nA'])+'nA. CoreSNS has888 neurons; two existingvalve and two newMN error-relay equations are separately modeled and unit-verified (892 total neuronal states), not hidden inside the reported888-core census.']
 if 'VALIDATION_STATUS' in status:
  validation=status['VALIDATION_STATUS'];text+=['','Freshnominal30s exact: '+str(validation.get('fresh_exact',False))+'. Independentvalidation:']
  for c in validation.get('validation',[]):text.append('- '+c['condition']+': '+('PASS' if c['accepted'] else 'FAIL')+'; releasedtime '+format(c['moving_s'],'.3f')+'s'+('; '+c['failure'] if c['failure'] else ''))
  text+=['','Frozenvalidation is3/4, with weakerforce-map0.78 failing the unchanged pectineusdomain at18.072s. This rejects generalrobustness; no parameter/model reselection follows. No full-controller±1%basin gate or longstrong30s durability is claimed.']
 if 'VALVE_ONLY_CONTROLS' in status:
  text+=['','MatchedMN-off/valve10 controls:']
  for c in status['VALVE_ONLY_CONTROLS']['rows']:text.append('- '+c['condition']+': '+('PASS' if c['accepted'] else 'FAIL')+'; '+format(c['moving_s'],'.3f')+'s released.')
  text+=['','The active frozen5nA case passes the knownstrong20s study while its matchedMN-off control fails at17.294s. This supports necessity of the ADDED pathway for that particular finite challenge, not necessity of the entireSNS or universal compensation.']
 strong_path=DEST/'STRONG_FRESH_AUDIT.json'
 if strong_path.exists():text+=['','Knownstrong case also passes a separatefreshprocess20s repeat; all original/training full-array hashes match exactly: '+str(json.loads(strong_path.read_text())['all_training_array_hashes_exact'])+'. This adds reproducibility evidence without reselection. MN_RECRUITMENT_STRONG.png/.pdf/.svg shows a short rightleg corrective-feedback window; its pulses are NOT an RG/PF gait rhythm or spiking neuron action potentials.']
 (DEST/'EXPERIMENT_RECORD.md').write_text('\n'.join(text),encoding='utf-8');circuit=dict(core='Frozen normalizedw2lvar standing SNS',existing_valve_branch='Two reserve cells→bounded10kPa reference increment',new_nodes=[dict(name='RECRUIT_pect_'+s,tau_ms=6.,capacitance_nF=6.,leak_uS=1.,input_cap_nA=10.) for s in ('r','l')],new_edges=[dict(source='RECRUIT_pect_'+s,target='MN_pect_'+s,kind='positive conductance, reversal40mV, capped sharedImax') for s in ('r','l')],recruit_log_channels=['geometric_minimum_tension_N','predicted_minimum_tension_N','measured_tension_N','original_MN_requested_tension_N','normalized_error','relay_voltage_mV','added_MN_current_nA','current_cap_nA'],physical_state_projected=False,hidden_synaptic_cancellation=False,oracle_gain=False)
 (DEST/'CIRCUIT.json').write_text(json.dumps(circuit,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in status.items() if k not in ('VALIDATION_STATUS','VALVE_ONLY_CONTROLS')}),flush=True)

def launch(args,name):
 with (DEST/(name+'.log')).open('a',encoding='utf-8') as log:p=subprocess.run(args,stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
 if p.returncode:raise RuntimeError(str(DEST/(name+'.log')))

def repeat_strong():
 # Confirm the selected strong-case trajectory in a separate process: training
 # used a persistent worker, so nominal-only repeats do not settle that issue.
 status=json.loads((DEST/'pressure/VALIDATION_STATUS.json').read_text());fresh=status.get('nominal30s',[])
 if not(status.get('alltraining_candidate') and status.get('fresh_exact') and len(fresh)==2 and all(c['accepted'] and all(c['gates'].values()) for c in fresh)):raise ValueError('Fresh nominal prerequisites required before selected strong repeat')
 import mn_recruitment_experiment as E
 cfg=E.registration();selection=json.loads((DEST/'pressure/SELECTED.json').read_text());spec=next(s for s in cfg['training'] if s['name']=='train_strong');folder=DEST/'pressure/fresh_training_strong';r=E.evaluate_case(selection['params'],spec,folder,duration=20.,full=True);reference=Path(json.loads((DEST/'pressure/rustuna/BEST.json').read_text())['trial']['diagnostics']['case_cache'])/'train_strong/TELEMETRY_MANIFEST.json';checks={}
 # Selection is frozen regardless of this confirmation or validation outcomes.
 for item in json.loads(reference.read_text())['archives']:
  p=folder/Path(item['original_path']).name
  with np.load(p) as z:
   for key,e in item['arrays'].items():checks[p.name+':'+key]=hashlib.sha256(z[key].tobytes()).hexdigest()==e['sha256']
 exact=all(checks.values());(DEST/'STRONG_FRESH_AUDIT.json').write_text(json.dumps(dict(accepted=r['accepted'],all_training_array_hashes_exact=exact,checks=checks,summary=r,scope='Frozen selected knowntraining challenge repeated in freshprocess; no parameter/validation reselection. Full selectedfinalist evidence retained.'),indent=2),encoding='utf-8')
 if not exact:raise ValueError('Fresh strong repeat differs from training; preserve as anti-lottery failure')

def main():
 report();deadline=time.monotonic()+1900.
 while not (DEST/'pressure/rustuna/SUMMARY.json').exists():
  if time.monotonic()>deadline:raise TimeoutError('No second physics process started; optimizer incomplete')
  time.sleep(4.)
 trials=[json.loads(s) for s in (DEST/'pressure/rustuna/TRIALS.jsonl').read_text().splitlines()];eligible=[t for t in trials if t['diagnostics']['all_training_accepted'] and len(t['diagnostics']['conditions'])==3 and all(c['accepted'] and all(c['gates'].values()) for c in t['diagnostics']['conditions'])]
 if len(trials)!=4:raise ValueError('Complete registered trial budget required')
 statuspath=DEST/'pressure/VALIDATION_STATUS.json'
 if not eligible:statuspath.write_text(json.dumps(dict(mode='pressure',alltraining_candidate=False,new_validation_started=False,scope='No alltraining passer; no fresh/validation/control process launched.'),indent=2),encoding='utf-8')
 elif not statuspath.exists():launch([sys.executable,str(HERE/'mn_recruitment_experiment.py'),'--validate'],'validation')
 status=json.loads(statuspath.read_text());fresh=status.get('nominal30s',[]);controls_allowed=bool(eligible and status.get('alltraining_candidate') and status.get('fresh_exact') and len(fresh)==2 and all(c['accepted'] and all(c['gates'].values()) for c in fresh));controlpath=DEST/'pressure/VALVE_ONLY_CONTROLS.json'
 if controls_allowed and not controlpath.exists():
  from mn_recruitment_experiment import child,registration
  cfg=registration();rows=[]
  for s in cfg['training']:
   if s['name']=='train_nominal':r=dict(json.loads((DEST/'zero_default_parity/CASE_SUMMARY.json').read_text()),reused_completed_exact_baseline=True)
   else:r=child('off_'+s['name'])
   rows.append(dict(condition=s['name'],**r))
  controlpath.write_text(json.dumps(dict(rows=rows,scope='MN relay gain0; valveobserver10 and original feedback remain live. Same complete preparation. Do not claim necessity when live and control bothpass.'),indent=2),encoding='utf-8')
 if controls_allowed and not (DEST/'STRONG_FRESH_AUDIT.json').exists():launch([sys.executable,str(Path(__file__)),'--repeat-strong'],'strong_fresh')
 report();(DEST/'COMPLETE.json').write_text(json.dumps(dict(completed=True,prerequisites_explicitly_enforced=True),indent=2),encoding='utf-8')

if __name__=='__main__':
 if '--report' in sys.argv:report()
 elif '--repeat-strong' in sys.argv:repeat_strong()
 else:main()
