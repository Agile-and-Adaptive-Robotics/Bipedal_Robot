"""Persistent Rustuna scientific adapter with compact exploratory telemetry."""
import sys,json,hashlib,gc
from pathlib import Path
import numpy as np
from register_reserve_uncertainty import DEST,ROOT
import reserve_valve_runner as runner
from compact_trial_traces import CompactTraces

KEYS={'pressure_cap_kpa','preview_s','sensor_gain'}

def registration():
 p=DEST/'PREREGISTRATION.json';cfg=json.loads(p.read_text())
 snapshot=DEST/'source_snapshot/MANIFEST.json'
 if snapshot.exists():
  for path,entry in json.loads(snapshot.read_text()).items():
   if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Experiment source changed: '+path)
 for f in cfg['fixtures']:
  if hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()!=f['sha256']:raise ValueError('Registered fixture changed')
 return cfg

def validator(params,initial=False):
 reference=(ROOT/'pressure_reserve_valve/w2lvar_winner_fresh.release_state.npz') if initial else (DEST/'baseline/w2lvar_baseline.release_state.npz')
 def check(path):
  result={}
  with np.load(path) as a,np.load(reference) as b:
   if set(a.files)!=set(b.files):raise ValueError('Prepared state fields differ')
   for key in b.files:
    if initial and key=='pressure_history':
     ca=int(a['pressure_step_count']);cb=int(b['pressure_step_count']);na=a[key].shape[1];nb=b[key].shape[1]
     result[key+'_chronological_overlap']=bool(np.array_equal(a[key][:,[(ca-i)%na for i in range(1,nb+1)]],b[key][:,[(cb-i)%nb for i in range(1,nb+1)]]))
    elif initial and key=='pressure_step_count':result[key+'_minus_capacity']=int(a[key])-a['pressure_history'].shape[1]==int(b[key])-b['pressure_history'].shape[1]
    elif key=='transformed_reserve_parameters':
     expected=b[key].copy();expected[0]=float(params['pressure_cap_kpa']>0);expected[1]=params['preview_s'];expected[5]=1./params['sensor_gain'];result[key+'_only_declared_parameters']=bool(np.array_equal(a[key],expected))
    elif key=='transformed_valve_pressure_cap_kpa':result[key+'_declared_parameter']=float(a[key])==params['pressure_cap_kpa']
    else:result[key]=bool(np.array_equal(a[key],b[key]))
  if not all(result.values()):raise ValueError('Preparation mismatch: '+str([k for k,v in result.items() if not v]))
  return result
 return check

def summarize(c,root):
 r=json.loads(Path(c['result']).read_text());a=c['assessment'];cfg=registration();g=cfg['gates'];mesh=r.get('actual_visual_foot_mesh_clearance_range_m');meshmin=float(np.min(mesh)) if mesh is not None else -1e9
 up=r.get('mean_upward_rig_support_bw');down=r.get('mean_downward_rig_force_bw');up=1e9 if up is None else up;down=1e9 if down is None else down
 gates=dict(original_standing=bool(a['accepted']),knees=bool(r.get('max_knee_angle_deg') and max(r['max_knee_angle_deg'])<=0),apparatus_up=up<=g['maximum_mean_upward_rig_bw'],apparatus_down=down<=g['maximum_mean_downward_rig_bw'],mesh=meshmin>=g['minimum_visual_mesh_clearance_m'])
 prefix=Path(c['result'].removesuffix('.result.json'))
 with np.load(prefix.with_suffix('.valve.npz')) as z:vals=z['samples'];delta=vals[:,:,2];mean=float(np.mean(delta)) if len(delta) else 0.;sat=float(np.mean(delta>=.98*c['pressure_cap_kpa'])) if c['pressure_cap_kpa']>0 and len(delta) else 0.
 score=10000.+100.*max(0.,19.4-r['completed_samples']*.002) if not all(gates.values()) else float(a.get('sway_m',0)*1000+2*a.get('tilt_deg',0)+200*up+100*down+.25*mean+.5*sat)
 result=dict(accepted=all(gates.values()),gates=gates,objective=score,moving_s=r['completed_samples']*.002,failure=r['failure'],tilt_deg=a.get('tilt_deg'),sway_m=a.get('sway_m'),mesh_min_m=meshmin,upward_bw=up,downward_bw=down,mean_added_reference_kpa=mean,saturation_fraction=sat,neural_restore=c['complete_neural_restore_verified'],result=c['result'],comparison=str(prefix.with_suffix('.comparison.json')))
 (root/'CASE_SUMMARY.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def evaluate_case(params,spec,folder,full=False,initial=False):
 folder.mkdir(parents=True,exist_ok=True);summary=folder/'CASE_SUMMARY.json'
 if summary.exists():r=json.loads(summary.read_text());r['reused_completed_case']=True;return r
 config=registration();kwargs=dict(cap=params['pressure_cap_kpa'],preview=params['preview_s'],condition=spec['name'],encoding='graded',tag='baseline' if initial else 'case',duration=20.,reserve_target=1.,output_dir=folder,actuator_configuration=spec['path'],preparation_configuration=config['common_preparation']['path'],common_release_validator=validator(params,initial),sensor_gain=params['sensor_gain'],mechanical_push_bw=spec['pulse_bw'])
 if full:c=runner.run(**kwargs);r=summarize(c,folder)
 else:
  with CompactTraces(folder) as traces:c=runner.run(**kwargs);r=summarize(c,folder)
  # Only the first instance of a causal route/condition failure retains full
  # exploratory arrays. Every trial still retains exact gate/scalar summaries.
  mechanism='none'
  if r['failure']:
   mechanism=spec['name']+'_'+('inverse_domain' if 'declared contraction domain' in r['failure'] else 'other')
   if 'name\': \'' in r['failure']:mechanism+='_'+r['failure'].split("name': '",1)[1].split("'",1)[0]
  registry=DEST/'RETAINED_CAUSAL_FAILURES.json';seen=json.loads(registry.read_text()) if registry.exists() else {}
  retain=bool(r['failure'] and mechanism not in seen and sum(v.get('full_trace_retained',False) for v in seen.values())<3)
  traces.materialize(full=retain)
  if r['failure'] and mechanism not in seen:seen[mechanism]=dict(folder=str(folder),full_trace_retained=retain);registry.write_text(json.dumps(seen,indent=2),encoding='utf-8')
  del traces;gc.collect()
 return r

def preflight():
 cfg=registration();params=dict(pressure_cap_kpa=1.,preview_s=.1,sensor_gain=1.);spec=next(f for f in cfg['fixtures'] if f['name']=='train_nominal');r=evaluate_case(params,spec,DEST/'baseline',full=True,initial=True)
 old=ROOT/'pressure_reserve_valve/w2lvar_winner_fresh';new=DEST/'baseline/w2lvar_baseline';parity={}
 for suffix in ('.npz','.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.relay.npz','.valve.npz'):
  with np.load(old.with_suffix(suffix)) as a,np.load(new.with_suffix(suffix)) as b:
   assert set(a.files)==set(b.files);parity[suffix]={k:bool(np.array_equal(a[k],b[k])) for k in a.files}
 exact=all(v for group in parity.values() for v in group.values());io_path=DEST/'COMPACT_NOMINAL_GATE.json'
 if io_path.exists():io=json.loads(io_path.read_text())
 else:
  folder=DEST/'compact_nominal_gate';folder.mkdir(exist_ok=True)
  with CompactTraces(folder) as traces:
   c=runner.run(cap=1.,preview=.1,condition='train_nominal',encoding='graded',tag='case',duration=20.,reserve_target=1.,output_dir=folder,actuator_configuration=spec['path'],preparation_configuration=cfg['common_preparation']['path'],common_release_validator=validator(params),sensor_gain=1.,mechanical_push_bw=0.);summary=summarize(c,folder)
  arrays={}
  for path,content in traces.cache.items():
   suffix=path.name.removeprefix('w2lvar_case');reference=new.with_suffix(suffix)
   with np.load(reference) as z:
    assert set(content)==set(z.files);arrays[suffix]={key:bool(np.array_equal(value,z[key])) for key,value in content.items()}
  passed=all(v for group in arrays.values() for v in group.values());traces.materialize(full=False);io=dict(exact=passed,all_original_gate_arrays=arrays,summary=summary);io_path.write_text(json.dumps(io,indent=2),encoding='utf-8');del traces;gc.collect()
 report=dict(ok=True,ready_for_optimization=bool(r['accepted'] and exact and io['exact']),baseline=r,nominal_extended_history_physical_exact=parity,compact_RAM_gate_physical_exact=io,scope='24ms history capacity only; unchanged nominal pressure/NN/physics trajectories. Training fixtures are now exposed; registered validation seeds remain unrun.')
 (DEST/'BASELINE_GATE.json').write_text(json.dumps(report,indent=2),encoding='utf-8');return report

def handle_request(request):
 cfg=registration();manifest=request.get('manifest',{})
 if manifest.get('registration_sha256')!=hashlib.sha256((DEST/'PREREGISTRATION.json').read_bytes()).hexdigest():raise ValueError('Registration changed')
 if request['op']=='preflight':return preflight()
 if request['op']!='evaluate' or set(request['params'])!=KEYS:raise ValueError('Unknown/unused optimization keys')
 params=request['params'];gate=json.loads((DEST/'BASELINE_GATE.json').read_text())
 if not gate['ready_for_optimization']:raise ValueError('Nominal baseline/parity gate not met')
 key=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()[:16];root=DEST/'training_cases'/key;root.mkdir(parents=True,exist_ok=True);rows=[]
 for spec in cfg['fixtures']:
  if spec['split']!='train':continue
  if spec['name']=='train_nominal' and params==dict(pressure_cap_kpa=1.,preview_s=.1,sensor_gain=1.):r=dict(gate['baseline'],reused_completed_case=True)
  else:r=evaluate_case(params,spec,root/spec['name'])
  rows.append(dict(condition=spec['name'],**r))
 failed=sum(not r['accepted'] for r in rows)
 value=10000.*failed+sum((19.4-r['moving_s'])*100 for r in rows if not r['accepted'])+max((r['objective'] for r in rows if r['accepted']),default=0.)+sum(r['mean_added_reference_kpa']*.1+r['saturation_fraction']*.2 for r in rows)
 diagnostics=dict(all_training_accepted=failed==0,failed_conditions=failed,conditions=rows,params=params,case_cache=str(root),evaluation_seed=request['evaluation_seed'],seed_scope='Optimizer seed only; fixed paired-route fixture seeds are preregistered and unchanged across candidates')
 (root/'TRIAL_SUMMARY.json').write_text(json.dumps(dict(value=value,diagnostics=diagnostics),indent=2),encoding='utf-8');return dict(ok=True,value=float(value),diagnostics=diagnostics)
