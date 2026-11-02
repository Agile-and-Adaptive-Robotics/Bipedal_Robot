"""Bounded persistent Rustuna adapter; estimator weights never edit mechanics."""
import sys,json,hashlib,gc
from pathlib import Path
import numpy as np
from register_capacity_comparison import DEST,ROOT
import reserve_valve_runner as runner
from capacity_estimator_neural import CapacityReserveForward,SENSOR_CHANNELS
from compact_trial_traces import CompactTraces

def registration():
 cfg=json.loads((DEST/'PREREGISTRATION.json').read_text())
 for path,entry in json.loads((DEST/'source_snapshot/MANIFEST.json').read_text()).items():
  if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Frozen source changed: '+path)
 for spec in cfg['training']+cfg['validation']+[cfg['common_preparation']]:
  if hashlib.sha256(Path(spec['path']).read_bytes()).hexdigest()!=spec['sha256']:raise ValueError('Registered fixture changed')
 return cfg

def validator(params,seed,default=False):
 reference=ROOT/'pressure_reserve_uncertainty/baseline/w2lvar_baseline.release_state.npz' if default else DEST/'instrumentation_preparation/w2lvar_case.release_state.npz'
 def check(path):
  result={}
  with np.load(path) as a,np.load(reference) as b:
   if not set(b.files)<=set(a.files):raise ValueError('Prepared state fields missing')
   for key in b.files:
    if key=='transformed_reserve_parameters':
     expected=b[key].copy();expected[0]=float(params['pressure_cap_kpa']>0);expected[1]=params['preview_s'];expected[5]=1.;result[key]=bool(np.array_equal(a[key],expected))
    elif key=='transformed_valve_pressure_cap_kpa':result[key]=float(a[key])==params['pressure_cap_kpa']
    elif key=='transformed_capacity_parameters':result[key]=bool(np.array_equal(a[key],[params.get('capacity_gain',1.),3,.02,.05,seed]))
    elif key=='transformed_capacity_mode':result[key]=str(a[key])==params['estimator_mode']
    elif key=='transformed_capacity_noise_table_hash':
     rng=np.random.default_rng(seed);noise=np.clip(rng.standard_normal((30000,2,2)),-3.,3.);noise[:,:,0]*=.02;noise[:,:,1]*=.05;result[key]=str(a[key])==hashlib.sha256(noise.tobytes()).hexdigest()
    else:result[key]=bool(np.array_equal(a[key],b[key]))
   if not default and set(a.files)!=set(b.files):raise ValueError('Unexpected extra prepared fields')
  if not all(result.values()):raise ValueError('Common preparation mismatch: '+str([k for k,v in result.items() if not v]))
  return result
 return check

def summarize(c,folder,duration):
 r=json.loads(Path(c['result']).read_text());a=c['assessment'];g=registration()['gates'];mesh=r.get('actual_visual_foot_mesh_clearance_range_m');meshmin=float(np.min(mesh)) if mesh is not None else -1e9;up=r.get('mean_upward_rig_support_bw');down=r.get('mean_downward_rig_force_bw');up=1e9 if up is None else up;down=1e9 if down is None else down
 gates=dict(original_standing=bool(a['accepted']),knees=bool(r.get('max_knee_angle_deg') and max(r['max_knee_angle_deg'])<=0),apparatus_up=up<=g['maximum_mean_upward_rig_bw'],apparatus_down=down<=g['maximum_mean_downward_rig_bw'],mesh=meshmin>=g['minimum_visual_mesh_clearance_m'])
 prefix=Path(c['result'].removesuffix('.result.json'))
 with np.load(prefix.with_suffix('.valve.npz')) as z:delta=z['samples'][:,:,2];mean=float(np.mean(delta)) if len(delta) else 0.;sat=float(np.mean(delta>=.98*c['pressure_cap_kpa'])) if c['pressure_cap_kpa']>0 and len(delta) else 0.
 score=10000.+100.*max(0.,duration-.6-r['completed_samples']*.002) if not all(gates.values()) else float(a.get('sway_m',0)*1000+2*a.get('tilt_deg',0)+200*up+100*down+.25*mean+.5*sat)
 result=dict(accepted=all(gates.values()),gates=gates,objective=score,moving_s=r['completed_samples']*.002,study_duration_s=duration,failure=r['failure'],tilt_deg=a.get('tilt_deg'),sway_m=a.get('sway_m'),mesh_min_m=meshmin,upward_bw=up,downward_bw=down,mean_added_reference_kpa=mean,saturation_fraction=sat,neural_restore=c['complete_neural_restore_verified'],result=c['result'])
 (folder/'CASE_SUMMARY.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def evaluate_case(params,spec,folder,duration=20.,full=False,default=False,prepare_only=False):
 folder.mkdir(parents=True,exist_ok=True);summary=folder/'CASE_SUMMARY.json'
 if summary.exists():return dict(json.loads(summary.read_text()),reused_completed_case=True)
 cfg=registration();created=[];original=runner.ValveReserveForward
 def factory(inner,names,context,cap,preview,target):
  wrapper=CapacityReserveForward(inner,names,context,cap,preview,target,mode=params['estimator_mode'],capacity_gain=params.get('capacity_gain',1.),delay_steps=0 if default else 3,pressure_noise_kpa=0. if default else .02,force_noise_N=0. if default else .05,noise_seed=0 if default else spec['sensor_noise_seed']);created.append(wrapper);return wrapper
 runner.ValveReserveForward=factory
 try:
  with CompactTraces(folder) as traces:
   checker=validator(params,spec['sensor_noise_seed'],default=default) if not prepare_only else validator(params,spec['sensor_noise_seed'],default=True)
   c=runner.run(cap=params['pressure_cap_kpa'],preview=params['preview_s'],condition=spec['name'],encoding='graded',tag='case',duration=duration,reserve_target=1.,output_dir=folder,actuator_configuration=spec['path'],preparation_configuration=cfg['common_preparation']['path'],common_release_validator=checker,sensor_gain=1.,mechanical_push_bw=spec['pulse_bw'])
   wrapper=created[0];np.savez_compressed(folder/'ESTIMATOR_SENSORS.npz',samples=np.asarray(wrapper.sensor_logs),channels=np.array(SENSOR_CHANNELS),routes=np.array(wrapper.routes))
   r=summarize(c,folder,duration)
   if default:
    reference=ROOT/'pressure_reserve_uncertainty/baseline/w2lvar_baseline';parity={}
    for suffix in ('.npz','.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.relay.npz','.valve.npz'):
     with np.load(reference.with_suffix(suffix)) as a,np.load((folder/'w2lvar_case').with_suffix(suffix)) as b:
      assert set(a.files)==set(b.files);parity[suffix]={key:bool(np.array_equal(a[key],b[key])) for key in a.files}
    if not all(v for group in parity.values() for v in group.values()):raise ValueError('Default estimator numerical parity failed')
    (folder/'DEFAULT_PARITY.json').write_text(json.dumps(parity,indent=2),encoding='utf-8')
  retain=full
  registry=DEST/'RETAINED_FAILURES.json';seen=json.loads(registry.read_text()) if registry.exists() else {};mechanism=params['estimator_mode']+'_'+spec['name']
  if r['failure'] and mechanism not in seen and len(seen)<2:retain=True;seen[mechanism]=str(folder);registry.write_text(json.dumps(seen,indent=2),encoding='utf-8')
  traces.materialize(full=retain)
  if prepare_only and not retain:
   path=folder/'w2lvar_case.release_state.npz';traces.original_save(path,**traces.cache[path.resolve()])
  del traces;gc.collect();return r
 finally:runner.ValveReserveForward=original

def preflight():
 gate=DEST/'BASELINE_GATE.json'
 if gate.exists():return json.loads(gate.read_text())
 cfg=registration();spec=cfg['training'][0];params=dict(estimator_mode='pressure',capacity_gain=1.,pressure_cap_kpa=1.,preview_s=.1)
 r=evaluate_case(params,spec,DEST/'default_parity',default=True)
 if not r['accepted']:raise ValueError('Unchanged default nominal failed')
 instrument=evaluate_case(params,spec,DEST/'instrumentation_preparation',duration=.602,prepare_only=True)
 result=dict(ok=True,ready_for_optimization=True,default_baseline=r,default_physical_exact=True,instrumentation_preparation=instrument,scope='Instrumentation0.602s probe only establishes common complete prepared state; not a standing result. All trained/heldout observers use new declared6ms extra sensor transport and noise. Prior default trace remains bitexact.')
 gate.write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def handle_request(request):
 cfg=registration();manifest=request.get('manifest',{})
 if manifest.get('registration_sha256')!=hashlib.sha256((DEST/'PREREGISTRATION.json').read_bytes()).hexdigest():raise ValueError('Registration changed')
 if request['op']=='preflight':return preflight()
 if request['op']!='evaluate':raise ValueError('Unknown operation')
 params=request['params'];mode=params['estimator_mode'];keys={'estimator_mode'}|set(cfg['modes'][mode]['search_space'])
 if set(params)!=keys:raise ValueError('Unknown or unused optimization keys')
 if not json.loads((DEST/'BASELINE_GATE.json').read_text())['ready_for_optimization']:raise ValueError('Baseline gate unmet')
 key=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()[:16];folder=DEST/mode/'training_cases'/key;folder.mkdir(parents=True,exist_ok=True);rows=[]
 for spec in cfg['training']:rows.append(dict(condition=spec['name'],**evaluate_case(params,spec,folder/spec['name'])))
 failed=sum(not row['accepted'] for row in rows);value=10000.*failed+sum(max(0.,19.4-row['moving_s'])*100 for row in rows if not row['accepted'])+max((row['objective'] for row in rows if row['accepted']),default=0.)+sum(row['mean_added_reference_kpa']*.1+row['saturation_fraction']*.2 for row in rows)
 diagnostics=dict(all_training_accepted=failed==0,failed_conditions=failed,conditions=rows,params=params,case_cache=str(folder),evaluation_seed=request['evaluation_seed'],seed_scope='Fixed registered fixture/sensor seeds; optimizer random seed is not an actuator or observation oracle')
 result=dict(ok=True,value=float(value),diagnostics=diagnostics);(folder/'TRIAL_SUMMARY.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

if __name__=='__main__':print(json.dumps(preflight(),indent=2))
