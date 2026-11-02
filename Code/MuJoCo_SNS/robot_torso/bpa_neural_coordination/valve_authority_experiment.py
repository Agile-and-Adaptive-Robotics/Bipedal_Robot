"""Distinct bounded output-authority study; original estimator/physics frozen."""
import os
for key in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[key]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
import json,sys,hashlib,gc,argparse,subprocess
from pathlib import Path
import numpy as np
import capacity_comparison_adapter as base
from capacity_estimator_neural import CapacityReserveForward
from compact_trial_traces import CompactTraces
from register_capacity_comparison import ROOT
OLD=ROOT/'pressure_reserve_capacity';DEST=ROOT/'valve_authority';HERE=Path(__file__).resolve().parent
ORIGINAL_SUMMARIZE=base.summarize
CONTROL='live';REPLAY=None

class AuthorityForward(CapacityReserveForward):
 def __init__(self,*args,**kwargs):
  args=list(args);cap=float(args[3])
  if not (cap==0. or 0<cap<=10.):raise ValueError('New declared neural output cap0..10kPa only')
  args[3]=min(cap,5.);super().__init__(*args,**kwargs);self.pressure_cap=cap;self.control_mode=CONTROL;self.replay=None if REPLAY is None else REPLAY.copy();self.replay_hash='none' if self.replay is None else hashlib.sha256(self.replay.tobytes()).hexdigest();self.replay_count=0
 def _step(self,x):
  index=self.steps;old=self.relay_v.copy();result=super()._step(x)
  if self.control_mode=='replay' and self.context.get('active',False) and self.pressure_cap>0:
   if self.replay is None or index>=len(self.replay):raise ValueError('No replay padding/extrapolation')
   signal=self.replay[index];self.relay_v[:]=old+self.dt/self.tau*(10.*signal-old);self.logs[-1][:,3]=signal;self.replay_count+=1
  return result
 def export_state(self):
  d=super().export_state();d.update(authority_control_mode=np.array(self.control_mode),authority_replay_hash=np.array(self.replay_hash),authority_replay_count=np.array(self.replay_count));return d
 def import_state(self,d):
  assert str(d['authority_control_mode'])==self.control_mode and str(d['authority_replay_hash'])==self.replay_hash
  super().import_state(d);self.replay_count=int(d['authority_replay_count'])
 def circuit(self):
  d=super().circuit();d.update(new_output_authority_only=True,neural_reference_cap_kpa=self.pressure_cap,reserve_control=self.control_mode,replay_scope='Only reserve-neuron normalized input replayed from frozen nominal. Original neural feedback stays live. Actual-pressure lag/physics unchanged.');return d

def registration():
 cfg=json.loads((DEST/'PREREGISTRATION.json').read_text())
 for path,entry in json.loads((DEST/'source_snapshot/MANIFEST.json').read_text()).items():
  if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Registered source changed:'+path)
 for f in cfg['training']+cfg['validation']+[cfg['common_preparation']]:
  if hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()!=f['sha256']:raise ValueError('Fixture changed')
 return cfg

def validator(params,seed,default=False):
 reference=OLD/'instrumentation_preparation/w2lvar_case.release_state.npz'
 def check(path):
  checks={}
  with np.load(path) as a,np.load(reference) as b:
   for key in b.files:
    if key=='transformed_reserve_parameters':
     expected=b[key].copy();expected[0]=float(params['pressure_cap_kpa']>0);expected[1]=.18;expected[5]=1.;checks[key]=bool(np.array_equal(a[key],expected))
    elif key=='transformed_valve_pressure_cap_kpa':checks[key]=float(a[key])==params['pressure_cap_kpa']
    elif key=='transformed_capacity_parameters':checks[key]=bool(np.array_equal(a[key],[.75,3,.02,.05,seed]))
    elif key=='transformed_capacity_noise_table_hash':
     rng=np.random.default_rng(seed);noise=np.clip(rng.standard_normal((30000,2,2)),-3,3);noise[:,:,0]*=.02;noise[:,:,1]*=.05;checks[key]=str(a[key])==hashlib.sha256(noise.tobytes()).hexdigest()
    else:checks[key]=bool(np.array_equal(a[key],b[key]))
   extras=set(a.files)-set(b.files);expected={'transformed_authority_control_mode','transformed_authority_replay_hash','transformed_authority_replay_count'}
   checks['only_declared_extra_fields']=extras==expected;checks['replay_counter_zero']=int(a['transformed_authority_replay_count'])==0;checks['control_mode_declared']=str(a['transformed_authority_control_mode'])==CONTROL;checks['replay_hash_declared']=str(a['transformed_authority_replay_hash'])==('none' if REPLAY is None else hashlib.sha256(REPLAY.tobytes()).hexdigest())
  if not all(checks.values()):raise ValueError('Common state mismatch:'+str([k for k,v in checks.items() if not v]))
  return checks
 return check

def summarize(c,folder,duration):
 r=ORIGINAL_SUMMARIZE(c,folder,duration);physical=json.loads(Path(c['result']).read_text());r['actual_pressure_peak_kpa']=physical.get('pressure_peak_kpa');r['gates']['physical_pressure_ceiling']=r['actual_pressure_peak_kpa'] is not None and r['actual_pressure_peak_kpa']<=620.+1e-9;r['gates']['native_force_parity']=physical.get('native_pressure_force_max_error_N',1.)==0.;r['accepted']=all(r['gates'].values())
 if not r['accepted'] and r['objective']<10000:r['objective']+=10000.
 r['control_mode']=CONTROL;(folder/'CASE_SUMMARY.json').write_text(json.dumps(r,indent=2),encoding='utf-8');return r

def setup():
 base.DEST=DEST;base.registration=registration;base.validator=validator;base.CapacityReserveForward=AuthorityForward;base.summarize=summarize

def evaluate_case(params,spec,folder,duration=20.,full=False):
 setup();return base.evaluate_case(params,spec,folder,duration=duration,full=full)

def preflight():
 cfg=registration();gate=DEST/'BASELINE_GATE.json'
 if gate.exists():return json.loads(gate.read_text())
 # The same cap5 source is compared to its previously retained full-array hashes.
 p=dict(estimator_mode='pressure',capacity_gain=.75,pressure_cap_kpa=5.,preview_s=.18);r=evaluate_case(p,cfg['training'][0],DEST/'cap5_parity',full=True);reference=OLD/'pressure/training_cases/9cac6942f10c6400/train_nominal/TELEMETRY_MANIFEST.json';checks={}
 for item in json.loads(reference.read_text())['archives']:
  name=Path(item['original_path']).name
  if name in ('ESTIMATOR_SENSORS.npz',) or name.endswith(('.final_neural_state.npz','.release_state.npz')):continue
  path=DEST/'cap5_parity'/name
  with np.load(path) as z:
   for key,entry in item['arrays'].items():checks[name+':'+key]=hashlib.sha256(z[key].tobytes()).hexdigest()==entry['sha256']
 if not all(checks.values()):raise ValueError('Cap5 inherited path parity failed')
 result=dict(ok=True,ready_for_optimization=r['accepted'],cap5_old_physical_exact=checks,baseline=r,scope='New 6..10kPa neural reference authority only. Physical pressure/rate/geometry/force law unchanged; estimator gamma0.75/preview180ms frozen from prior training.');gate.write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def handle_request(request):
 cfg=registration();m=request.get('manifest',{})
 if m.get('registration_sha256')!=hashlib.sha256((DEST/'PREREGISTRATION.json').read_bytes()).hexdigest():raise ValueError('Registration changed')
 if request['op']=='preflight':return preflight()
 p=request['params']
 if request['op']!='evaluate' or set(p)!={'estimator_mode','capacity_gain','pressure_cap_kpa','preview_s'} or p['capacity_gain']!=.75 or p['preview_s']!=.18 or p['estimator_mode']!='pressure' or not 6<=p['pressure_cap_kpa']<=10:raise ValueError('Unregistered or unused fit knobs')
 if not json.loads((DEST/'BASELINE_GATE.json').read_text())['ready_for_optimization']:raise ValueError('Parity gate absent')
 key=hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest()[:16];folder=DEST/'pressure/training_cases'/key;folder.mkdir(parents=True,exist_ok=True);rows=[dict(condition=s['name'],**evaluate_case(p,s,folder/s['name'])) for s in cfg['training']];failed=sum(not r['accepted'] for r in rows);value=10000.*failed+sum(max(0.,19.4-r['moving_s'])*100 for r in rows if not r['accepted'])+max((r['objective'] for r in rows if r['accepted']),default=0.)+sum(.1*r['mean_added_reference_kpa']+.2*r['saturation_fraction'] for r in rows)
 result=dict(ok=True,value=value,diagnostics=dict(params=p,all_training_accepted=failed==0,failed_conditions=failed,conditions=rows,case_cache=str(folder),evaluation_seed=request['evaluation_seed']));(folder/'TRIAL_SUMMARY.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def child(name):
 folder=DEST/'pressure'/name;folder.mkdir(parents=True,exist_ok=True);p=folder/'CASE_SUMMARY.json'
 if not p.exists():
  with (folder/'RUN.log').open('w',encoding='utf-8') as log:process=subprocess.run([sys.executable,str(Path(__file__)),'--case',name],stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
  if process.returncode:raise RuntimeError(str(folder/'RUN.log'))
 return json.loads(p.read_text())

def validate():
 import capacity_comparison_validation as v
 v.DEST=DEST;v.registration=registration;v.evaluate_case=evaluate_case;v.child=lambda mode,name:child(name);v.run('pressure')

def run_one(name):
 global CONTROL,REPLAY
 cfg=registration();selection=json.loads((DEST/'pressure/SELECTED.json').read_text());fresh=name.startswith('fresh_nominal_');matched=name.startswith(('replay_','zero_'));duration=30. if fresh else 20.;spec=cfg['training'][0] if fresh else next(s for s in cfg['training']+cfg['validation'] if s['name']==name.removeprefix('replay_').removeprefix('zero_'))
 p=selection['params'].copy()
 if name.startswith('zero_'):p['pressure_cap_kpa']=0.;CONTROL='zero'
 elif name.startswith('replay_'):
  CONTROL='replay'
  with np.load(DEST/'pressure/fresh_nominal_a/w2lvar_case.relay.npz') as z:REPLAY=z['samples'][:,:,3].copy()
 print(json.dumps(evaluate_case(p,spec,DEST/'pressure'/name,duration=duration,full=fresh),indent=2))

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--validate',action='store_true');a.add_argument('--case');args=a.parse_args()
 if args.case:run_one(args.case)
 elif args.validate:validate()
 else:print(json.dumps(preflight(),indent=2))
