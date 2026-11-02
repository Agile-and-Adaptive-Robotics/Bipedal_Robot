"""Opt-in predicted geometric/actual-tension error → positive MN recruitment."""
import os
for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
os.environ['CONDA_PREFIX']=r'D:\Anaconda\envs\myo'
import json,sys,hashlib,argparse,subprocess,time
from pathlib import Path
import numpy as np
import capacity_comparison_adapter as B
import valve_authority_experiment as V
from pressure_reserve_neural import STATE_KEYS,SPARSE_KEYS
from register_capacity_comparison import ROOT
DEST=ROOT/'source_domain_recruitment';HERE=Path(__file__).resolve().parent
GAIN=0.

class Recruitment:
 def __init__(self,inner,names):
  self.inner=inner;self.outer=None;self.cap=float(GAIN);self.relay=np.zeros(2);self.steps=0;self.logs=[];self.buffer=np.zeros(inner.c.input_connectivity.shape[1]);self.names=list(names);self.mn=np.array([self.names.index('MN_pect_'+s) for s in ('r','l')]);self.ports=[]
  if not 0<=self.cap<=8.:raise ValueError('Declared positive MN current cap0..8nA')
  for i in self.mn:
   p=np.flatnonzero(inner.c.input_connectivity[i])
   if len(p)!=1 or inner.c.input_connectivity[i,p[0]]!=1 or np.count_nonzero(inner.c.input_connectivity[:,p[0]])!=1:raise ValueError('Dedicated MN input transport required')
   self.ports.append(int(p[0]))
  self.ports=np.array(self.ports)
 def __getattr__(self,key):return getattr(self.inner,key)
 def forward(self,x):
  np.copyto(self.buffer,x);row=np.zeros((2,8));o=self.outer
  if self.cap>0 and o.context.get('active',False) and o.initialized:
   a=o.context['adapter'];idx=np.array([a.bank.names.index(n) for n in o.routes]);p=a.bank.parameters[idx];observed=o.buffer[(o.buffer_count-1-o.delay_steps)%len(o.buffer)];length=observed[:,2];velocity=observed[:,3];tension=o.previous_force.copy();minimum=p[:,4]*np.maximum(0.,length-p[:,0]-p[:,2]);predicted=p[:,4]*np.maximum(0.,length+(.18+o.delay_steps*.002)*np.maximum(velocity,0.)-p[:,0]-p[:,2]);signal=np.where(predicted>0,np.clip((predicted+5.-tension)/5.,0.,1.),0.);old=self.relay.copy();current=np.clip((self.cap/40.)*np.clip(old/10.,0,1)*(40.-self.inner.c.V[self.mn]),0.,self.cap);self.buffer[self.ports]+=current;self.relay+=.002/.006*(10.*signal-self.relay);target=np.clip(self.inner.c.V[self.mn]/B.runner.F.R.E_HI,0,1)*a.bank.scales[idx];row[:]=np.column_stack((minimum,predicted,tension,target,signal,old,current,np.full(2,self.cap)))
  result=self.inner.forward(x if self.cap==0 else self.buffer);self.logs.append(row);self.steps+=1;return result
 def export_state(self):
  d=dict(recruit_parameters=np.array([self.cap,.002,.006,.18,5.,5.,40.]),recruit_ports=self.ports.copy(),recruit_names=np.array(['pect_r','pect_l']),recruit_relay=self.relay.copy(),recruit_input_buffer=self.buffer.copy(),recruit_steps=np.array(self.steps),recruit_logs=np.asarray(self.logs).reshape(-1,2,8),recruit_rng_state=np.empty(0,dtype=np.uint8))
  for k in STATE_KEYS:
   if hasattr(self.inner.c,k):d['backend_'+k]=np.asarray(getattr(self.inner.c,k)).copy()
  for k in SPARSE_KEYS:
   if hasattr(self.inner,k):d['sparse_'+k]=np.asarray(getattr(self.inner,k)).copy()
  return d
 def import_state(self,d):
  assert np.array_equal(d['recruit_parameters'],[self.cap,.002,.006,.18,5.,5.,40.]) and len(d['recruit_rng_state'])==0
  self.relay[:]=d['recruit_relay'];self.buffer[:]=d['recruit_input_buffer'];self.steps=int(d['recruit_steps']);self.logs=[r.copy() for r in d['recruit_logs']]
  for k in STATE_KEYS:
   if 'backend_'+k in d:setattr(self.inner.c,k,d['backend_'+k].copy())
  for k in SPARSE_KEYS:
   if 'sparse_'+k in d:np.copyto(getattr(self.inner,k),d['sparse_'+k])

class RecruitmentOuter(V.AuthorityForward):
 def __init__(self,inner,names,*args,**kwargs):
  recruit=Recruitment(inner,names);super().__init__(recruit,names,*args,**kwargs);self.recruit=recruit;recruit.outer=self
 def export_state(self):
  d=super().export_state()
  for key in list(d):
   if key.startswith(('inner_backend_','inner_sparse_','inner_recruit_')):d[key[6:]]=d.pop(key)
  return d
 def import_state(self,state):
  d={('inner_'+k if k.startswith(('backend_','sparse_','recruit_')) else k):v for k,v in state.items()};super().import_state(d)
 def forward(self,x=None):
  if self.context.get('active',False) and not self.restore_verified and max(float(np.max(self.relay_v)),float(np.max(self.recruit.relay)))>1e-8:
   initial=self.export_state();a=self._step(x).copy();end=self.export_state();self.import_state(initial);b=self._step(x);repeat=self.export_state();assert np.array_equal(a,b)
   for k in end:assert np.array_equal(end[k],repeat[k]),'Recruitment complete state mismatch:'+k
   self.restore_verified=True;return b
  return self._step(x)
 def circuit(self):
  d=super().circuit();d.update(recruitment_architecture='Two mirrored6ms leaky error relays → positive capped homonymous MN conductance currents; valve10 baseline remains live',mn_current_cap_nA=self.recruit.cap,error='clip((K*max(0,L+(180ms+transport)*positivevelocity-rest-fixed)+5N-Tmeasured)/5N,0,1), gated by predicted positiveextension',force_or_voltage_projected=False,original_synapses_cancelled=False,oracle_gain_used=False,eligible_routes=['pect_r','pect_l']);return d

def registration():
 cfg=json.loads((DEST/'PREREGISTRATION.json').read_text())
 for p,e in json.loads((DEST/'source_snapshot/MANIFEST.json').read_text()).items():
  if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=e['sha256']:raise ValueError('Frozen source changed:'+p)
 for f in cfg['training']+cfg['validation']+[cfg['common_preparation']]:
  if hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()!=f['sha256']:raise ValueError('Fixture changed')
 return cfg

def validator(params,seed,default=False):
 reference=DEST/'PREPARATION_REFERENCE.npz'
 initial=not reference.exists()
 if initial:reference=ROOT/'valve_authority/cap5_parity/w2lvar_case.release_state.npz'
 def check(path):
  checks={}
  with np.load(path) as a,np.load(reference) as b:
   for key in b.files:
    if key=='transformed_valve_pressure_cap_kpa':checks[key]=float(a[key])==10.
    elif key=='transformed_recruit_parameters':checks[key]=bool(np.array_equal(a[key],[GAIN,.002,.006,.18,5.,5.,40.]))
    elif key=='transformed_capacity_parameters':checks[key]=bool(np.array_equal(a[key],[.75,3,.02,.05,seed]))
    elif key=='transformed_capacity_noise_table_hash':
     rng=np.random.default_rng(seed);noise=np.clip(rng.standard_normal((30000,2,2)),-3,3);noise[:,:,0]*=.02;noise[:,:,1]*=.05;checks[key]=str(a[key])==hashlib.sha256(noise.tobytes()).hexdigest()
    else:checks[key]=bool(np.array_equal(a[key],b[key]))
   if initial:
    extras=set(a.files)-set(b.files);checks['only_explicit_recruit_states']=bool(extras) and all(k.startswith('transformed_recruit_') for k in extras);checks['new_recruit_resting']=bool(np.all(a['transformed_recruit_relay']==0) and np.all(a['transformed_recruit_logs']==0))
   else:checks['same_complete_fields']=set(a.files)==set(b.files)
  if not all(checks.values()):raise ValueError('Recruitment common preparation mismatch:'+str([k for k,v in checks.items() if not v]))
  return checks
 return check

def summarize(c,folder,duration):
 r=V.ORIGINAL_SUMMARIZE(c,folder,duration);physical=json.loads(Path(c['result']).read_text());prefix=Path(c['result'].removesuffix('.result.json'))
 with np.load(prefix.with_suffix('.apparatus.npz')) as z:rig=z['samples'][:,10]/physical['body_weight_N']
 r['lateral_rig_rms_bw']=float(np.sqrt(np.mean(rig**2)));r['lateral_rig_peak_bw']=float(np.max(np.abs(rig)));r['actual_pressure_peak_kpa']=physical.get('pressure_peak_kpa');r['gates'].update(lateral_rig_rms=r['lateral_rig_rms_bw']<=.05,lateral_rig_peak=r['lateral_rig_peak_bw']<=.15,pressure_ceiling=r['actual_pressure_peak_kpa']<=620.+1e-9,native_force_parity=physical.get('native_pressure_force_max_error_N')==0.)
 with np.load(prefix.with_suffix('.final_neural_state.npz')) as z:logs=z['recruit_logs'];currents=logs[:,:,6];r['mean_added_MN_current_nA']=float(np.mean(currents));r['max_added_MN_current_nA']=float(np.max(currents));r['gates']['positive_capped_MN_current']=bool(np.min(currents)>=0 and np.max(currents)<=GAIN+1e-12)
 r['accepted']=all(r['gates'].values());r['objective']+=.1*r['mean_added_MN_current_nA']+100*r['lateral_rig_rms_bw']
 if not r['accepted'] and r['objective']<10000:r['objective']+=10000.
 if not (DEST/'PREPARATION_REFERENCE.npz').exists():
  with np.load(prefix.with_suffix('.release_state.npz')) as z:np.savez_compressed(DEST/'PREPARATION_REFERENCE.npz',**{k:z[k] for k in z.files})
 (folder/'CASE_SUMMARY.json').write_text(json.dumps(r,indent=2),encoding='utf-8');return r

def evaluate_case(params,spec,folder,duration=20.,full=False):
 global GAIN
 GAIN=params['mn_current_cap_nA'];B.DEST=DEST;B.registration=registration;B.validator=validator;B.CapacityReserveForward=RecruitmentOuter;B.summarize=summarize;return B.evaluate_case(params,spec,folder,duration=duration,full=full)

def preflight():
 path=DEST/'BASELINE_GATE.json'
 if path.exists():return json.loads(path.read_text())
 cfg=registration();p=cfg['fixed_parameters']|dict(mn_current_cap_nA=0.);r=evaluate_case(p,cfg['training'][0],DEST/'zero_default_parity');old=json.loads((ROOT/'valve_authority/pressure/training_cases/41a77bc5159f9d6c/train_nominal/TELEMETRY_MANIFEST.json').read_text());new=json.loads((DEST/'zero_default_parity/TELEMETRY_MANIFEST.json').read_text());lookup={Path(a['original_path']).name:a for a in new['archives']};checks={}
 for a in old['archives']:
  name=Path(a['original_path']).name
  if name.endswith(('.release_state.npz','.final_neural_state.npz')):continue
  for key,e in a['arrays'].items():checks[name+':'+key]=lookup[name]['arrays'][key]['sha256']==e['sha256']
 if not all(checks.values()):raise ValueError('Zero recruitment physical parity failed')
 result=dict(ok=True,ready_for_optimization=r['accepted'],baseline=r,all_physical_hashes_exact=checks,scope='Only new MN pathway off; inherited valve10 remains live and byteexact. Added lateral gates must pass before fitting.');path.write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def handle_request(request):
 cfg=registration()
 if request.get('manifest',{}).get('registration_sha256')!=hashlib.sha256((DEST/'PREREGISTRATION.json').read_bytes()).hexdigest():raise ValueError('Registration changed')
 if request['op']=='preflight':return preflight()
 p=request['params']
 if request['op']!='evaluate' or set(p)!=set(cfg['fixed_parameters'])|{'mn_current_cap_nA'} or any(p[k]!=v for k,v in cfg['fixed_parameters'].items()) or not 1<=p['mn_current_cap_nA']<=8.:raise ValueError('Unregistered fit knobs')
 if not json.loads((DEST/'BASELINE_GATE.json').read_text())['ready_for_optimization']:raise ValueError('Parity/lateral baseline gate failed')
 key=hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest()[:16];folder=DEST/'pressure/training_cases'/key;folder.mkdir(parents=True,exist_ok=True);rows=[dict(condition=s['name'],**evaluate_case(p,s,folder/s['name'])) for s in cfg['training']];failed=sum(not r['accepted'] for r in rows);value=10000.*failed+sum(max(0.,19.4-r['moving_s'])*100 for r in rows if not r['accepted'])+max((r['objective'] for r in rows if r['accepted']),default=0.)+sum(.1*r['mean_added_reference_kpa']+.1*r['mean_added_MN_current_nA'] for r in rows);result=dict(ok=True,value=value,diagnostics=dict(params=p,all_training_accepted=failed==0,failed_conditions=failed,conditions=rows,case_cache=str(folder),evaluation_seed=request['evaluation_seed']));(folder/'TRIAL_SUMMARY.json').write_text(json.dumps(result,indent=2),encoding='utf-8');return result

def child(name):
 folder=DEST/'pressure'/name;folder.mkdir(parents=True,exist_ok=True);p=folder/'CASE_SUMMARY.json'
 if not p.exists():
  with (folder/'RUN.log').open('w',encoding='utf-8') as log:proc=subprocess.run([sys.executable,str(Path(__file__)),'--case',name],stdout=log,stderr=subprocess.STDOUT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
  if proc.returncode:raise RuntimeError(str(folder/'RUN.log'))
 return json.loads(p.read_text())

def validate():
 import capacity_comparison_validation as C
 C.DEST=DEST;C.registration=registration;C.evaluate_case=evaluate_case;C.child=lambda mode,name:child(name);C.run('pressure')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--validate',action='store_true');p.add_argument('--case');a=p.parse_args()
 if a.validate:validate()
 elif a.case:
  cfg=registration();selection=json.loads((DEST/'pressure/SELECTED.json').read_text());fresh=a.case.startswith('fresh_nominal_');off=a.case.startswith('off_');spec=cfg['training'][0] if fresh else next(s for s in cfg['training']+cfg['validation'] if s['name']==a.case.removeprefix('off_'));params=selection['params'].copy()
  if off:params['mn_current_cap_nA']=0.
  print(json.dumps(evaluate_case(params,spec,DEST/'pressure'/a.case,duration=30. if fresh else 20.,full=fresh),indent=2))
 else:print(json.dumps(preflight(),indent=2))
