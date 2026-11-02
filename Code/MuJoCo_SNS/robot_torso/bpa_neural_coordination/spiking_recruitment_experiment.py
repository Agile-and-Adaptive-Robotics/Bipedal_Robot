"""Only the two source-domain error relays become genuine LIF populations."""
import json,hashlib,argparse,sys,math
from pathlib import Path
import numpy as np
from numba import njit
import mn_recruitment_experiment as E
BASE=E.DEST;DEST=BASE.parent/'source_domain_recruitment_spiking'

@njit(cache=True)
def lif_step(target,v,ref,decoded,counts,drive,gain):
 dt=.0001;tau_m=.005;refractory=.001;rate_max=200.;tau_syn=-.002/math.log(2./3.)
 leak=math.exp(-dt/tau_m);decay=math.exp(-dt/tau_syn);counts[:]=0
 for p in range(2):
  for j in range(v.shape[1]):
   rate=rate_max*min(1.,max(0.,target[p]))*(.8+.4*j/(v.shape[1]-1));drive[p,j]=0.
   if rate>0:drive[p,j]=1./(-math.expm1(-(1./rate-refractory)/tau_m))
 for sub in range(20):
  for p in range(2):
   count=0
   for j in range(v.shape[1]):
    if ref[p,j]>0:ref[p,j]=max(0.,ref[p,j]-dt);continue
    v[p,j]=drive[p,j]+(v[p,j]-drive[p,j])*leak
    if v[p,j]>=1.:v[p,j]=0.;ref[p,j]=refractory;count+=1
   decoded[p]=decoded[p]*decay+gain*count/(v.shape[1]*rate_max*tau_syn)
   counts[p]+=count

class Population:
 def __init__(self,gain=1.):
  self.gain=gain;self.v=np.tile(np.linspace(.001,.999,128),(2,1));self.ref=np.zeros((2,128));self.drive=np.zeros((2,128));self.decoded=np.zeros(2);self.counts=np.zeros(2,dtype=np.int64);self.total=np.zeros(2,dtype=np.int64);self.steps=0;self.events=[]
 def step(self,signal):
  lif_step(np.asarray(signal),self.v,self.ref,self.decoded,self.counts,self.drive,self.gain);self.total+=self.counts;self.steps+=1;self.events.append(self.counts.copy());return self.decoded
 def export(self):
  return dict(parameters=np.array([self.gain,128,.0001,20,.005,.001,200.,-.002/math.log(2/3),.8,1.2]),v=self.v.copy(),ref=self.ref.copy(),drive=self.drive.copy(),decoded=self.decoded.copy(),counts=self.counts.copy(),total=self.total.copy(),steps=np.array(self.steps),events=np.asarray(self.events,dtype=np.int64).reshape(-1,2),rng_state=np.empty(0,dtype=np.uint8))
 def restore(self,d):
  assert np.array_equal(d['parameters'],self.export()['parameters']) and len(d['rng_state'])==0
  for k in ('v','ref','drive','decoded','counts','total'):np.copyto(getattr(self,k),d[k])
  self.steps=int(d['steps']);self.events=[r.copy() for r in d['events']]

class SpikingRecruitment(E.Recruitment):
 def __init__(self,*args):
  super().__init__(*args);self.pop=Population(json.loads((DEST/'OFFLINE_CHECKS.json').read_text())['decoder_gain'])
 def forward(self,x):
  np.copyto(self.buffer,x);row=np.zeros((2,8));o=self.outer
  if self.cap>0 and o.context.get('active',False) and o.initialized:
   a=o.context['adapter'];idx=np.array([a.bank.names.index(n) for n in o.routes]);p=a.bank.parameters[idx];obs=o.buffer[(o.buffer_count-1-o.delay_steps)%len(o.buffer)];length=obs[:,2];velocity=obs[:,3];tension=o.previous_force.copy();minimum=p[:,4]*np.maximum(0.,length-p[:,0]-p[:,2]);predicted=p[:,4]*np.maximum(0.,length+(.18+o.delay_steps*.002)*np.maximum(velocity,0.)-p[:,0]-p[:,2]);signal=np.where(predicted>0,np.clip((predicted+5.-tension)/5.,0.,1.),0.);old=self.relay.copy();current=np.clip((self.cap/40.)*np.clip(old/10.,0,1)*(40.-self.inner.c.V[self.mn]),0.,self.cap);self.buffer[self.ports]+=current
   self.relay[:]=10.*self.pop.step(signal);target=np.clip(self.inner.c.V[self.mn]/E.B.runner.F.R.E_HI,0,1)*a.bank.scales[idx];row[:]=np.column_stack((minimum,predicted,tension,target,signal,old,current,np.full(2,self.cap)))
  result=self.inner.forward(x if self.cap==0 else self.buffer);self.logs.append(row);self.steps+=1;return result
 def export_state(self):
  d=super().export_state();d.update({'recruit_lif_'+k:v for k,v in self.pop.export().items()});return d
 def import_state(self,d):
  super().import_state(d);self.pop.restore({k.removeprefix('recruit_lif_'):v for k,v in d.items() if k.startswith('recruit_lif_')})

class SpikingOuter(E.RecruitmentOuter):
 def __init__(self,inner,names,*args,**kw):
  recruit=SpikingRecruitment(inner,names);E.V.AuthorityForward.__init__(self,recruit,names,*args,**kw);self.recruit=recruit;recruit.outer=self
 def circuit(self):
  d=super().circuit();d.update(spiking_scope='Only two new error relays: 128 LIF cells per side. Core888 and valve2 stay graded.',spiking_cells=256,graded_cells=890,spike_parameters=self.recruit.pop.export()['parameters'].tolist(),spike_counts=self.recruit.pop.total.tolist(),extra_transport_latency_s=0.,same_previous_sample_output_latency_s=.002);return d

def offline():
 DEST.mkdir(parents=True,exist_ok=True);path=DEST/'OFFLINE_CHECKS.json'
 if path.exists():return json.loads(path.read_text())
 rows=[]
 for level in (.05,.1,.2,.4,.6,.8,1.):
  p=Population();ys=[]
  for t in range(800):ys.append(p.step([level,level]).copy())
  rows.append(dict(input=level,raw_mean=float(np.mean(ys[300:]))))
 gain=sum(r['input']*r['raw_mean'] for r in rows)/sum(r['raw_mean']**2 for r in rows)
 for r in rows:r['decoded_mean']=gain*r['raw_mean'];r['absolute_error']=abs(r['decoded_mean']-r['input'])
 pop=Population(gain);graded=np.zeros(2);samples=[]
 for k in range(3000):
  t=k*.002;signal=np.full(2, .5+.45*math.sin(2*math.pi*2*t) if t<3 else (.9 if int(t*10)%2 else .1));old=graded.copy();decoded=pop.decoded.copy();pop.step(signal);graded+=(signal-graded)/3.;samples.append([t,signal[0],old[0],decoded[0]])
 a=np.asarray(samples);rmse=float(np.sqrt(np.mean((a[250:,2]-a[250:,3])**2)))
 saved=pop.export();first=[]
 for k in range(64):first.append(pop.step([.37,.72]).copy())
 end=pop.export();pop.restore(saved);second=[pop.step([.37,.72]).copy() for k in range(64)];again=pop.export();exact=bool(np.array_equal(first,second) and all(np.array_equal(end[k],again[k]) for k in end))
 zero=Population(gain)
 for k in range(200):zero.step([0.,0.])
 z=zero.export();zero_ok=bool(not np.any(z['total']) and not np.any(z['decoded']))
 result=dict(decoder_gain=gain,calibration='Only continuous electrical input mapping; no physical-case fitting',static_rows=rows,dynamic_rmse_fraction=rmse,state_restore_all_arrays_exact=exact,zero_input_no_spikes_or_current=zero_ok,passed=max(r['absolute_error'] for r in rows)<.04 and rmse<.15 and exact and zero_ok,policy=dict(cells_per_side=128,dt_s=.0001,rate_max_hz=200,tau_m_s=.005,refractory_s=.001,tau_filter_s=-.002/math.log(2/3),rng='None; deterministic initial phases',parameters_mirrored=True))
 np.savez_compressed(DEST/'OFFLINE_SIGNALS.npz',samples=a,channels=['t','input','graded_previous_output','spiking_previous_output']);path.write_text(json.dumps(result,indent=2));return result

def run(name):
 checks=offline();assert checks['passed'];cfg=E.registration();spec=cfg['training'][2 if name.startswith('strong') else 0];params=json.loads((BASE/'pressure/SELECTED.json').read_text())['params'];E.GAIN=params['mn_current_cap_nA']
 def validator(p,seed,default=False):
  check=E.validator(p,seed,default)
  def wrapped(path):
   # The only added prepared state belongs to resting deterministic LIF cells.
   with np.load(path) as a,np.load(BASE/'PREPARATION_REFERENCE.npz') as b:
    result={k:bool(np.array_equal(a[k],b[k])) for k in b.files if k not in ('transformed_capacity_parameters','transformed_capacity_noise_table_hash','transformed_recruit_parameters')};result['transformed_recruit_parameters']=np.array_equal(a['transformed_recruit_parameters'],[5.,.002,.006,.18,5.,5.,40.])
    for k in ('transformed_capacity_parameters','transformed_capacity_noise_table_hash'):
     if k.endswith('parameters'):result[k]=np.array_equal(a[k],[.75,3,.02,.05,seed])
     else:
      rng=np.random.default_rng(seed);noise=np.clip(rng.standard_normal((30000,2,2)),-3,3);noise[:,:,0]*=.02;noise[:,:,1]*=.05;result[k]=str(a[k])==hashlib.sha256(noise.tobytes()).hexdigest()
    extra=set(a.files)-set(b.files);pop=Population(checks['decoder_gain']).export();result['only_lif_added']=all(k.startswith('transformed_recruit_lif_') for k in extra) and len(extra)==len(pop)
    result['complete_lif_resting_state']=all(np.array_equal(a['transformed_recruit_lif_'+k],v) for k,v in pop.items())
   if not all(result.values()):raise ValueError('Spiking preparation mismatch:'+str([k for k,v in result.items() if not v]))
   return result
  return wrapped
 E.B.DEST=DEST;E.B.registration=E.registration;E.B.validator=validator;E.B.CapacityReserveForward=SpikingOuter;E.B.summarize=E.summarize
 folder=DEST/name;result=E.B.evaluate_case(params,spec,folder,duration=20. if name.startswith('strong') else 8.,full=True);print(json.dumps(result,indent=2));return result

def preflight():
 checks=offline();assert checks['passed']
 class Raw:
  input_connectivity=np.eye(2);V=np.array([0.,0.])
 class Inner:
  def __init__(self):self.c=Raw()
  def forward(self,x):return x+1.
 before=E.GAIN;E.GAIN=0.;old=E.Recruitment(Inner(),['MN_pect_r','MN_pect_l']);new=SpikingRecruitment(Inner(),['MN_pect_r','MN_pect_l']);exact=True
 for k in range(32):
  x=np.array([math.sin(k),math.cos(k)]);exact=exact and np.array_equal(old.forward(x),new.forward(x))
 E.GAIN=before
 with np.load(DEST/'OFFLINE_SIGNALS.npz') as z:a=z['samples']
 errors=[]
 for voltage in (-5.,0.,2.5,5.):
  graded=np.clip(5*np.clip(a[:,2],0,1)*(40-voltage)/40,0,5);spikes=np.clip(5*np.clip(a[:,3],0,1)*(40-voltage)/40,0,5);errors.append(float(np.sqrt(np.mean((graded[250:]-spikes[250:])**2))))
 r=dict(offline_pass=checks['passed'],zero_gain_forward_exact=exact,zero_gain_no_events=not np.any(new.pop.total),dynamic_current_rmse_nA=errors,current_gate_nA=.75,passed=exact and max(errors)<.75,physical_fit=False,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Two error relays only; frozen5nA MN cap.128LIF each, deterministic frequency spread0.8..1.2 mirrored, population mean200Hz maximum. No additional transport beyond same2ms output hold.')
 (DEST/'PHYSICAL_PREFLIGHT.json').write_text(json.dumps(r,indent=2));return r

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--case');a=p.parse_args()
 if a.case:
  assert preflight()['passed'];run(a.case)
 else:print(json.dumps(preflight(),indent=2))
