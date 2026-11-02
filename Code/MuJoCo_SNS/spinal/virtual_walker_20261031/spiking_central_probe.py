"""Neural-only Ca-burster threshold/reset feasibility gate, before mechanics.

Twelve RG/PF Ca cells become populations of actual spiking Ca/LIF cells;
other source cells remain graded. This is not the passing sensory hybrid.
"""
import hashlib,json,math,time
from pathlib import Path
import numpy as np
from numba import njit
import native_kernel as K
from import_virtual_walker import HERE
from recruitment import apply_recruitment

@njit(cache=True)
def integrate(v,m,h,params,edges,tonic,stimuli,encoded,pf_mask,settings,duration,dt):
 steps=round(duration/dt);stride=round(.002/dt);P=len(encoded);N=16  # Explicit finite source-cell populations.
 lv=np.empty((P,N));lm=np.empty((P,N));lh=np.empty((P,N));ref=np.zeros((P,N));adapt=np.zeros((P,N));decoded=np.empty(P);total_spikes=np.zeros(P,dtype=np.int64)
 for p in range(P):
  i=encoded[p];decoded[p]=max(0.,(v[i]-params[i,0])/20.)
  for j in range(N):lv[p,j]=v[i]+.2*(j/(N-1)-.5);lm[p,j]=m[i];lh[p,j]=h[i]
 gs=np.zeros(len(v));gi=gs.copy();applied=gs.copy();newv=v.copy();trace=np.empty((steps//stride,len(v)));events=np.zeros((steps//stride,P),dtype=np.int64);counts=np.zeros(P,dtype=np.int64);soma_trace=np.empty((steps//stride,P,2));gate_trace=np.empty((steps//stride,P));adapt_trace=np.empty((steps//stride,P))
 gscale,threshold,reset,syn_tau,rate_max=settings[:5];jump=settings[5] if len(settings)>5 else 0.;tau_adapt=settings[6] if len(settings)>6 else .35;pf_ca=settings[7] if len(settings)>7 else 1.;pf_jump=settings[8] if len(settings)>8 else 1.;pf_tau=settings[9] if len(settings)>9 else 1.
 for tick in range(steps):
  t=tick*dt;applied[:]=0.
  for s in stimuli:
   if s[1]<=t<s[2]:applied[int(s[0])]+=s[3]
  K.neural_steps(v,m,h,newv,gs,gi,params,edges,tonic,applied,dt*1000.,1)
  for p in range(P):
   i=encoded[p];rest,tau,size,gmax,vm,sm,tm,vh,sh,th,eca=params[i];spikes=0
   for j in range(N):
    oldv=lv[p,j];gca=gmax*gscale*(pf_ca if pf_mask[p] else 1.)*lm[p,j]*lh[p,j];adapt[p,j]*=math.exp(-dt/(tau_adapt*(pf_tau if pf_mask[p] else 1.)))
    z=math.exp(-sm*(oldv-vm));inf=1/(1+z);tc=inf*math.sqrt(z)*tm;lm[p,j]+=(inf-lm[p,j])*(-math.expm1(-dt*1000./tc))
    z=.5*math.exp(-sh*(oldv-vh));inf=1/(1+z);tc=inf*math.sqrt(z)*th;lh[p,j]+=(inf-lh[p,j])*(-math.expm1(-dt*1000./tc))
    if ref[p,j]>0:ref[p,j]=max(0.,ref[p,j]-dt);lv[p,j]=reset;continue
    conductance=1.+gs[i]+gca;decay=math.exp(-conductance*dt*1000./tau);current=(tonic[i]+applied[i])/size+gi[i]+gca*(eca-rest)-adapt[p,j]
    lv[p,j]=rest+(oldv-rest)*decay+current*(1-decay)/conductance
    if lv[p,j]>=threshold:lv[p,j]=reset;ref[p,j]=.001;adapt[p,j]+=jump*(pf_jump if pf_mask[p] else 1.);spikes+=1
   decoded[p]=decoded[p]*math.exp(-dt/syn_tau)+spikes/(N*rate_max*syn_tau);v[i]=rest+20.*decoded[p];counts[p]+=spikes;total_spikes[p]+=spikes
  if tick%stride==stride-1:
   row=tick//stride;trace[row]=v;events[row]=counts;counts[:]=0
   for p in range(P):soma_trace[row,p]=lv[p,:2];gate_trace[row,p]=np.mean(lh[p]);adapt_trace[row,p]=np.mean(adapt[p])
 return trace,events,total_spikes,soma_trace,gate_trace,adapt_trace

STATE={}
def prepare(manifest):
 cfg=json.loads((HERE/'GENTLE_HIP_0P6_FLEX_0P5.json').read_text());source=Path(cfg['source']);assert hashlib.sha256(source.read_bytes()).hexdigest()==cfg['expected_source_sha256']
 bank=json.loads(Path(manifest['bank_audit']).read_text())['muscles'];c=K.NativeController(source,bank,ankle_pf_crossed=True);apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);assert len(c.v)==108 and len(c.edges)==216
 encoded=np.array([i for i,n in enumerate(c.names) if c.params[i,3]>0 and (' RG ' in n or ' PF ' in n)],dtype=np.int64);assert len(encoded)==12;assert c.dt==.0002
 STATE.update(controller=c,encoded=encoded,pf_mask=np.array([' PF ' in c.names[i] for i in encoded]),duration=10.,output=Path(manifest['output_dir']))

def metrics(trace):
 c=STATE['controller'];t=np.arange(1,len(trace)+1)*.002;out={};onsets={}
 for i in STATE['encoded']:
  name=c.names[i];v=np.clip(trace[:,i],-60.,-40.);tail=v[t>=2.];lo=float(np.min(tail));hi=float(np.max(tail));rises=[];active=False;last=-1e6
  if hi-lo>=3.:
   lower=lo+.3*(hi-lo);upper=lo+.7*(hi-lo)
   for k in np.flatnonzero(t>=2.):
    if active and v[k]<lower:active=False
    elif not active and v[k]>=upper:
     active=True
     if t[k]-last>=.25:rises.append(float(t[k]));last=t[k]
  period=float(np.median(np.diff(rises))) if len(rises)>=3 else None;onsets[name]=rises;out[name]=dict(cycles=len(rises),period_s=period,range_mV=hi-lo)
 L=onsets['L RG flx'];R=np.array(onsets['R RG flx']);errors=[]
 for a,b in zip(L[:-1],L[1:]):
  candidates=R[(R>=a)&(R<b)]
  if len(candidates):errors.append(abs(float((candidates[0]-a)/(b-a))-.5))
 phase=float(np.mean(errors)) if errors else .5
 phases={}
 for name,rises in onsets.items():
  samples=[(r-a)/(b-a) for a,b in zip(L[:-1],L[1:]) for r in rises if a<=r<b]
  phases[name]=float(np.angle(np.mean(np.exp(2j*np.pi*np.array(samples))))/(2*np.pi)%1.) if samples else None
 return dict(cells=out,RG_phase_error=phase,phases=phases)

def motor_profiles(trace):
 c=STATE['controller'];t=np.arange(1,len(trace)+1)*.002;v=np.clip(trace[:,c.names.index('L RG flx')],-60.,-40.);tail=v[t>=2.];lo=float(np.min(tail));hi=float(np.max(tail));active=False;rises=[]
 for k in np.flatnonzero(t>=2.):
  if active and v[k]<lo+.3*(hi-lo):active=False
  elif not active and v[k]>=lo+.7*(hi-lo):
   active=True
   if not rises or t[k]-rises[-1]>=.25:rises.append(float(t[k]))
 ids=[i for i,n in enumerate(c.names) if ' MN ' in n and not n.endswith(' RE')];profiles=[]
 for a,b in zip(rises[:-1],rises[1:]):profiles.append(np.stack([np.interp(np.linspace(a,b,100,endpoint=False),t,trace[:,i]) for i in ids]))
 return np.mean(profiles,axis=0) if profiles else np.full((len(ids),100),np.nan),[c.names[i] for i in ids]

def handle_request(request):
 if request['op']=='preflight':
  prepare(request['manifest']);c=STATE['controller'];v=c.v.copy();m=c.m.copy();h=c.h.copy();nv=v.copy();gs=np.zeros(len(v));gi=gs.copy();applied=gs.copy();trace=[]
  for tick in range(round(STATE['duration']/c.dt)):
   t=tick*c.dt;applied[:]=0.
   for i,a,b,current in c.stimuli:
    if a<=t<b:applied[i]+=current
   K.neural_steps(v,m,h,nv,gs,gi,c.params,c.edges,c.tonic,applied,c.dt*1000.,1)
   if tick%10==9:trace.append(v.copy())
  baseline=metrics(np.array(trace));STATE['baseline']=baseline;STATE['target_profiles'],STATE['motor_names']=motor_profiles(np.array(trace));np.savez_compressed(STATE['output']/'GRADED_NEURAL_BASELINE.npz',t=np.arange(1,len(trace)+1)*.002,v=trace,names=c.names)
  settings=np.array([1.,-45.,-60.,.02,200.]);integrate(c.v.copy(),c.m.copy(),c.h.copy(),c.params,c.edges,c.tonic,np.array(c.stimuli),STATE['encoded'],STATE['pf_mask'],settings,.02,c.dt)
  ready=all(x['cycles']>=5 and x['period_s'] is not None for x in baseline['cells'].values()) and baseline['RG_phase_error']<.1
  return dict(ok=True,ready_for_optimization=ready,baseline=baseline,scope='Neural-only source rhythm matching; not biomechanical training or walking.',source_sha256=hashlib.sha256(Path(c_source(request)).read_bytes()).hexdigest())
 if request['op']=='evaluate':
  c=STATE['controller'];search=request['params'];output_fit=set(search)=={'pf_output_scale'};phase_keys=('pf_ca_scale','pf_adapt_jump_scale','pf_adapt_tau_scale');phase_fit=set(search)==set(phase_keys);p=request['manifest']['central_parameters'] if output_fit or phase_fit else search;keys=('gCa_scale','threshold_mV','reset_mV','syn_tau_s','readout_rate_max')+(('adapt_jump_nA','adapt_tau_s') if 'adapt_jump_nA' in p else ());assert set(p)==set(keys);settings=np.array([p[k] for k in keys]+([search[k] for k in phase_keys] if phase_fit else []));tick=time.perf_counter();edges=c.edges.copy()
  if output_fit or phase_fit:
   output_scale=search['pf_output_scale'] if output_fit else request['manifest']['pf_output_scale']
   for k,e in enumerate(edges):
    pre,post=c.names[int(e[0])],c.names[int(e[1])]
    if ' PF ' in pre and not pre.endswith(' IN') and ' MN ' in post and not post.endswith(' RE'):edges[k,2]*=output_scale
  if p['reset_mV']>=p['threshold_mV']-3.:return dict(ok=True,value=1000.,diagnostics={'reason':'reset must be at least3mV below threshold'})
  v,events,counts,somas,gates,adaptation=integrate(c.v.copy(),c.m.copy(),c.h.copy(),c.params,edges,c.tonic,np.array(c.stimuli),STATE['encoded'],STATE['pf_mask'],settings,STATE['duration'],c.dt);finite=bool(np.all(np.isfinite(v)) and np.all(np.isfinite(somas)) and np.all(np.isfinite(gates)));q=metrics(v) if finite else {};value=500.
  if finite:
   value=10*q['RG_phase_error']
   for name,a in q['cells'].items():
    target=STATE['baseline']['cells'][name]['period_s'];value+=10.*max(0,5-a['cycles'])+(10. if a['period_s'] is None else abs(math.log(a['period_s']/target)))
  phase_errors={n:abs((p0-STATE['baseline']['phases'][n]+.5)%1.-.5) if p0 is not None and STATE['baseline']['phases'][n] is not None else .5 for n,p0 in q.get('phases',{}).items()}
  gate=finite and all(a['cycles']>=5 and a['period_s'] is not None and abs(a['period_s']/STATE['baseline']['cells'][n]['period_s']-1.)<=.2 for n,a in q.get('cells',{}).items()) and q.get('RG_phase_error',.5)<.1 and bool(np.all(counts>0)) and bool(phase_errors) and max(phase_errors.values())<=.1
  profiles,motor_names=motor_profiles(v);profile_error=float(np.mean(((profiles-STATE['target_profiles'])/np.maximum(1.,np.ptp(STATE['target_profiles'],axis=1,keepdims=True)))**2))
  if output_fit:value=profile_error+(0. if gate else 500.)
  if phase_fit:value+=5.*sum(phase_errors.values())+profile_error
  if value<STATE.get('best',float('inf')):
   STATE['best']=value;np.savez_compressed(STATE['output']/'BEST_NEURAL_TRACE.npz',t=np.arange(1,len(v)+1)*.002,v=v,spike_counts=events,names=c.names,encoded_indices=STATE['encoded'],parameters_json=json.dumps(p),searched_parameters_json=json.dumps(search),soma_membrane_mV=somas,h_gate_mean=gates,adaptation_current_nA=adaptation,motor_profiles=profiles,motor_target_profiles=STATE['target_profiles'],motor_names=motor_names)
  return dict(ok=True,value=float(value),diagnostics=dict(neural_gate=bool(gate),finite=finite,rhythm=q,phase_errors=phase_errors,motor_profile_normalized_MSE=profile_error,spikes_per_population=counts.tolist(),wall_s=time.perf_counter()-tick,scope='192 actual Ca/LIF spiking cells replace12RG/PF cells;96 remaining cells graded. Optional spike-triggered inhibitory adaptation explicit. DesiredMN is the graded reference, not measuredEMG. No mechanics.'))
 raise ValueError('Unknown operation')

def c_source(request):return json.loads((HERE/'GENTLE_HIP_0P6_FLEX_0P5.json').read_text())['source']

if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--manifest',required=True);parser.add_argument('--parameters',required=True);parser.add_argument('--out',required=True);args=parser.parse_args();manifest=json.loads(Path(args.manifest).read_text());manifest['output_dir']=args.out;Path(args.out).mkdir(parents=True,exist_ok=True);parameters=json.loads(Path(args.parameters).read_text());parameters=parameters.get('trial',parameters).get('params',parameters)
 preflight=handle_request(dict(op='preflight',manifest=manifest));result=handle_request(dict(op='evaluate',manifest=manifest,params=parameters));record=dict(preflight=preflight,result=result,parameters=parameters,adapter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Fresh neural-only validation; effective presynaptic activation is clipped to source synapse threshold/saturation before phase measurement.');(Path(args.out)/'VALIDATION_RESULT.json').write_text(json.dumps(record,indent=2));print(json.dumps(result,indent=2))
