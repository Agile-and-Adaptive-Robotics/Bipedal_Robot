"""Opt-in adaptive Ca/LIF RG/PF populations; other source cells stay graded."""
import math,json
from types import FunctionType,SimpleNamespace
import numpy as np
from numba import njit
import native_kernel as K
BaseController=K.NativeController

@njit(cache=True)
def central_tick(v,m,h,newv,gs,gi,params,edges,tonic,applied,encoded,pf_mask,settings,lv,lm,lh,ref,adapt,decoded,counts,dt):
 K.neural_steps(v,m,h,newv,gs,gi,params,edges,tonic,applied,dt*1000.,1)
 gscale,threshold,reset,syn_tau,rate_max,jump,tau_adapt,pf_ca,pf_jump,pf_tau=settings[:10];rg_reset=settings[10] if len(settings)>10 else reset;N=lv.shape[1]
 for p in range(len(encoded)):
  i=encoded[p];rest,tau,size,gmax,vm,sm,tm,vh,sh,th,eca=params[i];spikes=0;cell_reset=reset if pf_mask[p] else rg_reset
  for j in range(N):
   oldv=lv[p,j];gca=gmax*gscale*(pf_ca if pf_mask[p] else 1.)*lm[p,j]*lh[p,j];adapt[p,j]*=math.exp(-dt/(tau_adapt*(pf_tau if pf_mask[p] else 1.)))
   z=math.exp(-sm*(oldv-vm));inf=1/(1+z);tc=inf*math.sqrt(z)*tm;lm[p,j]+=(inf-lm[p,j])*(-math.expm1(-dt*1000./tc))
   z=.5*math.exp(-sh*(oldv-vh));inf=1/(1+z);tc=inf*math.sqrt(z)*th;lh[p,j]+=(inf-lh[p,j])*(-math.expm1(-dt*1000./tc))
   if ref[p,j]>0:ref[p,j]=max(0.,ref[p,j]-dt);lv[p,j]=cell_reset;continue
   conductance=1.+gs[i]+gca;decay=math.exp(-conductance*dt*1000./tau);current=(tonic[i]+applied[i])/size+gi[i]+gca*(eca-rest)-adapt[p,j]
   lv[p,j]=rest+(oldv-rest)*decay+current*(1-decay)/conductance
   if lv[p,j]>=threshold:lv[p,j]=cell_reset;ref[p,j]=.001;adapt[p,j]+=jump*(pf_jump if pf_mask[p] else 1.);spikes+=1
  decoded[p]=decoded[p]*math.exp(-dt/syn_tau)+spikes/(N*rate_max*syn_tau);v[i]=rest+20.*decoded[p];counts[p]+=spikes

class SpikingCentralController(BaseController):
 population_size=16;configuration=None;last_instance=None;reference_calculations=False;sensory_spiking=False;audit_state=False;balance_configuration=None;clock_scale=1.
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs);config=self.configuration
  assert self.clock_scale>0 and np.isfinite(self.clock_scale)
  assert self.clock_scale==1 or self.balance_configuration is None,'Cadence comparison isolates the frozen core; combined task-relay clock not yet validated'
  from hip_knee_interlock import OPTIONAL_KEYS
  keys=('gCa_scale','threshold_mV','reset_mV','syn_tau_s','readout_rate_max','adapt_jump_nA','adapt_tau_s','pf_ca_scale','pf_adapt_jump_scale','pf_adapt_tau_scale');assert set(config)-set(OPTIONAL_KEYS)-{'rg_reset_mV','ankle_pf_output_multiplier','hip_flexor_ii_multiplier','hip_extensor_pf_output_multiplier','ankle_dorsiflexor_pf_output_multiplier','hip_extensor_pf_activation_floor','knee_flexor_pf_output_multiplier'}==set(keys)|{'pf_output_scale'}
  assert all(type(value) in (int,float) and np.isfinite(value) for value in config.values()),'Central configuration values must be finite JSON numbers'
  self.settings=np.array([config[k] for k in keys]+([config['rg_reset_mV']] if 'rg_reset_mV' in config else []));self.encoded=np.array([i for i,n in enumerate(self.names) if self.params[i,3]>0 and (' RG ' in n or ' PF ' in n)],dtype=np.int64);assert len(self.v)==108 and len(self.encoded)==12
  self.pf_mask=np.array([' PF ' in self.names[i] for i in self.encoded]);P=len(self.encoded);N=self.population_size;self.lv=np.tile(self.v[self.encoded,None],(1,N))+.2*(np.arange(N)/(N-1)-.5);self.lm=np.tile(self.m[self.encoded,None],(1,N));self.lh=np.tile(self.h[self.encoded,None],(1,N));self.ref=np.zeros((P,N));self.adapt=np.zeros((P,N));self.decoded=np.maximum(0.,(self.v[self.encoded]-self.params[self.encoded,0])/20.);self.counts=np.zeros(P,dtype=np.int64);self.total_counts=self.counts.copy();self.logs=[];self.trace=[];self.output_edits=[]
  self.kernel=central_tick
  if self.reference_calculations:
   reference_globals=dict(central_tick.py_func.__globals__);reference_globals['K']=SimpleNamespace(neural_steps=K.neural_steps.py_func)
   self.kernel=FunctionType(central_tick.py_func.__code__,reference_globals,'central_tick_scalar_reference')
  self.tick(self.v.copy(),self.m.copy(),self.h.copy(),self.newv.copy(),self.gs.copy(),self.gi.copy(),self.lv.copy(),self.lm.copy(),self.lh.copy(),self.ref.copy(),self.adapt.copy(),self.decoded.copy(),self.counts.copy());SpikingCentralController.last_instance=self
  self.sensors=None
  self.ii_gain_edits=[]
  self.state_audit_result={}
  if self.sensory_spiking:
   from spiking_sensory import SensorPopulationBank
   self.sensors=SensorPopulationBank(self)
   assert not np.intersect1d(self.encoded,self.sensors.indices).size
   assert not any(idx in self.sensors.indices for idx,_,_,_ in self.stimuli+self.test_pulses),'Sensor-target current pulses need a separate encoder-current policy'
  self.balance=None
  if self.balance_configuration is not None:
   from balance_pathway import BodyBalanceRelay,LoadTransferTerminationRelay
   from stance_recruitment import StanceRecruitmentRelay
   self.balance=({'contralateral_load_termination':LoadTransferTerminationRelay,'phase_stance_recruitment':StanceRecruitmentRelay}.get(self.balance_configuration.get('pathway_kind'),BodyBalanceRelay))(self,self.balance_configuration)
  from hip_knee_interlock import build
  build(self)
 def apply_output_gain(self):
  if self.output_edits:return
  ii_scale=self.configuration.get('hip_flexor_ii_multiplier',1.);assert ii_scale>0
  if ii_scale!=1.:
   for e in self.edges:
    pre,post=self.names[int(e[0])],self.names[int(e[1])]
    if any((pre,post)==(side+' Hip flexor II relay',side+' Hip MN flx') for side in ('L','R')):
     before=float(e[2]);e[2]*=ii_scale;self.ii_gain_edits.append(dict(source=pre,target=post,before_gmax_uS=before,multiplier=ii_scale,effective_gmax_uS=float(e[2]),reversal_mV=float(e[5])))
   assert len(self.ii_gain_edits)==2
  for k,e in enumerate(self.edges):
   pre,post=self.names[int(e[0])],self.names[int(e[1])]
   if ' PF ' in pre and not pre.endswith(' IN') and ' MN ' in post and not post.endswith(' RE'):
    before=float(e[2]);scale=self.configuration['pf_output_scale']
    if ' Ank MN ' in post:scale*=self.configuration.get('ankle_pf_output_multiplier',1.)
    # Source Ank MN ext actuates the positive-dorsiflexion route on BOTH sides.
    if post in ('L Ank MN ext','R Ank MN ext'):scale*=self.configuration.get('ankle_dorsiflexor_pf_output_multiplier',1.)
    if post in ('L Knee MN flx','R Knee MN flx'):scale*=self.configuration.get('knee_flexor_pf_output_multiplier',1.)
    if post in ('L Hip MN ext','R Hip MN ext'):scale*=self.configuration.get('hip_extensor_pf_output_multiplier',1.)
    assert scale>0;self.edges[k,2]*=scale;edit=dict(source=pre,target=post,before_gmax_uS=before,scale=scale,effective_gmax_uS=float(self.edges[k,2]))
    floor=self.configuration.get('hip_extensor_pf_activation_floor',0.);assert 0<=floor<=.35
    if floor and post in ('L Hip MN ext','R Hip MN ext'):
     old=float(e[3]);e[3]=old+floor*(e[4]-old);edit.update(activation_floor=floor,before_threshold_mV=old,effective_threshold_mV=float(e[3]),saturation_mV=float(e[4]))
    self.output_edits.append(edit)
  assert len(self.output_edits)==12
 def tick(self,v,m,h,newv,gs,gi,lv,lm,lh,ref,adapt,decoded,counts):
  self.kernel(v,m,h,newv,gs,gi,self.params,self.edges,self.tonic,self.applied,self.encoded,self.pf_mask,self.settings,lv,lm,lh,ref,adapt,decoded,counts,self.dt)
 def export_state(self):
  """Complete local neural continuation state; external physics/time stay caller-owned."""
  state={'version':np.array(2),'configuration':np.array(json.dumps(self.configuration,sort_keys=True)),'output_edits':np.array(json.dumps(self.output_edits,sort_keys=True)),'ii_gain_edits':np.array(json.dumps(self.ii_gain_edits,sort_keys=True)),'sensor_present':np.array(self.sensors is not None)}
  for key in ('v','m','h','newv','gs','gi','applied','edges','tonic','lv','lm','lh','ref','adapt','decoded','counts','total_counts','last_test_impulse_nA_s'):state['central_'+key]=getattr(self,key).copy()
  for owner,prefix,fields in ((self,'central',('t','counts','decoded')),(self.sensors,'sensor',('t','counts','decoded','current'))):
   if owner is None:continue
   for i,key in enumerate(fields):state[prefix+'_log_'+key]=np.array([row[i] for row in owner.logs])
   for i,key in enumerate(('t','v','aux1','aux2') if owner is self else ('t','v','spikes')):state[prefix+'_trace_'+key]=np.array([row[i] for row in owner.trace])
  if self.sensors is not None:
   for key in ('v','ref','decoded','spikes','counts','total_counts','target','current'):state['sensor_'+key]=getattr(self.sensors,key).copy()
   state['sensor_proxy_present']=np.array(hasattr(self.sensors,'proxy'))
   if hasattr(self.sensors,'proxy'):state['sensor_proxy']=self.sensors.proxy.copy()
  state['balance_present']=np.array(self.balance is not None)
  state['clock_scale']=np.array(self.clock_scale)
  if self.balance is not None:
   for key,value in self.balance.export_state().items():state['balance_'+key]=value
  if self.heteronymous_load_indices:
   state['primitive_load_current_nA']=self.heteronymous_load_current_nA.copy();state['primitive_load_log_t']=np.array([v[0] for v in self.heteronymous_load_logs]);state['primitive_load_log_current_nA']=np.array([v[1] for v in self.heteronymous_load_logs])
  return state
 def import_state(self,state):
  assert int(state['version'])==2 and str(state['configuration'])==json.dumps(self.configuration,sort_keys=True)
  assert bool(state['sensor_present'])==(self.sensors is not None)
  for key in ('v','m','h','newv','gs','gi','applied','edges','tonic','lv','lm','lh','ref','adapt','decoded','counts','total_counts','last_test_impulse_nA_s'):np.copyto(getattr(self,key),state['central_'+key])
  self.output_edits=json.loads(str(state['output_edits']))
  self.ii_gain_edits=json.loads(str(state['ii_gain_edits']))
  for owner,prefix,fields in ((self,'central',('t','counts','decoded')),(self.sensors,'sensor',('t','counts','decoded','current'))):
   if owner is None:continue
   owner.logs=[tuple(float(values[0]) if i==0 else values[i].copy() for i in range(len(fields))) for values in zip(*(state[prefix+'_log_'+key] for key in fields))]
   trace_fields=('t','v','aux1','aux2') if owner is self else ('t','v','spikes');owner.trace=[tuple(float(values[0]) if i==0 else values[i].copy() for i in range(len(trace_fields))) for values in zip(*(state[prefix+'_trace_'+key] for key in trace_fields))]
  if self.sensors is not None:
   for key in ('v','ref','decoded','spikes','counts','total_counts','target','current'):np.copyto(getattr(self.sensors,key),state['sensor_'+key])
   if bool(state['sensor_proxy_present']):self.sensors.proxy=state['sensor_proxy'].copy()
   elif hasattr(self.sensors,'proxy'):del self.sensors.proxy
  assert bool(state.get('balance_present',False))==(self.balance is not None)
  assert float(state.get('clock_scale',1.))==self.clock_scale
  if self.balance is not None:self.balance.import_state({k[len('balance_'):]:v for k,v in state.items() if k.startswith('balance_') and k!='balance_present'})
  if self.heteronymous_load_indices:
   np.copyto(self.heteronymous_load_current_nA,state['primitive_load_current_nA']);self.heteronymous_load_logs=[(float(t),v.copy()) for t,v in zip(state['primitive_load_log_t'],state['primitive_load_log_current_nA'])]
 def observe_load_feedback(self,model,measured):
  """Explicit engineering normal-load afferent; previous actual physics frame."""
  if not self.heteronymous_load_indices:return
  from import_virtual_walker import mujoco
  bw=model.body_subtreemass[mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root')]*9.81
  f=np.zeros(2) if measured is None else np.maximum(0,np.asarray(measured['foot_force_world'])[:,1])/bw
  self.heteronymous_load_current_nA[:]=8*np.clip((f-.025)/.075,0,1)
 def step(self,time,ia,ii,ib,contacts,interval,feedback=True):
  if self.clock_scale!=1.:time/=self.clock_scale;interval/=self.clock_scale
  if self.audit_state and not self.state_audit_result and time>=2.:
   initial=self.export_state()
   for j in range(32):self._advance(time+j*interval,ia,ii,ib,contacts,interval,feedback)
   first=self.export_state();self.import_state(initial)
   for j in range(32):self._advance(time+j*interval,ia,ii,ib,contacts,interval,feedback)
   second=self.export_state();checks={k:bool(v.shape==second[k].shape and v.dtype==second[k].dtype and v.tobytes()==second[k].tobytes()) for k,v in first.items()};self.import_state(initial)
   self.state_audit_result=dict(all_arrays_byte_identical=all(checks.values()),arrays=checks,central_spikes_during_probe=int(first['central_total_counts'].sum()-initial['central_total_counts'].sum()),sensory_spikes_during_probe=int(first['sensor_total_counts'].sum()-initial['sensor_total_counts'].sum()) if self.sensors is not None else 0,scope='32 neural controller intervals with actual instantaneous physical inputs held fixed; all local neural states/history restored before the real physics interval. Physics, Hill states, and external time are caller-owned. No RNG.')
   assert all(checks.values()),'Neural continuation-state audit failed'
  result=self._advance(time,ia,ii,ib,contacts,interval,feedback)
  if self.clock_scale!=1.:self.last_test_impulse_nA_s*=self.clock_scale
  return result
 def _advance(self,time,ia,ii,ib,contacts,interval,feedback=True):
  self.apply_output_gain()  # After the caller's source recruitment calibration.
  self.applied[:]=0.;self.last_test_impulse_nA_s[:]=0.;self.counts[:]=0
  if feedback:
   for target,idx,kind,coeff,body in self.inputs:
    if kind=='Ia':value=ia[idx]
    elif kind=='II':value=ii[idx]
    elif kind=='Ib':value=ib[idx]
    elif kind=='Tension':value=ib[idx]/self.ibgain[idx]
    elif kind in ('ContactCount','Contacts'):value=contacts.get(body,0.)
    else:raise ValueError(kind)
    self.applied[target]+=K.polynomial(coeff,value)*1e9
  if self.heteronymous_load_indices and feedback:self.applied[self.heteronymous_load_indices]+=self.heteronymous_load_current_nA
  ticks=round(interval/self.dt);assert abs(ticks*self.dt-interval)<1e-12
  if self.sensors is not None:self.sensors.prepare()
  if self.balance is not None:self.balance.prepare()
  for sub in range(ticks):
   t=time+sub*self.dt
   for idx,start,end,amplitude in self.stimuli:
    if start<=t<end:self.applied[idx]+=amplitude
   for idx,start,end,amplitude in self.test_pulses:
    if start<=t<end:self.applied[idx]+=amplitude;self.last_test_impulse_nA_s[idx]+=amplitude*self.dt
   if self.sensors is not None:self.sensors.advance(t)
   if self.balance is not None:self.balance.advance()
   self.tick(self.v,self.m,self.h,self.newv,self.gs,self.gi,self.lv,self.lm,self.lh,self.ref,self.adapt,self.decoded,self.counts)
   if self.sensors is not None:self.v[self.sensors.indices]=self.sensors.proxy
   if self.balance is not None:self.v[self.balance.indices]=self.balance.proxy
   if 2.<=t<2.25:self.trace.append((t,self.lv[:,:2].copy(),self.adapt.mean(axis=1),self.lh.mean(axis=1)))
   for idx,start,end,amplitude in self.test_pulses:
    if start<=t<end:self.applied[idx]-=amplitude
   for idx,start,end,amplitude in self.stimuli:
    if start<=t<end:self.applied[idx]-=amplitude
  if self.sensors is not None:self.sensors.finish(time)
  if self.balance is not None:self.balance.finish(time)
  if self.heteronymous_load_indices:self.heteronymous_load_logs.append((time,self.heteronymous_load_current_nA.copy() if feedback else np.zeros(2)))
  self.total_counts+=self.counts;self.logs.append((time,self.counts.copy(),self.decoded.copy()));assert np.isfinite(self.v).all() and np.isfinite(self.lv).all();return self.v
 def rg_current_terms(self):
  rows=super().rg_current_terms()
  for row in rows:
   idx=self.names.index(row['target']);p=np.flatnonzero(self.encoded==idx)
   if not len(p):continue
   p=int(p[0]);mean=float(self.lv[p].mean());ca_scale=self.settings[0]*(self.settings[7] if self.pf_mask[p] else 1.);row['scope']='Instantaneous mean per spiking cell after reset; proxy voltage is not soma voltage';row['decoded_proxy_mV']=row['voltage_mV'];row['voltage_mV']=mean;row['leak_nA']=float(self.params[idx,0]-mean);row['calcium_nA']=float(np.mean(self.params[idx,3]*ca_scale*self.lm[p]*self.lh[p]*(self.params[idx,10]-self.lv[p])));row['adaptation_nA']=-float(self.adapt[p].mean())
   for edge in row['incoming']:edge['current_nA']=float(edge['g_active']*(edge['reversal_mV']-mean))
   row['synaptic_nA']=sum(x['current_nA'] for x in row['incoming']);row['net_nA']=sum(row[k] for k in ('leak_nA','calcium_nA','adaptation_nA','tonic_nA','physical_applied_nA','synaptic_nA'))
  return rows
 def metadata(self):
  meta=dict(architecture='AdaptiveCaLIF_RG_PF',scope='Twelve RG/PF source cells replaced by actual Ca/LIF populations with explicit inhibitory adaptation. '+('Sixteen II/contact channels also replaced by 512 LIF cells; 80 other source cells remain graded. Only spike-filtered outputs drive the replaced central/sensory chemical pathways.' if self.sensors is not None else 'Other 96 source cells, including Ca-bearing Renshaw cells and MN, remain graded.')+' Proxy voltage is decoded spike activity, not soma voltage.',spiking_populations=len(self.encoded),central_spiking_populations=len(self.encoded),sensory_spiking_populations=16 if self.sensors is not None else 0,sensory_cells_per_population=32 if self.sensors is not None else 0,total_spiking_populations=len(self.encoded)+(16 if self.sensors is not None else 0),cells_per_population=self.population_size,actual_spiking_cells=len(self.encoded)*self.population_size+(512 if self.sensors is not None else 0),remaining_graded_cells=len(self.v)-len(self.encoded)-(16 if self.sensors is not None else 0),actual_total_cells=len(self.v)-len(self.encoded)+len(self.encoded)*self.population_size+(496 if self.sensors is not None else 0),state_snapshot_verification=self.state_audit_result,sensory_spiking=self.sensors is not None,sensory_spike_counts=None if self.sensors is None else self.sensors.total_counts.tolist(),sensory_spiking_names=[] if self.sensors is None else [self.names[i] for i in self.sensors.indices],homonymous_II_gain_edits=self.ii_gain_edits,spike_count_per_population=self.total_counts.tolist(),spiking_names=[self.names[i] for i in self.encoded],configuration=self.configuration,output_gain_edits=self.output_edits,startup_stimuli=self.stimuli)
  if self.balance is not None:
   termination=self.balance.cfg.get('pathway_kind')=='contralateral_load_termination';scope=' Contralateral force-qualified transient load events excite inherited PF-flexor INs to terminate opposite extensor bursts; no direct RG/MN/muscle command.' if termination else ' Separate phase-gated body-balance populations excite configured joint MNs and inherited inhibitory Ia cells.'
   stance=self.balance.cfg.get('pathway_kind')=='phase_stance_recruitment'
   if stance:scope=' Separate source-RG-phase and actual-load-qualified spiking PF stance populations supply fixed mirrored excitatory/inhibitory chemical recruitment to MNs; existing transient PF axons stay unchanged. Offline source-law fit, no online allocator.'
   meta.update(architecture='AdaptiveCaLIF_RG_PF_LoadTransferTermination' if termination else 'AdaptiveCaLIF_RG_PF_BodyBalance',actual_spiking_cells=meta['actual_spiking_cells']+32,remaining_graded_cells=meta['remaining_graded_cells']-2,actual_total_cells=meta['actual_total_cells']+30,body_balance_populations=2,total_spiking_populations=meta['total_spiking_populations']+2,body_balance_configuration=self.balance.cfg,body_balance_spikes=self.balance.total_counts.tolist(),runtime_logical_channels=len(self.v),runtime_chemical_edges=len(self.edges),scope=meta['scope']+scope+f' Adds32 actual LIF cells and{len(self.balance.edge_specs)} chemical edges; source graph remains108/216.')
   if stance:meta['architecture']='AdaptiveCaLIF_RG_PF_StanceRecruitment'
  if self.heteronymous_edges:meta.update(architecture=meta['architecture']+'_HipII_KneeInterlock',heteronymous_edges=self.heteronymous_edits,extra_hipII_spiking_PF_names=[self.names[i] for i in self.heteronymous_primitive_indices],runtime_logical_channels=len(self.v),runtime_chemical_edges=len(self.edges),scope=meta['scope']+f' Adds {self.heteronymous_nodes} logical channels, including {len(self.heteronymous_primitive_indices)} adapting-LIF PF populations; remaining added cells are graded inhibitory INs. Spike-driven hip-extensor-II inhibits knee-extensor MNs and optionally recruits knee-flexor MNs directly or through the new spiking PF primitive. Extra primitive output gmax is an absolute separately calibrated value, not multiplied by the original twelve-PF output scale. No source Ia RG/PF collateral changes or mechanical commands.')
  if self.heteronymous_load_indices:meta.update(primitive_load_feedback=dict(input='Previous actual ipsilateral foot normal GRF, not geometric ContactCount or tendon Ib',current_map_nA='8*clip((normalF/BW−0.025)/0.075,0,1)',input_ceiling_nA=8.,readout='Two added graded inhibitory INs inhibit only the extra spiking PF populations',source_original_contact_routes_unchanged=True))
  if self.clock_scale!=1.:meta.update(neural_clock_multiplier=self.clock_scale,neural_physical_substep_s=self.dt*self.clock_scale,source_startup_stimuli_biological_time=meta['startup_stimuli'],startup_stimuli=[(i,start*self.clock_scale,end*self.clock_scale,current) for i,start,end,current in self.stimuli],clock_scope='Whole neural biological clock rescaled; source currents, conductances, reset potentials, and readout normalization unchanged. Physics and Hill clock unchanged. Not an RG-only kinetics intervention or computational speedup.')
  return meta
 def save(self,folder):
  if self.sensors is not None:self.sensors.save(folder)
  if self.balance is not None:self.balance.save(folder)
  if self.heteronymous_load_indices:np.savez_compressed(folder/'primitive_load_input.npz',t=np.array([v[0] for v in self.heteronymous_load_logs]),current_nA=np.array([v[1] for v in self.heteronymous_load_logs]),names=[self.names[i] for i in self.heteronymous_load_indices])
  np.savez_compressed(folder/'central_spikes.npz',t=np.array([x[0] for x in self.logs])*self.clock_scale,spike_counts=np.array([x[1] for x in self.logs]),decoded_activation=np.array([x[2] for x in self.logs]),names=[self.names[i] for i in self.encoded],trace_t=np.array([x[0] for x in self.trace])*self.clock_scale,trace_membrane_mV=np.array([x[1] for x in self.trace]),trace_adaptation_nA=np.array([x[2] for x in self.trace]),trace_h_mean=np.array([x[3] for x in self.trace]))
