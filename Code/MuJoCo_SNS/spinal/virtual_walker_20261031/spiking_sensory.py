"""Actual reset/refractory LIF sensory populations in the source108 graded circuit.

This is explicitly a sensory-spiking hybrid: source RG/PF/IN/MN dynamics stay
graded. Sixteen input-only II/contact encoders are replaced, not shadowed by
graded encoder dynamics. Filtered spike counts supply their outgoing synapses.
"""
import math
import numpy as np
from numba import njit
import native_kernel as K
BaseController=K.NativeController

@njit(cache=True)
def encode_step(target,v,refractory,decoded,spikes,counts,dt,rate_max,tau_m,tau_syn,refractory_s):
 leak=math.exp(-dt/tau_m);pop=v.shape[1]
 for p in range(v.shape[0]):
  decay=math.exp(-dt/tau_syn[p])
  rate=rate_max*min(1.,max(0.,target[p]));drive=0.
  if rate>0:drive=1./(-math.expm1(-(1./rate-refractory_s)/tau_m))
  count=0
  for j in range(pop):
   spikes[p,j]=0
   if refractory[p,j]>0:refractory[p,j]=max(0.,refractory[p,j]-dt);continue
   v[p,j]=drive+(v[p,j]-drive)*leak
   if v[p,j]>=1.:
    v[p,j]=0.;refractory[p,j]=refractory_s;spikes[p,j]=1;count+=1
  decoded[p]=decoded[p]*decay+count/(pop*rate_max*tau_syn[p])
  counts[p]+=count

class SensorPopulationBank:
 """The verified sixteen input encoders, reusable beside a spiking central core."""
 def __init__(self,controller):
  self.c=controller;self.indices=np.array([i for i,n in enumerate(controller.names) if n.endswith(' II relay') or n in ('L heel contact','L toe contact','R heel contact','R toe contact')],dtype=np.int64);assert len(self.indices)==16
  for i in self.indices:
   assert not np.any(controller.edges[:,1]==i) and controller.params[i,3]==0 and controller.params[i,2]>0 and controller.tonic[i]==0
   edges=controller.edges[controller.edges[:,0]==i];assert np.all(edges[:,3]==-60) and np.all(edges[:,4]==-40)
  self.v=np.tile(np.linspace(.01,.99,32),(16,1));self.ref=np.zeros((16,32));self.decoded=np.zeros(16);self.spikes=np.zeros((16,32),dtype=np.int8);self.counts=np.zeros(16,dtype=np.int64);self.total_counts=self.counts.copy();self.target=np.zeros(16);self.current=np.zeros(16);self.tau=controller.params[self.indices,1]*.001;self.logs=[];self.trace=[]
  encode_step(self.target,self.v.copy(),self.ref.copy(),self.decoded.copy(),self.spikes.copy(),self.counts.copy(),controller.dt,200.,.005,self.tau,.001)
 def prepare(self):
  self.counts[:]=0;self.current[:]=self.c.applied[self.indices];self.target[:]=np.maximum(0.,self.current/(20.*self.c.params[self.indices,2]));self.c.applied[self.indices]=0.
 def advance(self,time):
  encode_step(self.target,self.v,self.ref,self.decoded,self.spikes,self.counts,self.c.dt,200.,.005,self.tau,.001);self.proxy=self.c.params[self.indices,0]+20.*self.decoded;self.c.v[self.indices]=self.proxy
  if 2.<=time<2.25:self.trace.append((time,self.v[:,:2].copy(),self.spikes[:,:2].copy()))
 def finish(self,time):
  self.total_counts+=self.counts;self.logs.append((time,self.counts.copy(),self.decoded.copy(),self.current.copy()))
 def save(self,folder):
  np.savez_compressed(folder/'sensory_spikes.npz',t=np.array([r[0] for r in self.logs])*getattr(getattr(self,'c',self),'clock_scale',1.),spike_counts=np.array([r[1] for r in self.logs]),decoded_activation=np.array([r[2] for r in self.logs]),physical_input_nA=np.array([r[3] for r in self.logs]),names=[self.c.names[i] for i in self.indices],trace_t=np.array([r[0] for r in self.trace])*getattr(getattr(self,'c',self),'clock_scale',1.),trace_membrane_mV=np.array([r[1] for r in self.trace])*20.-70.,trace_spikes=np.array([r[2] for r in self.trace]))

class SpikingSensoryController(BaseController):
 population_size=32;rate_max=200.;tau_m=.005;refractory_s=.001;tau_syn_scale=1.;ablate_outputs=False;last_instance=None
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs);self.encoded=np.array([i for i,n in enumerate(self.names) if n.endswith(' II relay') or n in ('L heel contact','L toe contact','R heel contact','R toe contact')],dtype=np.int64)
  assert len(self.encoded)==16 and len(self.v)==108
  for i in self.encoded:
   assert not np.any(self.edges[:,1]==i),'Only input-only sensory cells can use this encoder replacement'
   assert self.params[i,3]==0 and self.params[i,2]>0 and self.tonic[i]==0
   e=self.edges[self.edges[:,0]==i];assert np.all(e[:,3]==-60) and np.all(e[:,4]==-40)
  P=len(self.encoded);N=self.population_size;self.lif_v=np.tile(np.linspace(.01,.99,N),(P,1));self.refractory=np.zeros((P,N));self.decoded=np.zeros(P);self.spikes=np.zeros((P,N),dtype=np.int8);self.counts=np.zeros(P,dtype=np.int64);self.targets=np.zeros(P);self.original_current=np.zeros(P);self.tau_syn=self.params[self.encoded,1]*.001*self.tau_syn_scale;self.logs=[];self.trace=[];SpikingSensoryController.last_instance=self
  encode_step(self.targets,self.lif_v.copy(),self.refractory.copy(),self.decoded.copy(),self.spikes.copy(),self.counts.copy(),self.dt,self.rate_max,self.tau_m,self.tau_syn,self.refractory_s)
 def step(self,time,ia,ii,ib,contacts,interval,feedback=True):
  if self.ablate_outputs and not hasattr(self,'ablated_edges'):
   indices=np.flatnonzero(np.isin(self.edges[:,0].astype(np.int64),self.encoded));self.ablated_edges=[dict(source=self.names[int(self.edges[k,0])],target=self.names[int(self.edges[k,1])],before_gmax_uS=float(self.edges[k,2]),after_gmax_uS=0.) for k in indices];self.edges[indices,2]=0.
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
  count=round(interval/self.dt);assert abs(count*self.dt-interval)<1e-12
  self.original_current[:]=self.applied[self.encoded];self.targets[:]=np.maximum(0.,self.original_current/(20.*self.params[self.encoded,2]));self.applied[self.encoded]=0.
  for sub in range(count):
   t=time+sub*self.dt
   for idx,start,end,amplitude in self.stimuli:
    if start<=t<end:self.applied[idx]+=amplitude
   for idx,start,end,amplitude in self.test_pulses:
    if start<=t<end:self.applied[idx]+=amplitude;self.last_test_impulse_nA_s[idx]+=amplitude*self.dt
   encode_step(self.targets,self.lif_v,self.refractory,self.decoded,self.spikes,self.counts,self.dt,self.rate_max,self.tau_m,self.tau_syn,self.refractory_s)
   proxy=self.params[self.encoded,0]+20.*self.decoded;self.v[self.encoded]=proxy
   K.neural_steps(self.v,self.m,self.h,self.newv,self.gs,self.gi,self.params,self.edges,self.tonic,self.applied,self.dt*1000,1)
   self.v[self.encoded]=proxy
   if 2.<=t<2.25:self.trace.append((t,self.lif_v[:,:2].copy(),self.spikes[:,:2].copy()))
   for idx,start,end,amplitude in self.test_pulses:
    if start<=t<end:self.applied[idx]-=amplitude
   for idx,start,end,amplitude in self.stimuli:
    if start<=t<end:self.applied[idx]-=amplitude
  self.logs.append((time,self.counts.copy(),self.decoded.copy(),self.original_current.copy()));assert np.all(self.v<40);return self.v
 def save(self,folder):
  np.savez_compressed(folder/'spiking_encoders.npz',t=np.array([r[0] for r in self.logs])*getattr(getattr(self,'c',self),'clock_scale',1.),spike_counts=np.array([r[1] for r in self.logs]),decoded_activation=np.array([r[2] for r in self.logs]),physical_input_nA=np.array([r[3] for r in self.logs]),names=[self.names[i] for i in self.encoded],trace_t=np.array([r[0] for r in self.trace])*getattr(getattr(self,'c',self),'clock_scale',1.),trace_membrane_mV=np.array([r[1] for r in self.trace])*20.-70.,trace_spikes=np.array([r[2] for r in self.trace]))
