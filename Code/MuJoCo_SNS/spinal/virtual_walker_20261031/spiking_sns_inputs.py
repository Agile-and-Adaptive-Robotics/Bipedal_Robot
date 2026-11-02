"""Opt-in source-graph pure-input LIF population encoders for SNS sparse models."""
import math
import numpy as np
from spiking_sensory import encode_step

class SpikeInputForward:
 def __init__(self,inner,names,population_size=16,rate_max=200.):
  self.inner=inner;raw=inner.c;self.names=names;self.population_size=population_size;self.rate_max=rate_max
  incoming=set(inner.post[inner.gmax!=0].tolist());indices=[];ports=[]
  for i,name in enumerate(names):
   if not name.upper().startswith(('IA_','II_','IB_','HEEL_','TOE_','TOEDF_','HUNT_AF_')) or i in incoming:continue
   p=np.flatnonzero(raw.input_connectivity[i])
   if len(p)!=1 or np.count_nonzero(raw.input_connectivity[:,p[0]])!=1:continue
   if raw.gated and np.any(raw.g_ion[:,i]!=0):continue
   indices.append(i);ports.append(int(p[0]))
  self.indices=np.array(indices,dtype=np.int64);self.ports=np.array(ports,dtype=np.int64);assert len(indices)>0,'No source pure-input sensory encoders found'
  P=len(indices);N=population_size;self.v=np.tile(np.linspace(.01,.99,N),(P,1));self.ref=np.zeros((P,N));self.decoded=np.zeros(P);self.spikes=np.zeros((P,N),dtype=np.int8);self.counts=np.zeros(P,dtype=np.int64);self.total_counts=np.zeros(P,dtype=np.int64);self.target=np.zeros(P);self.logs=[];self.steps=0;self.dt=.0002;self.substeps=round(raw.dt/self.dt);assert abs(self.substeps*self.dt-raw.dt)<1e-12
  decay=1.-raw.time_factor_membrane[self.indices]*raw.g_m[self.indices];assert np.all((decay>0)&(decay<1)),'Instantaneous/unstable source inputs need an explicit encoder policy'
  self.tau= -raw.dt/np.log(decay);self.input_gain=raw.input_connectivity[self.indices,self.ports];self.span=np.full(P,5.);self.rest=raw.V_rest[self.indices].copy()
  encode_step(self.target,self.v.copy(),self.ref.copy(),self.decoded.copy(),self.spikes.copy(),self.counts.copy(),self.dt,self.rate_max,.005,self.tau,.001)
 def __getattr__(self,name):return getattr(self.inner,name)
 def export_state(self):
  """Numeric continuation state of this encoder and its graded sparse backend.

  External mechanics, actuator states, and network/runner counters are owned by
  the caller. There is no RNG: initial population phases are deterministic.
  """
  raw=self.inner.c;assert not raw.spiking,'Snapshot policy is for the graded sparse backend'
  state={'encoder_version':np.array(1,dtype=np.int64),'encoder_steps':np.array(self.steps,dtype=np.int64),'encoder_names':np.array(self.names),'encoder_indices':self.indices.copy(),'encoder_ports':self.ports.copy(),'encoder_parameters':np.array([self.population_size,self.rate_max,self.dt,self.substeps]),'encoder_logs':np.asarray(self.logs,dtype=np.int64).reshape(-1,len(self.indices)),'encoder_rng_state':np.empty(0,dtype=np.uint8)}
  for key in ('v','ref','decoded','spikes','counts','total_counts','target'):state['encoder_'+key]=getattr(self,key).copy()
  for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last','outputs'):
   if hasattr(raw,key):state['backend_'+key]=np.asarray(getattr(raw,key)).copy()
  for key in ('last_voltage','input_values','i_app','output_values','sums','acc','values'):
   if hasattr(self.inner,key):state['sparse_'+key]=np.asarray(getattr(self.inner,key)).copy()
  return state
 def import_state(self,state):
  """Restore a snapshot only into the same source ordering and encoder policy."""
  assert int(state['encoder_version'])==1
  assert np.array_equal(state['encoder_names'],self.names)
  assert np.array_equal(state['encoder_indices'],self.indices) and np.array_equal(state['encoder_ports'],self.ports)
  assert np.array_equal(state['encoder_parameters'],[self.population_size,self.rate_max,self.dt,self.substeps]) and len(state['encoder_rng_state'])==0
  for key in ('v','ref','decoded','spikes','counts','total_counts','target'):
   value=state['encoder_'+key];current=getattr(self,key);assert value.shape==current.shape and value.dtype==current.dtype;np.copyto(current,value)
  self.steps=int(state['encoder_steps']);self.logs=[row.copy() for row in state['encoder_logs']]
  for key in ('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last','outputs'):
   if 'backend_'+key in state:setattr(self.inner.c,key,state['backend_'+key].copy())
  for key in ('last_voltage','input_values','i_app','output_values','sums','acc','values'):
   if 'sparse_'+key in state:np.copyto(getattr(self.inner,key),state['sparse_'+key])
  self.inner.c.V_last=self.inner.last_voltage
  if self.inner.output_onehot:self.inner.c.outputs=self.inner.output_values
 def forward(self,x=None):
  raw=self.inner.c;assert x is not None;u=np.array(x,copy=True)
  physical=u[self.ports]*self.input_gain;self.target[:]=np.maximum(0.,(physical+raw.i_b[self.indices])/(raw.g_m[self.indices]*self.span));u[self.ports]=0.
  raw.V[self.indices]=self.rest+self.span*self.decoded
  result=self.inner.forward(u)
  self.counts[:]=0
  for _ in range(self.substeps):encode_step(self.target,self.v,self.ref,self.decoded,self.spikes,self.counts,self.dt,self.rate_max,.005,self.tau,.001)
  raw.V[self.indices]=self.rest+self.span*self.decoded;self.total_counts+=self.counts;self.logs.append(self.counts.copy());self.steps+=1
  if self.inner.output_onehot:np.take(raw.V,self.inner.oq,out=self.inner.output_values);raw.outputs=self.inner.output_values
  else:raw.outputs=np.matmul(raw.output_voltage_connectivity,raw.V)
  return raw.outputs
 def metadata(self):return dict(scope='Pure-input sensory populations only; source central RG/PF/MN remain graded. Source input-only encoders replaced by reset/refractory LIF populations; spike-filter time constants match source discrete membrane decay.',spiking_source_names=[self.names[i] for i in self.indices],source_indices=self.indices.tolist(),cells_per_population=self.population_size,actual_spiking_cells=len(self.indices)*self.population_size,remaining_graded_cells=len(self.names)-len(self.indices),spike_counts=self.total_counts.tolist(),rate_max_hz=self.rate_max,neural_spike_dt_s=self.dt,synaptic_filter_tau_s=self.tau.tolist())
