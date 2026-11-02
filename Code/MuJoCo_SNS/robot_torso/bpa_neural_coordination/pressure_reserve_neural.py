"""Explicit bilateral engineering reserve relays, not biological spindle replicas.

Two leaky relay cells receive a pressure/route reserve predictor and excite the
homonymous pectineus MN through bounded conductance currents. Existing POST
input ports are used solely as the compiled additive-current transport; their
original tonic values are preserved and all extra currents are logged.
No motor voltage/output, force, pressure state, or geometric parameter is edited.
"""
import numpy as np

CHANNELS=['reserve_m','reserve_rate_m_s','predicted_reserve_m','sensor_normalized','relay_voltage_mV','added_MN_current_nA','pressure_kpa','route_length_m','route_velocity_m_s','pressure_rate_kpa_s']
STATE_KEYS=('V','V_last','b_gate','c_gate','b_gate_last','c_gate_last','outputs')
SPARSE_KEYS=('last_voltage','input_values','i_app','output_values','sums','acc','values')

class ReserveForward:
 def __init__(self,inner,names,context,current_cap_nA,preview_s):
  self.inner=inner;self.names=list(names);self.context=context;self.cap=float(current_cap_nA);self.preview=float(preview_s);self.dt=.002;self.tau=.006;self.derivative_tau=.010;self.force_scale=1.;self.reversal=40.
  if not 0<=self.cap<=.5 or not .02<=self.preview<=.2:raise ValueError('Outside predeclared bounded relay range')
  self.routes=['pect_r','pect_l'];self.mn=np.array([self.names.index('MN_'+a) for a in self.routes]);ports=[]
  for i in self.mn:
   p=np.flatnonzero(inner.c.input_connectivity[i])
   if len(p)!=1 or inner.c.input_connectivity[i,p[0]]!=1 or np.count_nonzero(inner.c.input_connectivity[:,p[0]])!=1:raise ValueError('Require dedicated unit-weight direct MN input port')
   ports.append(int(p[0]))
  self.ports=np.array(ports);self.relay_v=np.zeros(2);self.previous_pressure=np.zeros(2);self.pressure_rate=np.zeros(2);self.initialized=False;self.steps=0;self.logs=[];self.restore_verified=False;self.input_buffer=np.zeros(inner.c.input_connectivity.shape[1])
 def __getattr__(self,key):return getattr(self.inner,key)
 def _step(self,x):
  row=np.zeros((2,len(CHANNELS)));np.copyto(self.input_buffer,x)
  if self.context.get('active',False) and self.cap>0:
   adapter=self.context['adapter'];data=self.context['data'];idx=np.array([adapter.bank.names.index(a) for a in self.routes]);ids=adapter.ids[idx];p=adapter.bank.parameters[idx];pressure=adapter.bank.pressure[idx];length=data.actuator_length[ids];velocity=data.actuator_velocity[ids]
   if self.initialized:raw_rate=(pressure-self.previous_pressure)/self.dt;self.pressure_rate+=self.dt/self.derivative_tau*(raw_rate-self.pressure_rate)
   else:self.initialized=True
   self.previous_pressure[:]=pressure;factor=p[:,3]/(620.*p[:,4]);reserve=p[:,0]+p[:,2]+factor*pressure-length;rate=factor*self.pressure_rate-velocity;predicted=reserve+self.preview*np.minimum(rate,0.);signal=np.clip(-predicted*p[:,4]/self.force_scale,0.,1.)
   current=np.clip((self.cap/self.reversal)*np.clip(self.relay_v/10.,0,1)*(self.reversal-self.inner.c.V[self.mn]),0.,self.cap)
   self.input_buffer[self.ports]+=current
   old_v=self.relay_v.copy();self.relay_v+=self.dt/self.tau*(10.*signal-self.relay_v)
   row[:]=np.column_stack((reserve,rate,predicted,signal,old_v,current,pressure,length,velocity,self.pressure_rate))
  output=self.inner.forward(x if self.cap==0 else self.input_buffer);self.logs.append(row);self.steps+=1;return output
 def forward(self,x=None):
  if x is None:raise ValueError('Reserve relay needs explicit source inputs')
  if self.cap>0 and self.context.get('active',False) and not self.restore_verified and np.max(self.relay_v)>1e-8:
   initial=self.export_state();a=self._step(x).copy();end=self.export_state();self.import_state(initial);b=self._step(x);repeat=self.export_state()
   assert np.array_equal(a,b)
   for k in end:assert np.array_equal(end[k],repeat[k]),'Complete neural restore mismatch: '+k
   self.restore_verified=True;return b
  return self._step(x)
 def export_state(self):
  state=dict(reserve_version=np.array(1),reserve_parameters=np.array([self.cap,self.preview,self.dt,self.tau,self.derivative_tau,self.force_scale,self.reversal]),reserve_names=np.array(self.routes),reserve_ports=self.ports.copy(),reserve_relay_voltage=self.relay_v.copy(),reserve_previous_pressure=self.previous_pressure.copy(),reserve_pressure_rate=self.pressure_rate.copy(),reserve_initialized=np.array(self.initialized),reserve_steps=np.array(self.steps),reserve_logs=np.asarray(self.logs).reshape(-1,2,len(CHANNELS)),reserve_restore_verified=np.array(self.restore_verified),reserve_active=np.array(self.context.get('active',False)),reserve_rng_state=np.empty(0,dtype=np.uint8),reserve_input_buffer=self.input_buffer.copy())
  if callable(getattr(self.inner,'export_state',None)):
   for k,v in self.inner.export_state().items():state['inner_'+k]=v
  else:
   for k in STATE_KEYS:
    if hasattr(self.inner.c,k):state['backend_'+k]=np.asarray(getattr(self.inner.c,k)).copy()
   for k in SPARSE_KEYS:
    if hasattr(self.inner,k):state['sparse_'+k]=np.asarray(getattr(self.inner,k)).copy()
  return state
 def import_state(self,state):
  assert np.array_equal(state['reserve_parameters'],[self.cap,self.preview,self.dt,self.tau,self.derivative_tau,self.force_scale,self.reversal]) and np.array_equal(state['reserve_names'],self.routes) and len(state['reserve_rng_state'])==0
  for key,target in (('reserve_relay_voltage',self.relay_v),('reserve_previous_pressure',self.previous_pressure),('reserve_pressure_rate',self.pressure_rate),('reserve_input_buffer',self.input_buffer)):np.copyto(target,state[key])
  self.initialized=bool(state['reserve_initialized']);self.steps=int(state['reserve_steps']);self.logs=[row.copy() for row in state['reserve_logs']];self.restore_verified=bool(state['reserve_restore_verified']);self.context['active']=bool(state['reserve_active'])
  if callable(getattr(self.inner,'import_state',None)):self.inner.import_state({k[6:]:v for k,v in state.items() if k.startswith('inner_')})
  else:
   for k in STATE_KEYS:
    if 'backend_'+k in state:setattr(self.inner.c,k,state['backend_'+k].copy())
   for k in SPARSE_KEYS:
    if 'sparse_'+k in state:np.copyto(getattr(self.inner,k),state['sparse_'+k])
 def circuit(self):
  return dict(architecture='Two named engineering pressure-reserve leaky relays with homonymous excitatory MN conductance connections; original tonic input is preserved',neurons=[dict(name='RESERVE_'+a,tau_ms=6.,leak_uS=1.,capacitance_nF=6.,sensor_fullscale_nA=10.) for a in self.routes],edges=[dict(source='RESERVE_'+a,target='MN_'+a,gmax_uS=self.cap/40.,reversal_mV=40.,current_cap_nA=self.cap,compiled_additive_input_port=int(port)) for a,port in zip(self.routes,self.ports)],sensor='h=l0+b+Fmax_nominal*Pactual/(620K)-L; filtered pressure derivative10ms; hdot=Fmax_nominal*Pdot/(620K)-route_velocity; signal=clip(-K*(h+preview*min(hdot,0))/1N,0,1)',preview_s=self.preview,bilateral_tying=True,oracle_actual_force_gain_used=False)

def reserve_adapter(base,module,context):
 class Adapter(base):
  def __init__(self,*args,**kwargs):super().__init__(*args,**kwargs);context['adapter']=self
  def advance(self,data,commands):context['data']=data;return super().advance(data,commands)
  def configure_response(self,cfg):super().configure_response(cfg);context['active']=True
 module.NativeBPAForce=Adapter
