"""Shared neural pressure-capacity input weight or measured-tension observer.

Actual actuator parameters are never changed. The pressure mode's learned gain
acts only on its sensory/internal capacity estimate. Force mode uses the same
positive tension available to the existing model Ib channel, with declared new
transport delay and noise. Neither receives the actual force-map uncertainty.
"""
import json,hashlib
import numpy as np
from pressure_valve_neural import ValveReserveForward
from pressure_reserve_neural import CHANNELS

SENSOR_CHANNELS=['measured_force_N','force_rate_N_s','pressure_capacity_reserve_m','force_series_reserve_m','neural_reference_increment_kpa','sensor_noise_pressure_kpa','sensor_noise_force_N','additional_transport_delay_s']

class CapacityReserveForward(ValveReserveForward):
 def __init__(self,inner,names,context,pressure_cap_kpa,preview_s,reserve_target_N=1.,mode='pressure',capacity_gain=1.,delay_steps=0,pressure_noise_kpa=0.,force_noise_N=0.,noise_seed=0):
  super().__init__(inner,names,context,pressure_cap_kpa,preview_s,reserve_target_N)
  if mode not in ('pressure','force') or not .7<=capacity_gain<=1. or not 0<=delay_steps<=10:raise ValueError('Invalid estimator architecture parameters')
  if mode=='force' and capacity_gain!=1.:raise ValueError('Force estimator must not consume unused capacity gain')
  self.mode=mode;self.capacity_gain=float(capacity_gain);self.delay_steps=int(delay_steps);self.pressure_noise=float(pressure_noise_kpa);self.force_noise=float(force_noise_N);self.noise_seed=int(noise_seed)
  rng=np.random.default_rng(self.noise_seed);noise=np.clip(rng.standard_normal((30000,2,2)),-3.,3.);noise[:,:,0]*=self.pressure_noise;noise[:,:,1]*=self.force_noise;self.noise=noise;self.noise_hash=hashlib.sha256(noise.tobytes()).hexdigest()
  self.buffer=np.zeros((self.delay_steps+1,2,4));self.buffer_count=0;self.buffer_initialized=False;self.noise_index=0;self.previous_force=np.zeros(2);self.force_rate=np.zeros(2);self.sensor_logs=[]
 def _step(self,x):
  row=np.zeros((2,len(CHANNELS)));sensor=np.zeros((2,len(SENSOR_CHANNELS)));np.copyto(self.input_buffer,x);self.increment[:]=0.;observed=None;p=None
  if 'adapter' in self.context and 'data' in self.context:
   adapter=self.context['adapter'];data=self.context['data'];idx=np.array([adapter.bank.names.index(a) for a in self.routes]);ids=adapter.ids[idx];p=adapter.bank.parameters[idx]
   raw=np.column_stack((adapter.bank.pressure[idx],adapter.bank.force[idx],data.actuator_length[ids],data.actuator_velocity[ids]))
   if not self.buffer_initialized:self.buffer[:]=raw;self.buffer_initialized=True
   self.buffer[self.buffer_count%len(self.buffer)]=raw;observed=self.buffer[(self.buffer_count-self.delay_steps)%len(self.buffer)].copy();self.buffer_count+=1
  if self.context.get('active',False) and self.pressure_cap>0:
   if observed is None:raise ValueError('Estimator has no physical measurement')
   if self.noise_index>=len(self.noise):raise ValueError('Declared sensor-noise sequence exhausted')
   error=self.noise[self.noise_index];self.noise_index+=1;observed[:,0]+=error[:,0];observed[:,1]+=error[:,1]
   pressure,force,length,velocity=observed.T
   if self.initialized:
    self.pressure_rate+=self.dt/self.derivative_tau*((pressure-self.previous_pressure)/self.dt-self.pressure_rate)
    self.force_rate+=self.dt/self.derivative_tau*((force-self.previous_force)/self.dt-self.force_rate)
   else:self.initialized=True
   self.previous_pressure[:]=pressure;self.previous_force[:]=force
   factor=self.capacity_gain*p[:,3]/(620.*p[:,4]);hp=p[:,0]+p[:,2]+factor*pressure-length;hf=p[:,0]+p[:,2]+force/p[:,4]-length
   if self.mode=='pressure':reserve=hp;rate=factor*self.pressure_rate-velocity
   else:reserve=hf;rate=self.force_rate/p[:,4]-velocity
   predicted=reserve+(self.preview+self.delay_steps*self.dt)*np.minimum(rate,0.);signal=np.clip((self.reserve_target-predicted*p[:,4])/self.force_scale,0.,1.)
   old=self.relay_v.copy();self.increment[:]=self.pressure_cap*np.clip(old/10.,0.,1.);self.relay_v+=self.dt/self.tau*(10.*signal-self.relay_v)
   row[:]=np.column_stack((reserve,rate,predicted,signal,old,np.zeros(2),pressure,length,velocity,self.pressure_rate))
   sensor[:]=np.column_stack((force,self.force_rate,hp,hf,self.increment,error[:,0],error[:,1],np.full(2,self.delay_steps*self.dt)))
  output=self.inner.forward(x);self.logs.append(row);self.sensor_logs.append(sensor);self.steps+=1;return output
 def export_state(self):
  state=super().export_state();state.update(capacity_mode=np.array(self.mode),capacity_parameters=np.array([self.capacity_gain,self.delay_steps,self.pressure_noise,self.force_noise,self.noise_seed]),capacity_sensor_buffer=self.buffer.copy(),capacity_buffer_count=np.array(self.buffer_count),capacity_buffer_initialized=np.array(self.buffer_initialized),capacity_noise_index=np.array(self.noise_index),capacity_noise_table_hash=np.array(self.noise_hash),capacity_noise_policy=np.array('SeededPCG64 precomputed clippedGaussian; no runtimeRNG draws; full sequence reconstructed from parameters/hash'),capacity_previous_force=self.previous_force.copy(),capacity_force_rate=self.force_rate.copy(),capacity_sensor_logs=np.asarray(self.sensor_logs).reshape(-1,2,len(SENSOR_CHANNELS)))
  return state
 def import_state(self,state):
  assert str(state['capacity_mode'])==self.mode and np.array_equal(state['capacity_parameters'],[self.capacity_gain,self.delay_steps,self.pressure_noise,self.force_noise,self.noise_seed]) and str(state['capacity_noise_table_hash'])==self.noise_hash
  super().import_state(state);np.copyto(self.buffer,state['capacity_sensor_buffer']);np.copyto(self.previous_force,state['capacity_previous_force']);np.copyto(self.force_rate,state['capacity_force_rate']);self.buffer_count=int(state['capacity_buffer_count']);self.buffer_initialized=bool(state['capacity_buffer_initialized']);self.noise_index=int(state['capacity_noise_index']);self.sensor_logs=[row.copy() for row in state['capacity_sensor_logs']]
 def circuit(self):
  result=super().circuit();result.update(estimator_mode=self.mode,learned_pressure_capacity_input_gain=self.capacity_gain if self.mode=='pressure' else None,actual_actuator_Fmax_modified=False,oracle_case_gain_used=False,force_source='Existing raw positive native/Bank tension, same physical source as modelIb',additional_measurement_delay_s=self.delay_steps*self.dt,source_engine_sensor_latency_s=.002,pressure_noise_std_kpa=self.pressure_noise,force_noise_std_N=self.force_noise,noise_clip_sigma=3.,noise_seed=self.noise_seed,measurement_alignment='Pressure,tension,length,velocity delayed together; common noiseless pre-release buffer; seeded noise starts afterrelease')
  return result
