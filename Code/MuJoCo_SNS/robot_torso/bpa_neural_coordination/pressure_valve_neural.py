"""Explicit reserve-neuron voltage → valve pressure-reference increment variant."""
import numpy as np
from pressure_reserve_neural import ReserveForward,CHANNELS
from inverse_feasibility_policy import BoundaryBPABank,commands,advance
from cad_interface import RoutedBPABank

class ValveReserveForward(ReserveForward):
 def __init__(self,inner,names,context,pressure_cap_kpa,preview_s,reserve_target_N=0.):
  super().__init__(inner,names,context,0.,preview_s)
  if not 0<=pressure_cap_kpa<=5.:raise ValueError('Declared valve increment bound is0..5kPa')
  if not 0<=reserve_target_N<=1.:raise ValueError('Declared reserve target is0..1N equivalent series reserve')
  self.reserve_target=float(reserve_target_N);self.pressure_cap=float(pressure_cap_kpa);self.cap=float(self.pressure_cap>0);self.increment=np.zeros(2);context['valve_increment']=self.increment
 def _step(self,x):
  row=np.zeros((2,len(CHANNELS)));np.copyto(self.input_buffer,x);self.increment[:]=0.
  if self.context.get('active',False) and self.pressure_cap>0:
   adapter=self.context['adapter'];data=self.context['data'];idx=np.array([adapter.bank.names.index(a) for a in self.routes]);ids=adapter.ids[idx];p=adapter.bank.parameters[idx];pressure=adapter.bank.pressure[idx];length=data.actuator_length[ids];velocity=data.actuator_velocity[ids]
   if self.initialized:raw=(pressure-self.previous_pressure)/self.dt;self.pressure_rate+=self.dt/self.derivative_tau*(raw-self.pressure_rate)
   else:self.initialized=True
   self.previous_pressure[:]=pressure;factor=p[:,3]/(620.*p[:,4]);h=p[:,0]+p[:,2]+factor*pressure-length;rate=factor*self.pressure_rate-velocity;predicted=h+self.preview*np.minimum(rate,0.);signal=np.clip((self.reserve_target-predicted*p[:,4])/self.force_scale,0.,1.)
   old=self.relay_v.copy();self.increment[:]=self.pressure_cap*np.clip(old/10.,0.,1.);self.relay_v+=self.dt/self.tau*(10.*signal-self.relay_v)
   row[:]=np.column_stack((h,rate,predicted,signal,old,np.zeros(2),pressure,length,velocity,self.pressure_rate))
  # Original neural input, motor voltage, and motor output are untouched.
  result=self.inner.forward(x);self.logs.append(row);self.steps+=1;return result
 def export_state(self):
  state=super().export_state();state['valve_pressure_cap_kpa']=np.array(self.pressure_cap);state['valve_reserve_target_N']=np.array(self.reserve_target);state['valve_increment_kpa']=self.increment.copy();state['valve_output_mode']=np.array('explicit_neuron_voltage_to_reference_pressure')
  return state
 def import_state(self,state):
  assert float(state['valve_pressure_cap_kpa'])==self.pressure_cap and float(state['valve_reserve_target_N'])==self.reserve_target
  super().import_state(state);np.copyto(self.increment,state['valve_increment_kpa'])
 def circuit(self):
  return dict(architecture='Two reserve leaky neurons feed explicit voltage-to-pressure increment converters; original network supplies nominal reference. No added existing-MN current.',neurons=[dict(name='RESERVE_'+a,tau_ms=6.,leak_uS=1.,capacitance_nF=6.,sensor_current_cap_nA=10.) for a in self.routes],output_maps=[dict(source='RESERVE_'+a,voltage_window_mV=[0.,10.],target='VALVE_PRESSURE_REFERENCE_'+a,increment_kpa=[0.,self.pressure_cap]) for a in self.routes],preview_s=self.preview,pressure_derivative_filter_s=self.derivative_tau,sensor_force_scale_N=self.force_scale,desired_reserve_force_equivalent_N=self.reserve_target,desired_reserve_length_m=self.reserve_target/400000.,bilateral_tying=True,actual_pressure_projection=False,force_floor=False,nominal_force_map_only=True)

class ValveBank(BoundaryBPABank):
 def initialize_valve(self,context):
  self.valve_context=context;self.valve_ids=np.array([self.names.index(a) for a in ('pect_r','pect_l')]);self.valve_logs=[];self.valve_steps=0
 def step(self,u,lengths):
  if self.policy is None:return RoutedBPABank.step(self,u,lengths)
  u=np.asarray(u);lengths=np.asarray(lengths)
  if u.shape!=self.pressure.shape or lengths.shape!=u.shape or not np.isfinite(u).all() or not np.isfinite(lengths).all():raise ValueError('Invalid valve-policy input')
  commands(u,lengths,self.parameters,self.scales,self.last_policy,self.selected,self.flags,self.invalid)
  for i in np.flatnonzero(self.invalid):self.policy_events.append(np.r_[float(self.policy_steps),float(i),self.last_policy[i]].copy())
  if not np.isfinite(self.selected).all():raise ValueError('No feasible nominal static pressure interval')
  nominal=self.selected[self.valve_ids].copy();requested=self.valve_context.get('valve_increment',np.zeros(2))
  if not np.isfinite(requested).all() or np.any(requested<0):raise ValueError('Invalid neural valve increment')
  self.selected[self.valve_ids]=np.minimum(self.parameters[self.valve_ids,8],nominal+requested)
  applied=self.selected[self.valve_ids]-nominal
  self.valve_logs.append(np.column_stack((nominal,requested,applied,self.selected[self.valve_ids])).copy())
  self.command[:]=self.selected
  advance(lengths,self.parameters,self.force_parameters,self.selected,self.history,self.step_count,self.delay,self.pressure,self.fill,self.vent,self.force,self.dt)
  self.step_count+=1;self.policy_steps+=1;self.valve_steps+=1
  if not np.isfinite(self.force).all():raise ValueError('BPA length is outside the declared contraction domain')
  return self.force

def install_valve_adapter(base,module,context):
 class Adapter(base):
  def __init__(self,*args,**kwargs):
   super().__init__(*args,**kwargs);self.bank.__class__=ValveBank;self.bank.initialize_valve(context);context['adapter']=self
  def advance(self,data,u):context['data']=data;return super().advance(data,u)
  def configure_response(self,cfg):super().configure_response(cfg);context['active']=True
 module.NativeBPAForce=Adapter
