"""Opt-in static inverse boundary handling; no forward-force/domain alteration.

At route L and series stiffness K: sleeve=L-fixed-T/K. Nonnegative
contraction requires T>=K*(L-fixed-rest). At zero contraction F=Fmax*P/620,
so the minimum feasible pressure is 620*max(0,K*(L-fixed-rest))/Fmax.
Forward equilibrium at that pressure and Pmax gives the attainable tension
interval. Shorter-than-full-stroke routes retain the original zero-force branch.
Bounds use the NOMINAL inverse map, not oracle knowledge of held-out force errors.
Static command feasibility does not guarantee dynamic feasibility with pressure lag.
"""
import numpy as np
from numba import njit
from cad_interface import RoutedBPABank
from actuator_interface import force_at_pressure,inverse_pressure,pressure_step

CHANNELS=['requested_tension_N','minimum_tension_N','maximum_tension_N','minimum_pressure_kpa','maximum_pressure_kpa','selected_static_tension_N','selected_pressure_kpa','static_tension_residual_N','raw_inverse_flag']

@njit(cache=True)
def feasible_interval(length,p):
 low=620.*max(0.,p[4]*(length-p[2]-p[0]))/p[3];high=p[8]
 if low>high:return low,high,np.nan,np.nan
 # Only directed floating-point rounding, never a domain tolerance or pressure margin.
 for _ in range(8):
  tlo=force_at_pressure(length,low,p)
  if np.isfinite(tlo):break
  low=np.nextafter(low,np.inf)
 thi=force_at_pressure(length,high,p)
 if low>high or not np.isfinite(tlo) or not np.isfinite(thi):return low,high,np.nan,np.nan
 return low,high,tlo,thi

@njit(cache=True)
def commands(u,lengths,parameters,scales,out,selected,flags,invalid):
 for i in range(len(u)):
  target=min(1.,max(0.,u[i]))*scales[i];p=parameters[i]
  raw,flag=inverse_pressure(target,lengths[i],p);plo,phi,tlo,thi=feasible_interval(lengths[i],p)
  flags[i]=flag;bad=not np.isfinite(tlo) or target<tlo or target>thi or flag!=0
  invalid[i]=bad
  if not np.isfinite(tlo):chosen=np.nan;tref=np.nan
  elif target<=tlo:chosen=plo;tref=tlo
  elif target>=thi:chosen=phi;tref=thi
  else:
   chosen=min(phi,max(plo,raw));tref=force_at_pressure(lengths[i],chosen,p)
  selected[i]=chosen
  out[i,0]=target;out[i,1]=tlo;out[i,2]=thi;out[i,3]=plo;out[i,4]=phi;out[i,5]=tref;out[i,6]=chosen;out[i,7]=tref-target;out[i,8]=flag

@njit(cache=True)
def advance(lengths,parameters,force_parameters,selected,history,count,delay,pressure,fill,vent,force,dt):
 for i in range(len(lengths)):
  history[i,count%history.shape[1]]=selected[i]
  command=history[i,(count-delay[i])%history.shape[1]] if count>=delay[i] else 0.
  pressure[i]=pressure_step(pressure[i],command,dt,fill[i],vent[i])
  force[i]=force_at_pressure(lengths[i],pressure[i],force_parameters[i])

class BoundaryBPABank(RoutedBPABank):
 def __init__(self,*args,**kwargs):
  super().__init__(*args,**kwargs);self.policy=None;self.policy_steps=0;self.policy_events=[]
  self.last_policy=np.zeros((len(self.names),len(CHANNELS)));self.selected=np.zeros(len(self.names));self.invalid=np.zeros(len(self.names),dtype=np.bool_)
 def step(self,u,lengths):
  if self.policy is None:return super().step(u,lengths)
  if self.mode!='inverse_force':raise ValueError('Boundary policy requires explicit inverse-force mode')
  u=np.asarray(u);lengths=np.asarray(lengths)
  if u.shape!=self.pressure.shape or lengths.shape!=u.shape or not np.isfinite(u).all() or not np.isfinite(lengths).all():raise ValueError('Invalid policy input')
  commands(u,lengths,self.parameters,self.scales,self.last_policy,self.selected,self.flags,self.invalid)
  for i in np.flatnonzero(self.invalid):self.policy_events.append(np.r_[float(self.policy_steps),float(i),self.last_policy[i]].copy())
  if self.policy=='abort' and np.any(self.invalid):
   i=int(np.flatnonzero(self.invalid)[0]);raise ValueError('Explicit inverse feasibility abort: '+self.names[i]+' '+str(dict(zip(CHANNELS,self.last_policy[i].tolist()))))
  if not np.isfinite(self.selected).all():raise ValueError('No feasible static pressure interval within declared pressure limit')
  self.command[:]=self.selected
  advance(lengths,self.parameters,self.force_parameters,self.selected,self.history,self.step_count,self.delay,self.pressure,self.fill,self.vent,self.force,self.dt)
  self.step_count+=1;self.policy_steps+=1
  if not np.isfinite(self.force).all():raise ValueError('BPA length is outside the declared contraction domain')
  return self.force

def install(policy):
 """Process-local adapter subclass; preparation stays identical until response switch."""
 if policy not in ('abort','nearest'):raise ValueError(policy)
 import native_bpa_force as module
 original=module.NativeBPAForce;instances=[]
 class Adapter(original):
  def __init__(self,*args,**kwargs):
   bank_original=module.RoutedBPABank;module.RoutedBPABank=BoundaryBPABank
   try:super().__init__(*args,**kwargs)
   finally:module.RoutedBPABank=bank_original
   instances.append(self)
  def configure_response(self,configuration):
   super().configure_response(configuration);self.bank.policy=policy
 module.NativeBPAForce=Adapter
 return module,original,instances
