"""Small independent neuronal/parity/state checks; no walker trials."""
import json
from types import SimpleNamespace
import numpy as np
from pathlib import Path
from capacity_estimator_neural import CapacityReserveForward
from pressure_valve_neural import ValveReserveForward
from register_capacity_comparison import DEST

class Inner:
 def __init__(self):self.c=SimpleNamespace(input_connectivity=np.eye(2),V=np.zeros(2),V_last=np.zeros(2))
 def forward(self,x):self.c.V_last[:]=self.c.V;self.c.V[:]=.9*self.c.V+.1*x;return self.c.V.copy()

def context():
 p=np.array([[.15,.05,.025,2000.,400000.,0.,0.,0.,620.]]*2)
 bank=SimpleNamespace(names=['pect_r','pect_l'],parameters=p,pressure=np.array([60.,65.]),force=np.array([10.,12.]))
 return dict(active=True,adapter=SimpleNamespace(bank=bank,ids=np.array([0,1])),data=SimpleNamespace(actuator_length=np.array([.175,.175]),actuator_velocity=np.array([0.,0.])))

def main():
 result={};names=['MN_pect_r','MN_pect_l'];x=np.array([.2,.4])
 for cap in (0.,1.):
  a=ValveReserveForward(Inner(),names,context(),cap,.1,1.);b=CapacityReserveForward(Inner(),names,context(),cap,.1,1.)
  for i in range(80):
   for w in (a,b):w.context['adapter'].bank.pressure[:]=[60.-.02*i,65.-.01*i]
   assert np.array_equal(a.forward(x),b.forward(x)) and np.array_equal(a.increment,b.increment) and np.array_equal(a.logs[-1],b.logs[-1])
  result['zero_output_exact' if cap==0 else 'default_relay_exact']=True
 for mode in ('pressure','force'):
  w=CapacityReserveForward(Inner(),names,context(),5.,.18,1.,mode=mode,capacity_gain=.8 if mode=='pressure' else 1.,delay_steps=3,pressure_noise_kpa=.02,force_noise_N=.05,noise_seed=511031)
  for i in range(80):w.forward(x)
  state=w.export_state()
  for i in range(32):w.forward(x)
  end=w.export_state();w.import_state(state)
  for i in range(32):w.forward(x)
  repeated=w.export_state();checks={key:bool(np.array_equal(end[key],repeated[key])) for key in end};assert all(checks.values());result[mode+'_all_state_repeat']=checks
 from sns_toolbox.networks import Network
 from sns_toolbox.neurons import NonSpikingNeuron
 net=Network(name='Independent capacity reserve relays')
 for n in ('R','L'):
  net.add_neuron(NonSpikingNeuron(membrane_capacitance=6.,membrane_conductance=1.),name=n);net.add_input(dest=n);net.add_output(source=n)
 model=net.compile(dt=2.,backend='numpy');v=np.zeros(2);maximum=0.
 for i in range(200):
  signal=np.array([.5+.4*np.sin(i*.17),.5+.4*np.cos(i*.23)]);out=model.forward(10*signal);v+=2./6.*(10*signal-v);maximum=max(maximum,float(np.max(np.abs(v-model.V))))
 assert maximum==0.;result['independent_SNS_voltage_max_error_mV']=maximum
 (DEST/'INDEPENDENT_UNIT_STATE_CHECKS.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(dict(zero_exact=True,default_exact=True,pressure_state_exact=True,force_state_exact=True,SNS_error=maximum)))

def physical_zero():
 from capacity_comparison_adapter import registration,validator
 import reserve_valve_runner as runner
 from compact_trial_traces import CompactTraces
 from register_capacity_comparison import ROOT
 cfg=registration();spec=cfg['training'][0];params=dict(estimator_mode='force',pressure_cap_kpa=0.,preview_s=.1);folder=DEST/'zero_physical_parity';folder.mkdir(exist_ok=True);original=runner.ValveReserveForward
 def factory(inner,names,context,cap,preview,target):return CapacityReserveForward(inner,names,context,cap,preview,target,mode='force',delay_steps=3,pressure_noise_kpa=.02,force_noise_N=.05,noise_seed=spec['sensor_noise_seed'])
 runner.ValveReserveForward=factory
 try:
  with CompactTraces(folder) as traces:
   c=runner.run(cap=0.,preview=.1,condition=spec['name'],encoding='graded',tag='case',duration=20.,reserve_target=1.,output_dir=folder,actuator_configuration=spec['path'],preparation_configuration=cfg['common_preparation']['path'],common_release_validator=validator(params,spec['sensor_noise_seed']),sensor_gain=1.,mechanical_push_bw=0.);checks={};old=ROOT/'pressure_reserve_valve/w2lvar_zero_current'
   for suffix in ('.motor_state.npz','.pressure.npz','.inputs.npz','.apparatus.npz','.apparatus_extended.npz','.relay.npz','.valve.npz'):
    with np.load(old.with_suffix(suffix)) as a,np.load((folder/'w2lvar_case').with_suffix(suffix)) as b:
     assert set(a.files)==set(b.files);checks[suffix]={key:bool(np.array_equal(a[key],b[key])) for key in a.files}
   assert all(v for group in checks.values() for v in group.values())
  traces.materialize(full=False);result=dict(exact=True,all_arrays=checks,baseline_accepted=c['assessment']['accepted'],scope='Output-off physical control exactly reproduces preserved domain-failing nearest-reference baseline, not an accepted standing solution.')
  (DEST/'ZERO_PHYSICAL_PARITY.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
 finally:runner.ValveReserveForward=original

if __name__=='__main__':
 import sys
 physical_zero() if '--physical-zero' in sys.argv else main()
