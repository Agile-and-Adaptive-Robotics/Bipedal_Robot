"""Boundary equations and native forward-domain checks before physics experiments."""
import json,ast,hashlib
from pathlib import Path
import numpy as np
from cad_interface import RoutedBPABank
from actuator_interface import force_at_pressure
from inverse_feasibility_policy import feasible_interval,commands,BoundaryBPABank
from normalized_prepared_transfer import ROOT,OUT

def main():
 config=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text());bank=RoutedBPABank(list(config['actuators']),config,.002);rows=[]
 for case in ('rot0p8_graded','rot0p8_spike','hold20_graded','hold20_spike'):
  r=json.loads((OUT/('w2lvar_'+case+'_v1.result.json')).read_text());d=ast.literal_eval(r['failure'].split('failing routes: ',1)[1])[0];i=bank.names.index(d['name']);p=bank.parameters[i];length=d['route_length_m'];lo,hi,tlo,thi=feasible_interval(length,p)
  assert np.isfinite(force_at_pressure(length,lo,p)) and np.isfinite(force_at_pressure(length,hi,p))
  assert np.isnan(force_at_pressure(length,d['actual_pressure_kpa'],p)), 'Historical forward guard must remain unchanged'
  np.testing.assert_allclose(tlo,p[4]*(length-p[2]-p[0]),atol=1e-8,rtol=0)
  for target in (0.,tlo*.99,tlo,(tlo+thi)*.5,thi,thi*1.01):
   out=np.zeros((1,9));chosen=np.zeros(1);flags=np.zeros(1,dtype=np.int64);invalid=np.zeros(1,dtype=np.bool_)
   commands(np.array([target/(thi*2)]),np.array([length]),p[None,:],np.array([thi*2]),out,chosen,flags,invalid)
   value=force_at_pressure(length,chosen[0],p);assert np.isfinite(value)
   np.testing.assert_allclose(value,np.clip(target,tlo,thi),atol=1e-8,rtol=1e-10)
  rows.append(dict(case=case,pressure_interval_kpa=[lo,hi],tension_interval_N=[tlo,thi],historical_low_pressure_still_invalid=True))
 # 86 real design routes, three length regimes each, preserving the zero-force branch.
 for p in bank.parameters:
  for length in (p[2]+p[0]*.7,p[2]+p[0]*.9,p[2]+p[0]+.00002):
   lo,hi,a,b=feasible_interval(length,p);assert np.isfinite(a) and np.isfinite(b) and 0<=a<=b
   for pressure in np.linspace(lo,hi,5):assert np.isfinite(force_at_pressure(length,pressure,p))
 p=bank.parameters[0];assert np.isnan(feasible_interval(p[2]+p[0]+p[3]*(p[8]/620)/p[4]+.001,p)[2])
 # Abort is atomic: invalid target must not alter pressure/history/force or counter.
 tiny=dict(config);tiny['actuators']={'pect_r':config['actuators']['pect_r']};b=BoundaryBPABank(['pect_r'],tiny,.002);b.policy='abort';initial=(b.pressure.copy(),b.history.copy(),b.force.copy(),b.step_count)
 try:b.step(np.array([0.]),np.array([b.parameters[0,0]+b.parameters[0,2]+.00005]));raise AssertionError('Expected strict abort')
 except ValueError as error:assert 'Explicit inverse feasibility abort' in str(error)
 assert np.array_equal(initial[0],b.pressure) and np.array_equal(initial[1],b.history) and np.array_equal(initial[2],b.force) and initial[3]==b.step_count
 # Actual86-route command/history/pressure/force parity when every request is feasible.
 old=RoutedBPABank(bank.names,config,.002);new=BoundaryBPABank(bank.names,config,.002);new.policy='nearest';lengths=old.parameters[:,2]+old.parameters[:,0]*.95
 limits=np.array([feasible_interval(length,p)[3] for length,p in zip(lengths,old.parameters)]);u=.5*np.minimum(limits,old.scales)/old.scales
 old.pressure[:]=100.;new.pressure[:]=100.
 for _ in range(20):
  old.step(u,lengths);new.step(u,lengths)
  for key in ('command','flags','history','pressure','force'):assert np.array_equal(getattr(old,key),getattr(new,key)),key
 assert not new.policy_events
 out=ROOT/'inverse_boundary_policy';out.mkdir(exist_ok=True);sources={str(Path(__file__).with_name(name)):hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest() for name in ('inverse_feasibility_policy.py','actuator_interface.py','cad_interface.py','native_bpa_force.py')}
 report=dict(machine='easteregg2',four_actual_failure_boundaries=rows,design_regime_checks=86*3,atomic_abort=True,feasible86route_20step_bit_exact_parity=True,forward_domain_unchanged=True,pressure_lag_not_projected=True,source_sha256=sources)
 (out/'BOUNDARY_UNIT_CHECKS.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
