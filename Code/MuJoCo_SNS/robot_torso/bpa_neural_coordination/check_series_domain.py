"""Regression: domain applies to the sleeve after force-dependent series stretch."""
import os,sys,json
from pathlib import Path
os.environ['CONDA_PREFIX']=str(Path(sys.executable).parent)
import numpy as np
from actuator_interface import packed_parameters,force_at_pressure,inverse_pressure
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from bpa_muscle import BPAMuscle

def main():
 cases=[];invalid=0
 for dia in (10,20,40):
  for count in (1,2):
   p=packed_parameters(dia,.3,.25,.06,count*2e5);p[3]*=count
   ref=BPAMuscle('series_domain',dia,.3,.225,fitting_length=.03,bpa_count=count,series_stiffness=count*2e5)
   for pressure in (100.,300.,620.):
    ref.set_activation(pressure/620.)
    maximum_stretch=p[3]*(pressure/620.)/p[4]
    for fraction in (.5,.95):
     length=.36+fraction*maximum_stretch
     force=float(force_at_pressure(length,pressure,p));expected=ref.force(length)
     actual_sleeve_length=length-.06-force/p[4]
     assert np.isfinite(force) and force>0 and .225-1e-12<=actual_sleeve_length<=.3+1e-12
     inv,flag=inverse_pressure(force,length,p)
     assert flag==0 and abs(inv-pressure)<1e-7
     assert abs(force-expected)<1e-6
     cases.append(dict(diameter_mm=dia,count=count,pressure_kpa=pressure,total_route_length_m=length,actual_sleeve_length_m=actual_sleeve_length,series_stretch_m=force/p[4],tension_N=force,reference_error_N=abs(force-expected),inverse_error_kpa=abs(inv-pressure)))
    assert np.isnan(force_at_pressure(.36+1.1*maximum_stretch,pressure,p));invalid+=1
   assert np.isnan(force_at_pressure(.3601,0.,p));invalid+=1
 # Forward-map error changes the series equilibrium itself, while the
 # controller's inverse map remains nominal. Scaling solved force afterward
 # would give the wrong cable stretch and is explicitly excluded.
 uncertain=[]
 for gain in (.8,1.2):
  p=packed_parameters(20,.3,.25,.06,4e5);p[3]*=2*gain
  ref=BPAMuscle('uncertain_series',20,.3,.225,fitting_length=.03,bpa_count=2*gain,series_stiffness=4e5);ref.set_activation(.5)
  length=.3604;actual=float(force_at_pressure(length,310.,p));expected=ref.force(length)
  assert np.isfinite(actual) and abs(actual-expected)<1e-6
  assert length-.06-actual/p[4]<=.3
  uncertain.append(dict(forward_force_gain=gain,reference_error_N=abs(actual-expected),tension_N=actual))
 result=dict(machine='easteregg2',passed=True,valid_series_extended_cases=len(cases),correctly_rejected_insufficient_pressure_cases=invalid,max_reference_error_N=max(c['reference_error_N'] for c in cases),max_inverse_error_kpa=max(c['inverse_error_kpa'] for c in cases),cases=cases,uncertain_forward_equilibrium_cases=uncertain,
 scope='No sleeve extrapolation: valid equilibrium requires Lsleeve=Lroute−Lfixed−T/K within the characterized contraction range; cable stretch alone may make total route longer than sleeve rest+fixed. Reference is existing BPAMuscle scalar equilibrium on valid sleeve states.')
 out=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa\series_domain_regression.json');out.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2))

if __name__=='__main__':main()
