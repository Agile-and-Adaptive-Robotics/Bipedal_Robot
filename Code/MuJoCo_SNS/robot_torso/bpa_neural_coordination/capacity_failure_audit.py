"""Post-hoc declared-profile authority audit; never used by neural controller."""
import json
from pathlib import Path
import numpy as np
from register_capacity_comparison import DEST
from cad_interface import RoutedBPABank

def main():
 cfg=json.loads((DEST/'PREREGISTRATION.json').read_text());retained=json.loads((DEST/'RETAINED_FAILURES.json').read_text());rows=[]
 for key,folder in retained.items():
  folder=Path(folder);summary=json.loads((folder/'CASE_SUMMARY.json').read_text());trial=json.loads((folder.parent/'TRIAL_SUMMARY.json').read_text());params=trial['diagnostics']['params'];spec=next(f for f in cfg['training'] if f['name']==folder.name);actual=json.loads(Path(spec['path']).read_text());bank=RoutedBPABank(list(actual['actuators']),actual,.002)
  with np.load(folder/'w2lvar_case.relay.npz') as z:relay=z['samples'];routes=z['routes']
  with np.load(folder/'w2lvar_case.valve.npz') as z:valve=z['samples']
  with np.load(folder/'ESTIMATOR_SENSORS.npz') as z:sensors=z['samples']
  for j,name in enumerate(routes):
   name=str(name);p=bank.parameters[bank.names.index(name)];gain=spec['gain'];length=relay[:,j,7];velocity=relay[:,j,8];pressure=relay[:,j,6]-sensors[:,j,5];valid=np.flatnonzero(length>0)
   if not len(valid):continue
   i=int(valid[-1]);required=620*p[4]*max(0.,length[i]-p[0]-p[2])/(p[3]*gain);required_rate=620*p[4]*velocity[i]/(p[3]*gain);fill=actual['actuators'][name]['tau_fill_s'];reference_required=required+fill*max(0.,required_rate);tau=4;aligned=max(0,i-tau);aligned=min(aligned,len(valve)-1);nominal=valve[aligned,j,0];available=nominal+params['pressure_cap_kpa']
   rows.append(dict(mode=params['estimator_mode'],condition=spec['name'],route=name,physical_failure=summary['failure'],observation_neural_step=i,delayed_actual_pressure_kpa=float(pressure[i]),instantaneous_required_pressure_kpa=float(required),required_pressure_ramp_kpa_s=float(required_rate),fill_tau_s=fill,optimistic_ramp_tracking_reference_kpa=float(reference_required),nominal_reference_aligned_kpa=float(nominal),available_reference_at_cap_kpa=float(available),optimistic_reference_authority_margin_kpa=float(available-reference_required),neural_requested_increment_now_kpa=float(sensors[i,j,4]),neural_voltage_now_mV=float(relay[i,j,4]),input_current_nA=float(10*relay[i,j,3]),scope='Post-hoc uses known experimental mapgain and true filltau only for diagnosis, never controller input. Alignment assumes common6ms transport+2ms native sensor history (four2ms samples). Positive routevelocity predicts rising boundary; first-order ramp criterion ignores actuator delay and body/controller interaction, so negative margin is a local diagnostic, not a closed-loop impossibility proof.'))
 result=dict(rows=rows,selection_used=False,controller_oracle_used=False,source='Retained full unique failures only; no new simulations')
 (DEST/'POSTHOC_AUTHORITY_AUDIT.json').write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
