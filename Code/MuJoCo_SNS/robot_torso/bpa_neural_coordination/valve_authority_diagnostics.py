"""Read-only physical/neural attribution from full or compact failure witnesses."""
import json,ast
from pathlib import Path
import numpy as np
from register_valve_authority import DEST

def load(folder,name):
 p=folder/name
 if p.exists():
  with np.load(p) as z:return {k:z[k] for k in z.files}
 manifest=json.loads((folder/'TELEMETRY_MANIFEST.json').read_text())
 for i,item in enumerate(manifest['archives']):
  if Path(item['original_path']).name==name:
   with np.load(folder/'COMPACT_TELEMETRY.npz') as z:return {key:z[f'a{i}_{key}'] for key in item['arrays'] if f'a{i}_{key}' in z.files}
 raise FileNotFoundError(name)

def main():
 cfg=json.loads((DEST/'PREREGISTRATION.json').read_text());rows=[]
 for line in (DEST/'pressure/rustuna/TRIALS.jsonl').read_text().splitlines():
  t=json.loads(line)
  for c in t['diagnostics']['conditions']:
   if not c['failure'] or 'failing routes: ' not in c['failure']:continue
   folder=Path(c['result']).parent;failure=load(folder,'w2lvar_case.final_actuator_mechanics_state.npz');motors=load(folder,'w2lvar_case.motor_state.npz');relay=load(folder,'w2lvar_case.relay.npz');valve=load(folder,'w2lvar_case.valve.npz');spec=next(s for s in cfg['training'] if s['name']==c['condition']);actual=json.loads(Path(spec['path']).read_text())
   for d in ast.literal_eval(c['failure'].split('failing routes: ',1)[1]):
    name=d['name'];p=actual['actuators'][name];i=motors['actuator_names'].tolist().index(name);side=relay['routes'].tolist().index(name);excitation=float(failure['ctrl'][i]);target=excitation*p['force_command_scale_N'];stiffness=p['series_stiffness_N_m'];minimum=max(0.,stiffness*(d['route_length_m']-d['zero_stretch_route_rest_length_m']));static_extra=(1-spec['gain'])*d['actual_pressure_kpa'];rows.append(dict(trial=t['number'],cap_kpa=t['params']['pressure_cap_kpa'],condition=c['condition'],route=name,moving_s=c['moving_s'],actual_pressure_kpa=d['actual_pressure_kpa'],actual_length_m=d['route_length_m'],minimum_geometrically_feasible_tension_N=minimum,neural_excitation=excitation,original_neural_requested_tension_N=target,requested_minus_minimum_N=target-minimum,last_neural_voltage_mV=float(relay['samples'][-1,side,4]),last_added_reference_kpa=float(valve['samples'][-1,side,2]),boundary_static_map_correction_kpa=static_extra,authority_remaining_beyond_static_map_kpa=t['params']['pressure_cap_kpa']-static_extra,scope='Exact failure-state excitation and route length, not legacyactivation. MinimumT uses virtualseries stiffness. Static correction uses declaredtrue mapgain only offline, never a controller input; excludes pressurelag/delay and closedloopbodymotion.'))
 report=dict(rows=rows,source='Existing full/compact witnesses; no additional simulations',no_oracle_supplied_to_controller=True,interpretation='At these stretched routes the original neural tension request can be below the geometrically feasible BPA minimum. Extra valve authority maintains the sleeve domain temporarily but does not by itself correct the underlying posture/recruitment mismatch. A finite fixed pressure offset also spends increasing authority on the static map deficit as route stretch grows. No extrapolation or domain relaxation.')
 (DEST/'FAILURE_ATTRIBUTION.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
