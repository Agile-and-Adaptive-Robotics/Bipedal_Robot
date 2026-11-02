"""Read-only explanation of the unchanged frozen-transfer domain failures."""
import ast,json
from pathlib import Path
import numpy as np
from normalized_prepared_transfer import ROOT,OUT

def main():
 cfg=json.loads((ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json').read_text());rows=[]
 for case in ('rot0p8_graded','rot0p8_spike','hold20_graded','hold20_spike'):
  p=OUT/('w2lvar_'+case+'_v1');r=json.loads(p.with_suffix('.result.json').read_text());details=ast.literal_eval(r['failure'].split('failing routes: ',1)[1])
  with np.load(p.with_suffix('.pressure.npz')) as z,np.load(p.with_suffix('.motor_state.npz')) as motor:
   for d in details:
    name=d['name'];i=list(z['names']).index(name);aid=list(motor['actuator_names']).index(name);a=cfg['actuators'][name];last=z['samples'][-1,:,i];flags=z['samples'][:,3,i];minimum=a['series_stiffness_N_m']*(d['route_length_m']-d['zero_stretch_route_rest_length_m']);u=float(motor['excitation'][-1,aid])
    rows.append(dict(case=case,route=name,completed_free_s=r['completed_samples']*.002,declared_domain_failure=d,last_valid_sample_time_s=float(z['time'][-1]),last_valid_excitation=u,last_valid_requested_tension_N=u*a['force_command_scale_N'],last_valid_actual_tension_N=float(last[2]),last_valid_pressure_reference_kpa=float(last[0]),last_valid_actual_pressure_kpa=float(last[1]),last_valid_inverse_flag=int(last[3]),count_invalid_inverse_flag3=int(np.count_nonzero(flags==3)),minimum_tension_at_failure_geometry_N=minimum,scopes=['Last-valid requested tension and failure-geometry lower bound have adjacent rather than identical timestamps.','Flag3 specifically means requested tension leaves sleeve overextended; inverse currently returns0kPa and the bank consumes that value as a physical pressure command.','This is not a zero neural-output event. Inverse feasibility and pressure-domain failure do not establish a mechanical fall. No pressure or force floor was added.']))
 report=dict(machine='easteregg2',rows=rows,next_lever='Explicit invalid-inverse handling/actuator-feasibility feedback should be considered separately from changing anatomy or blindly changing muscle strength. Any future opt-in interface policy must retain the unchanged forward domain guard, declare pressure-lag/margin assumptions, and be compared with these frozen failures. No policy or geometry was changed in this audit.')
 (OUT/'INVERSE_FEASIBILITY_AUDIT.json').write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps([(x['case'],x['last_valid_requested_tension_N'],x['minimum_tension_at_failure_geometry_N']) for x in rows]))

if __name__=='__main__':main()
