"""Declared virtual BPA sizing, not identification of mounted CAD hardware.

Both sides share each design. Route-length bounds come from the completed
standing trajectory and do not establish a walking-range or packaging design.
"""
import json,hashlib
from pathlib import Path
import gait_standing_family as F
import numpy as np,mujoco
from actuator_interface import packed_parameters,force_at_pressure

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main():
 checkpoint=json.loads((ROOT/'STANDING_CHECKPOINT.json').read_text())
 if not checkpoint['accepted_standing_candidate']:raise ValueError('Standing prerequisite not satisfied')
 source=Path(checkpoint['fresh_result']);z=np.load(source.with_suffix('').with_suffix('.npz'))
 base=mujoco.MjModel.from_xml_path(str(F.R.MODEL));data=mujoco.MjData(base)
 F.R.apply_start_pose(base,data);m=F.R.apply_harness(base,data);d=mujoco.MjData(m)
 forcepath=F.SPINAL/'physics_fidelity/shared_force_law_functional_v3.json'
 for row in json.loads(forcepath.read_text())['muscles'].values():
  for name in row['actuator_names']:
   for key in ('gainprm','biasprm','dynprm','lengthrange'):getattr(m,'actuator_'+key)[m.actuator(name).id]=row[key]
 names=[m.actuator(i).name for i in range(m.nu)];q=z['qfull']
 if q.shape[1]!=m.nq:raise ValueError('No reinterpretation of differing qpos layouts')
 lengths=[]
 for index in np.linspace(0,len(q)-1,161,dtype=int):
  d.qpos[:]=q[index];mujoco.mj_forward(m,d);lengths.append(d.actuator_length.copy())
 lengths=np.array(lengths);lookup={n:i for i,n in enumerate(names)};rows=[];actuators={};done=set()
 fixed=.05;contraction=.25;count=2;diameter=20
 for name in names:
  if name in done:continue
  partner=F.mirror(name);pair=[name,partner] if partner!=name else [name]
  ids=[lookup[n] for n in pair];done.update(pair)
  if all(n in F.R.PRUNE_MUSCLES for n in pair):continue
  minimum=float(lengths[:,ids].min());maximum=float(lengths[:,ids].max())
  lower=maximum-fixed;upper=(minimum-fixed)/(1-contraction)
  feasible=bool(lower>0 and upper>=lower and upper>.0075)
  row=dict(names=pair,route_min_m=minimum,route_max_m=maximum,assumed_total_fixed_fitting_length_m=fixed,rest_lower_bound_m=lower,rest_upper_bound_m=upper,feasible=feasible)
  if feasible:
   rest=lower+.15*(upper-lower)
   scale=float(np.mean(m.actuator_gainprm[ids,2]))
   spec=dict(diameter_mm=diameter,resting_length_m=rest,maximum_contraction_fraction=contraction,
    nonmuscle_length_m=fixed,series_stiffness_N_m=count*2e5,parallel_count=count,
    pressure_max_kpa=620.,tau_fill_s=.04,tau_vent_s=.06,delay_s=.004,
    force_command_scale_N=scale,
    provenance='VIRTUAL DESIGN, not mounted hardware identification: two20mm BPAs per route, total50mm fixed fittings,25%stroke,200kN/m series stiffness per parallel unit,40/60ms pressure constants,4ms delay are explicit assumptions. Rest length is derived jointly for both sides from this standing-only route envelope. Neural force-reference scaling is pooled biological Fmax; it is not a BPA capacity claim.')
   p=packed_parameters(diameter,rest,contraction,fixed,spec['series_stiffness_N_m']);p[3]*=count
   capacity=[float(force_at_pressure(l,620.,p)) for l in (minimum,maximum)]
   row.update(resting_length_m=rest,reference_scale_N=scale,maximum_pressure_tension_at_envelope_ends_N=capacity)
   for n in pair:actuators[n]=spec.copy()
  rows.append(row)
 disabled=sorted(F.R.PRUNE_MUSCLES&set(names));missing=sorted(set(names)-set(actuators)-set(disabled))
 out=ROOT/'virtual_bpa_design';out.mkdir(exist_ok=True)
 configuration=dict(command_mode='inverse_force',actuators=actuators,disabled_actuators=disabled,
  valid_complete_assignment=not missing,unconfigured_actuators=missing,
  scope='Virtual Gait actuator comparison only; no mounted CAD identification, packaging validation, or full gait-ROM feasibility. Pressure dynamics are hypothetical stress-test parameters.',initial_pressure_state='Explicit precharge from the shared biological initializer activation after the common0.6s pinned warmup.')
 (out/'VIRTUAL_GAIT_BPA_CONFIGURATION.json').write_text(json.dumps(configuration,indent=2),encoding='utf-8')
 report=dict(machine='easteregg2',source_result=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),poses=len(lengths),model=str(F.R.MODEL),model_sha256=hashlib.sha256(F.R.MODEL.read_bytes()).hexdigest(),force_law=str(forcepath),configured=len(actuators),disabled=len(disabled),unconfigured=missing,valid_complete_assignment=not missing,designs=rows,scope=configuration['scope'])
 (out/'DESIGN_FEASIBILITY.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
 print(json.dumps({k:v for k,v in report.items() if k!='designs'},indent=2))

if __name__=='__main__':main()
