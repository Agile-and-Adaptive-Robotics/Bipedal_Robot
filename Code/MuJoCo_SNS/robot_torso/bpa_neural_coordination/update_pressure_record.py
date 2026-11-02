"""Summarize existing pressure experiments without rerunning or manuscript edits."""
import json
from pathlib import Path
import numpy as np
from gait_standing_objective import score

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def compare(a,b,keys):
 za=np.load(a);zb=np.load(b);result={}
 for key in keys:
  x,y=za[key],zb[key]
  if x.shape!=y.shape:result[key]=dict(shape_equal=False);continue
  result[key]=dict(shape_equal=True,bit_exact=np.array_equal(x,y),max_abs_difference=float(np.max(abs(x-y))) if x.size else 0.)
 return result

def main():
 source=ROOT/'pressure_standing/syn6_nominal_v1'
 fresh=ROOT/'pressure_standing/syn6_nominal_fresh_v1'
 r=json.loads(fresh.with_suffix('.result.json').read_text())
 checks=compare(source.with_suffix('.npz'),fresh.with_suffix('.npz'),('t','q','qfull','ctrl','neuro'))
 checks.update({'apparatus_'+k:v for k,v in compare(source.with_suffix('.apparatus.npz'),fresh.with_suffix('.apparatus.npz'),('samples','neuro')).items()})
 pa=np.load(source.with_suffix('.pressure.npz'));pb=np.load(fresh.with_suffix('.pressure.npz'))
 pressure_exact=np.array_equal(pa['samples'],pb['samples'][:,:5,:])
 _,assessment=score(r,fresh)
 longpath=ROOT/'pressure_standing/syn6_nominal20s_exact_v1.result.json'
 long=json.loads(longpath.read_text()) if longpath.exists() else None
 longer_assessment=score(long,longpath.with_suffix('').with_suffix(''))[1] if long else None
 basin=json.loads((ROOT/'standing_basin/BASIN_SUMMARY.json').read_text())
 heldpath=ROOT/'pressure_heldouts/HELDOUT_SUMMARY.json';held=json.loads(heldpath.read_text()) if heldpath.exists() else None
 report=dict(machine='easteregg2',condition='Virtual pressure-actuated Gait standing',controller_configuration=str(ROOT/'rustuna_syn6_v1/BEST.json'),actuator_configuration=str(ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json'),fresh_result=str(fresh.with_suffix('.result.json')),fresh_assessment=assessment,
 fresh_comparisons=checks,pressure_first_five_channels_bit_exact=pressure_exact,fresh_process_bit_exact=pressure_exact and all(v.get('bit_exact',False) for v in checks.values()),
 pressure_peak_kpa=r['pressure_peak_kpa'],inverse_saturation_fraction=r['inverse_saturation_fraction'],native_force_error_N=r['native_pressure_force_max_error_N'],
 long_result=str(longpath) if long else None,long_assessment=longer_assessment,biological_standing_basin=dict(accepted=basin['accepted'],completed=basin['completed']),heldout_summary=str(heldpath),
 limitations=['86virtual BPA pairs with assumed sizing/dynamics,6previously pruned rotators disabled; not mounted CAD hardware.','Controller is a separately fitted syn6 standing mode, not the locked Hunt controller or a curriculum pass.','Pelvis translation/rotation and lumbar/hip-rotation apparatus springs remain.','Native study duration includes0.6s pinned initialization; free-physics/pressure clocks exclude it.','Combined slow/weak/direction-memory stressor exceeded declared contraction domain before recovery; not evidence of general robustness.'])
 checkpoint=ROOT/'BPA_STANDING_CHECKPOINT.json';checkpoint.write_text(json.dumps(report,indent=2),encoding='utf-8')
 lines=['## Virtual pressure-actuated Gait comparison','','The biological syn6 standing candidate passed all8of8±1%corners of its three shared fitted parameters. Its exact-JSON fresh repeat is bit-identical. This does not change the failed standing status of the separate locked Hunt controller.','','A separate virtual design assigns two20mm BPAs per active route, both sides sharing dimensions and force-reference scale. Rest lengths are derived from the observed standing-only route envelope;50mm fixed fittings,25%stroke,400kN/m total series stiffness,40/60ms fill/vent constants, and4ms delay are declared assumptions. There are86active configured routes; the existing6pruned rotators remain inactive. This is **not identified mounted CAD hardware**, and full gait range/packaging are not established.','','The native user-force callback returns actual negative tensile force, which is visible to native force sensors and SNS feedback. It never cancels MN synaptic currents or applies an external torque allocator. Pressure updates once per control tick; repeated MuJoCo forward calls do not advance it. The native force/power regression passed6opposing-route/direction cases.','','| Condition | Accepted standing envelope | Tilt | Sway | Upward/downward rig |','| --- | --- | ---: | ---: | ---: |']
 for name,value,a in [('Fresh nominal8s',r,assessment),('Exact-parameter20s',long,longer_assessment)]:
  if value:lines.append(f"| {name} | {a['accepted']} | {a.get('tilt_deg',float('nan')):.3f}° | {a.get('sway_m',float('nan'))*1000:.2f}mm | {value['mean_upward_rig_support_bw']:.5f}/{value['mean_downward_rig_force_bw']:.5f}BW |")
 lines+=['',f"Nominal pressure peak{r['pressure_peak_kpa']:.2f}kPa; inverse saturation fraction{r['inverse_saturation_fraction']:.7f}; native signed-force/tension disagreement{r['native_pressure_force_max_error_N']:.2g}N. Fresh full trajectory, neural states, apparatus forces, and pressure channels bit-identical: **{report['fresh_process_bit_exact']}**. Loaded visual mesh minimums remain about−3.4mm. Full apparatus generalized loads and spatial pelvis forces/torques are logged in later `.apparatus_extended.npz` files; limited upward force does not imply independence from the other supports.", '', 'The first0.6s holds pose while the neurons settle. The original runner uses `mj_forward` there, which does not integrate activation. The pressure bank starts at zero at release; no biological actuation is retained during moving BPA physics. The20s study therefore contains19.4s free integration. This initialization is part of the stated protocol, not a measured precharged pneumatic state.','','## Fixed held-out checks','','Pressure-map/delay uncertainty was fixed before these tests and was not used to fit neural gains. Nominal sensory replay freezes all state-dependent SNS inputs while retaining DRIVE/POSTURE/POST_ bias commands and local length-aware pressure inversion. It is a test of SNS feedback compensation, not removal of the actuator interface.']
 if held:
  lines+=['','| Test | Gate | Tilt | Sway | Failure |','| --- | --- | ---: | ---: | --- |']
  for row in held['results']:
   a=row['assessment'];lines.append(f"| {row['case']} | {a['accepted']} | {a.get('tilt_deg',float('nan')):.3f}° | {a.get('sway_m',float('nan'))*1000:.2f}mm | {row['failure'] or 'none'} |")
 lines+=['','The composite80/120ms pressure response,16msdelay,20%force-map loss, and12%loading-direction memory cases exited the declared contraction domain after43physical steps (86ms), with both live feedback and replay. These aborted cases cannot establish successful recovery or failed neural compensation. The2%BW torso pulse lasts150ms on the MuJoCo physics clock at3s; live/replay both completed. A small change in neural state is not by itself evidence of improved balance. Isolated force-map and delay cases distinguish the startup/route-domain limitation.','','Records: `'+str(checkpoint)+'`, `'+str(heldpath)+'`. No manuscript numbers, original winners, or mounted CAD parameters were changed.']
 record=ROOT/'EXPERIMENT_RECORD.md';text=record.read_text(encoding='utf-8');marker='## Virtual pressure-actuated Gait comparison'
 if marker in text:text=text.split(marker)[0].rstrip()+'\n'
 record.write_text(text+'\n'+'\n'.join(lines)+'\n',encoding='utf-8')
 print(json.dumps(report,indent=2))

if __name__=='__main__':main()
