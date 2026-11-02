"""Read-only numerical comparison and explicit standing checkpoint packaging."""
import json,sys
from pathlib import Path
import numpy as np
from gait_standing_objective import score

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main():
 best=json.loads((ROOT/'rustuna_syn6_v1/BEST.json').read_text())
 source=Path(best['trial']['diagnostics']['result'])
 fresh=ROOT/'standing/syn6_rustuna_winner_json_fresh_v1.result.json'
 a=json.loads(source.read_text());b=json.loads(fresh.read_text())
 checks={}
 for suffix,keys in (('.npz',('t','q','qfull','ctrl','neuro')),('.apparatus.npz',('samples','neuro'))):
  za=np.load(source.with_suffix('').with_suffix(suffix),allow_pickle=True)
  zb=np.load(fresh.with_suffix('').with_suffix(suffix),allow_pickle=True)
  for key in keys:
   if key not in za.files:continue
   x,y=za[key],zb[key];same=x.shape==y.shape and np.array_equal(x,y)
   delta=float(np.max(np.abs(x-y))) if x.shape==y.shape and x.size else None
   checks[suffix+':'+key]={'bit_exact':same,'max_abs_difference':delta,'shape':list(x.shape),'within_repeat_tolerance':bool(delta is not None and delta<=1e-9)}
 objective,assessment=score(b,fresh.with_suffix('').with_suffix(''))
 exact=all(v['bit_exact'] for v in checks.values());repeat=all(v['within_repeat_tolerance'] for v in checks.values())
 report=dict(machine='easteregg2',family='syn6',actuation='Native biological Hill muscles',source_result=str(source),fresh_result=str(fresh),parameters=best['trial']['params'],fresh_process_bit_exact=exact,comparisons=checks,assessment=assessment,
   measured_upward_support_bw=b['mean_upward_rig_support_bw'],measured_downward_rig_force_bw=b['mean_downward_rig_force_bw'],measured_net_rig_force_bw=b['mean_net_rig_force_bw'],mean_ground_load_bw=b['mean_ground_load_bw'],visual_foot_mesh_clearance_range_m=b['actual_visual_foot_mesh_clearance_range_m'],max_knee_angle_deg=b['max_knee_angle_deg'],
   repeat_atol=1e-9,repeat_criterion_scope='Tight absolute full-trajectory numeric agreement in each stored unit; not a bit-exact claim. Numerical backend parity is separately exactly0mV.',fresh_process_numeric_repeat=repeat,
   accepted_standing_candidate=bool(repeat and assessment['accepted']),limitations=['Eight-second standing mode, not locomotion or the locked Hunt checkpoint.','Pelvis lateral/rotational apparatus remains; not harness-free.','Visual mesh compresses/penetrates the ground by about5mm under load; contact fidelity remains qualified.','No BPA pressure model used in this case.','No held-out or longer-duration robustness claim.'],
   command='D:\\Anaconda\\envs\\myo\\python.exe '+str(Path(__file__).with_name('gait_standing_family.py'))+' --family syn6 --out "'+str(ROOT/'standing')+'" --tag NEW_UNIQUE_TAG --parameters-json "'+str(ROOT/'rustuna_syn6_v1/BEST.json')+'"')
 target=ROOT/'STANDING_CHECKPOINT.json';target.write_text(json.dumps(report,indent=2),encoding='utf-8')
 record=ROOT/'EXPERIMENT_RECORD.md';text=record.read_text(encoding='utf-8');marker='## Shared standing calibration and fresh confirmation'
 if marker in text:text=text.split(marker)[0].rstrip()+'\n'
 if marker not in text:
  text+='\n'+marker+'\n\nThe matched per-muscle Ia/II/Ib input ablations did not solve standing. A subsequent current decomposition showed tonic PF→MN and central load-group drive remain active during standing; per-muscle input suppression does not silence aggregate load pathways. Shared PF→MN ablation improved the static posture, especially for syn6.\n\nRustuna0.2.1 completed12trials with seed100031, reusing2completed observations rather than repeating them. Three consumed, bilateral variables were fitted: PF→MN conductance multiplier, Ib-group→MN multiplier, and standing-bias scale. The best parameters and exact fresh-process comparisons are in `'+str(target)+'`.\n\n'
  text+=f"Fresh8-second standing acceptance: **{report['accepted_standing_candidate']}**. Maximum tilt{assessment.get('tilt_deg',float('nan')):.3f}°, sway{assessment.get('sway_m',float('nan'))*1000:.2f}mm, measured foot load{b['mean_ground_load_bw']:.5f}BW, upward apparatus force{b['mean_upward_rig_support_bw']:.5f}BW, and downward apparatus force{b['mean_downward_rig_force_bw']:.5f}BW. Knee maxima{b['max_knee_angle_deg']}. Full numeric trajectory/state comparisons agree within1e−9stored units: **{repeat}**; bit-exact: **{exact}**. Exact per-array differences are retained in the checkpoint.\n\n"
  text+='An earlier CLI repeat rounded bare PowerShell numeric arguments to about15digits, yielding differences up to2.95e−11N rather than exact equality. That run is preserved separately. The adopted fresh repeat reads every parameter directly from BEST.json, avoiding this shell conversion; use `--parameters-json` or quote scalar arguments.\n\n'
  text+='This is a new constrained syn6 **standing-mode** calibration, not an unchanged NMF circuit, the Hunt controller, a full curriculum pass, or independent standing. The apparatus still constrains lateral/rotational motion. Visual mesh minimums are approximately−5mm under load; no zero-penetration/contact-fidelity claim is made. Longer-duration and held-out checks are next, before actuator substitution. No paper numbers or original winners were changed.\n'
  record.write_text(text,encoding='utf-8')
 print(json.dumps(report,indent=2))

if __name__=='__main__':main()
