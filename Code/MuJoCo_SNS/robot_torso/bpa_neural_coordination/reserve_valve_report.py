"""Circuit/output map and scoped report for the explicit neural valve variant."""
import json,hashlib,shutil
from pathlib import Path
import numpy as np
from reserve_valve_runner import DEST,ROOT

def draw(destination=None,training_ranges=False):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.patches import Rectangle,Circle,FancyArrowPatch
 plt.rcParams.update({'font.family':'Arial','font.size':10,'svg.fonttype':'none','pdf.fonttype':42})
 fig,ax=plt.subplots(figsize=(14,7));ax.set(xlim=(0,1400),ylim=(0,700));ax.invert_yaxis();ax.axis('off')
 def box(x,y,w,h,lines,color='#EEEEEE'):
  ax.add_patch(Rectangle((x,y),w,h,facecolor=color,edgecolor='#333333',linewidth=1.3))
  for i,s in enumerate(lines):ax.text(x+w/2,y+h/2+(i-(len(lines)-1)/2)*22,s,ha='center',va='center',fontsize=10)
 def arrow(points,color='#333333',dashed=False):
  for a,b in zip(points[:-2],points[1:-1]):ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=1.4,ls='--' if dashed else '-')
  ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='-|>',mutation_scale=12,color=color,linewidth=1.4,linestyle='--' if dashed else '-'))
 ax.text(40,35,'Explicit pressure-reserve neurons → valve-reference increments',fontsize=16,fontweight='bold')
 box(45,100,240,110,['Frozen standing configuration','Normalized w2lvar core','Core weights retained'])
 box(335,100,280,110,['Original MN output map','Force → nominal pressure','Nearest-feasible reference'])
 ax.add_patch(Circle((730,155),28,facecolor='white',edgecolor='#333333',linewidth=1.3));ax.text(730,155,'+',fontsize=24,ha='center',va='center')
 box(820,100,240,110,['Final reference cap: 620 kPa','History + delay','Unchanged fill/vent dynamics','Actual pressure state'])
 box(1100,100,250,110,['86 virtual BPA routes','Series equilibrium + model','Unchanged domain guard'])
 arrow([(285,155),(335,155)]);arrow([(615,155),(702,155)]);arrow([(758,155),(820,155)]);arrow([(1060,155),(1100,155)])
 box(45,350,240,110,['Measured model signals','Actual pressure P','Route length L and velocity v'])
 predictor=['Closing-reserve predictor','h = l0 + b + Fmax P/(620K) − L','10 ms pressure-rate filter','Preview: 80–180 ms' if training_ranges else '100 ms closing preview','Reserve target: 1 N','Error gain: 0.5–2' if training_ranges else 'Error gain: 1']
 box(335,325,280,160,predictor)
 for y,side in ((355,'r'),(505,'l')):
  ax.add_patch(Circle((730,y),45,facecolor='#DCEAF6',edgecolor='#4477AA',linewidth=1.8));ax.text(730,y,'RESERVE\npect_'+side+'\nτ = 6 ms',ha='center',va='center',fontsize=10)
 box(820,325,240,225,['Neural voltage → reference map','0–10 mV → 0–Pcap kPa','Training Pcap = 1–5 kPa' if training_ranges else 'Selected Pcap = 1 kPa','Identical map on both sides','Reference increment only'],'#E5F0E5')
 arrow([(285,405),(335,405)]);arrow([(615,405),(650,405),(650,355),(685,355)],'#4477AA');arrow([(650,405),(650,505),(685,505)],'#4477AA')
 arrow([(775,355),(820,355)],'#4477AA');arrow([(775,505),(820,505)],'#4477AA');arrow([(940,325),(940,275),(730,275),(730,183)],'#228833')
 ax.text(1030,263,'ΔP only for pect_r / pect_l',ha='center',fontsize=10,color='#228833')
 arrow([(940,210),(940,295),(165,295),(165,350)],'#888888',True);ax.text(520,286,'Actual pressure feedback',ha='center',fontsize=10,color='#666666')
 arrow([(1225,210),(1225,585),(165,585),(165,460)],'#888888',True);ax.text(695,577,'Model route-length and velocity feedback',ha='center',fontsize=10,color='#666666')
 arrow([(1225,100),(1225,70),(165,70),(165,100)],'#888888',True);ax.text(680,62,'Existing proprioceptive and contact feedback retained',ha='center',fontsize=10,color='#666666')
 ax.text(40,630,'New: two leaky relay neurons. Sensor and voltage-to-pressure maps are engineering transducers; parameters are mirrored.',fontsize=10)
 ax.text(40,657,'No added current to original MNs, no actual-pressure projection, and no force floor. All other reference channels have ΔP = 0.',fontsize=10)
 ax.text(40,684,'Supported standing of a virtual human model; assumed actuator design, not CAD hardware identification. Six pruned routes disabled.',fontsize=10)
 fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.01)
 target=Path(destination) if destination is not None else DEST
 for ext in ('svg','pdf','png'):fig.savefig(target/('CIRCUIT.'+ext),dpi=180,facecolor='white')
 plt.close(fig)

def main():
 rows=[]
 for file in sorted(DEST.glob('*.comparison.json')):
  c=json.loads(file.read_text());r=json.loads(Path(c['result']).read_text());rows.append(dict(case=Path(c['result']).name.removesuffix('.result.json'),condition=c['condition'],encoding=c['encoding'],pressure_cap_kpa=c['pressure_cap_kpa'],reserve_target_N=c.get('reserve_target_N',0.),accepted=c['assessment']['accepted'],moving_s=r['completed_samples']*.002,failure=r['failure'],tilt_deg=c['assessment'].get('tilt_deg'),sway_mm=1000*c['assessment'].get('sway_m',0.),footload_bw=r.get('mean_ground_load_bw'),upward_bw=r.get('mean_upward_rig_support_bw'),downward_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),knee_max_deg=r.get('max_knee_angle_deg'),maximum_increment_kpa=c['maximum_applied_increment_kpa'],neural_restore=c['complete_neural_restore_verified'],release_exact=all(c['common_release_exact'].values()),result=c['result']))
 selected=json.loads((DEST/'SELECTED.json').read_text());audit=json.loads((DEST/'COMPLETE_STATE_AUDIT.json').read_text()) if (DEST/'COMPLETE_STATE_AUDIT.json').exists() else {}
 report=dict(machine='easteregg2',selected=selected,rows=rows,complete_state_audit=audit,unit_checks=json.loads((DEST/'VALVE_NEURAL_UNIT_CHECKS.json').read_text()))
 checkpoint=dict(machine='easteregg2',status='nominal_supported_standing_and_pulse_pass_generalization_rejected',pressure_cap_kpa=1.,reserve_target_N=1.,preview_s=.1,underlying_configuration=str(ROOT/'group_diagnosis/NORMALIZED_STANDING_CHECKPOINT.json'),selected_nominal_result=selected['result'],fresh_result=str(DEST/'w2lvar_winner_fresh.result.json'),heldout_plan=str(DEST/'PREDECLARED_HELDOUT.json'),no_postheldout_reselection=True,not_promoted_to_curriculum=True,global_current_checkpoint_unchanged=str(ROOT/'CURRENT_CHECKPOINT.json'),record=str(DEST/'EXPERIMENT_RECORD.md'))
 (DEST/'NOMINAL_CHECKPOINT.json').write_text(json.dumps(checkpoint,indent=2),encoding='utf-8')
 (DEST/'SUMMARY.json').write_text(json.dumps(report,indent=2),encoding='utf-8');(DEST/'CIRCUIT.json').write_text(json.dumps(selected['circuit'],indent=2),encoding='utf-8');draw()
 text='# Explicit reserve-to-valve neural reference — easteregg2\n\n'
 text+='The original normalized-w2lvar network controls nominal pressure. Two named mirrored leaky neurons, RESERVE_pect_r/l, supply an explicitly bounded reference increment through a linear 0–10 mV output map. This sum occurs before the original command history, pressure delay, fill/vent dynamics, and forward force/length guard. No original MN current, output voltage, actual pressure, force, or geometry is overwritten. Only the two pectineus reference channels receive the increment; other84routes retain the same control map.\n\n'
 text+='Each relay uses 1 µS leak, 6 nF capacitance (6 ms), and bounded 0–10 nA sensor input, integrated every2ms. Nominal-map reserve h=l0+b+Fmax Pactual/(620K)−L and hdot use actual pressure, measured route length/velocity, and a10ms pressure-rate filter. Closing prediction is h+100ms min(hdot,0). The declared sensor target is0or1N equivalent series reserve;1N corresponds to2.5µm at400kN/m. This is an explicit neural input-error setpoint, not an imposed actual force or pressure floor. The source uses no held-out force-gain oracle.\n\n'
 text+='Six nominal-only fits crossed Pcap0.5/1/2kPa with reserve targets0/1N. Selection was fixed beforehand: accepted cases only, smallestPcap then smallestsetpoint. Selected:1kPa cap,1N target,100ms preview. Cap1kPa with zero target failed; cap2kPa with either target passed. The failed direct-MN-current relay and no-output comparator remain preserved.\n\n'
 text+='Independent checks find0.0mV relay error versus SNS-Toolbox and0.0kPa output-map error. Existing-MN injected current is identicallyzero. Every actual pressure state for all86routes is reconstructed exactly from common release pressure/history, final reference commands, and original delayed first-order hydraulic equations (0.0kPa error). Thus the passing result does not bypass actual pressure dynamics. Zero output recovers the physical no-feedback comparator exactly. Complete genuine nonzero-output state-restore tests and fresh-process audits include neural/relay/filter/input-buffer states, hydraulic history, direction state, and mechanical integration state.\n\n'
 text+='| Case | Condition | Pcap (kPa) | Target (N) | Pass | Moving (s) | Max increment (kPa) |\n| --- | --- | ---: | ---: | --- | ---: | ---: |\n'
 for row in rows:text+=f"| {row['case']} | {row['condition']} | {row['pressure_cap_kpa']:.1f} | {row['reserve_target_N']:.0f} | {row['accepted']} | {row['moving_s']:.3f} | {row['maximum_increment_kpa']:.5f} |\n"
 text+='\nStudies are20s including0.6s pinned preparation (19.4s moving when complete). Strong lateral/rotational apparatus guidance remains. Signed apparatus loads, actual visual mesh penetration, and knee maxima are retained in SUMMARY.json and original results. These are virtual human-model/BPA standing experiments, not CAD hardware calibration, walking, or independent balance. Hydraulic times/constants and sleeve dimensions remain declared engineering assumptions.\n\n'
 text+='Held-outs are predeclared delay, 0.8 force map, forward pulse, and combined slow/weak/directional-memory cases with paired output-off controls. Actuator-parameter changes apply at physical time zero after common nominal preparation; paired on/off cases undergo the identical change. Equal preparation does not mean identical elastic energy after changing the force law. Selection is frozen after nominal fits, and no replacement is chosen after a held-out failure. The selected case passes nominal and pulse studies but fails delay, map, and combined tests; general robustness is rejected. All failures remain explicit. Directional-memory uncertainty in standing is not validation of a measured dynamic hysteresis law.\n\n'
 for row in rows:
  if row['failure']:text+='- '+row['case']+': '+row['failure']+'\n'
 text+='\n## Reproduction\n\nMachine easteregg2; `D:\\Anaconda\\envs\\myo\\python.exe`. Runner: `'+str(Path(__file__).with_name('reserve_valve_runner.py'))+' --fit`; bounded follow-ups: `'+str(Path(__file__).with_name('reserve_valve_checks.py'))+' --phase fresh|heldout|spike|audit`; independent neural/hydraulic checks: `'+str(Path(__file__).with_name('check_valve_neural.py'))+' --tag fit_p1_r1`. Existing completed cases are reused by batch drivers.\n'
 (DEST/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')
 repeat='# Repeat the bounded neural valve standing example\n\nMachine: **easteregg2**. This is the normalized-w2lvar virtual human/BPA model, not the CAD robot and not the Hunt walking curriculum. Nominal and pulse20s pass; delay/map/combined generalization fails. Retain all qualifiers in NOMINAL_CHECKPOINT.json.\n\nFrom PowerShell, use a new dot-free tag each time:\n\n```powershell\n& "D:\\Anaconda\\envs\\myo\\python.exe" "'+str(Path(__file__).with_name('reserve_valve_runner.py'))+'" --cap 1 --preview 0.1 --reserve-target 1 --condition nominal --duration 20 --tag manual_repeat_01\n```\n\nOutputs are `'+str(DEST)+'\\w2lvar_manual_repeat_01.*`. Existing completed tags are protected against overwrite. Replace `--condition nominal` with `delay`, `map`, `pulse`, or `harsh` only to reproduce the documented fixed-condition tests; failures are retained and do not count as acceptance. `--cap 0` is the explicit output-off comparator. The optional `--encoding spike` freezes the same configuration and changes only pure-input sensory encoding; it has not been run for this valve architecture in the present record.\n\nCircuit: CIRCUIT.png/.pdf/.svg and CIRCUIT.json. Complete-state audits: COMPLETE_STATE_AUDIT.json. Exact neuron/output/hydraulic checks: VALVE_NEURAL_UNIT_CHECKS.json. Source snapshots and checksums: source_snapshot\\MANIFEST.json. Strong apparatus guidance, visual mesh penetration, and all assumed actuator parameters remain explicit.\n'
 (DEST/'HOW_TO_REPEAT.md').write_text(repeat,encoding='utf-8')
 dest=DEST/'source_snapshot';dest.mkdir(exist_ok=True);manifest={}
 for name in ('pressure_valve_neural.py','reserve_valve_runner.py','reserve_valve_checks.py','check_valve_neural.py','reserve_valve_report.py','pressure_reserve_neural.py','inverse_feasibility_policy.py'):
  p=Path(__file__).with_name(name);digest=hashlib.sha256(p.read_bytes()).hexdigest();target=dest/(p.stem+'_'+digest[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  manifest[str(p)]=dict(sha256=digest,snapshot=str(target))
 (dest/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(dict(cases=len(rows),record=str(DEST/'EXPERIMENT_RECORD.md'),circuit=str(DEST/'CIRCUIT.png'))))
 index=ROOT/'EXPERIMENT_RECORD.md';begin='<!-- pressure-reserve-valve:start -->';end='<!-- pressure-reserve-valve:end -->'
 addition=begin+'\n## Explicit neural reserve-to-valve pressure reference\n\nRecord and inspected circuit: `'+str(DEST/'EXPERIMENT_RECORD.md')+'`, `'+str(DEST/'CIRCUIT.png')+'`. The distinct two-cell output map adds a bounded reference increment before unchanged hydraulic dynamics, with no current injected into original MNs. Six predeclared nominal fits select 1 kPa output cap, 1 N-equivalent reserve target, and 100 ms preview. Nominal20s and fresh repeat pass exactly; 0.05 BW/150ms forward-pulse20s also passes while output-off fails. Longer delay,0.8force-map,and combined uncertainty fail; the configuration is rejected for general robustness and is not reselected after held-outs. Independent SNS and voltage-output checks have zero error; all86actual pressure trajectories reconstruct exactly under the original hydraulic law. Physical/neuronal preparation and fresh full final states are audited. Strong apparatus guidance, up to4.71mm visual mesh penetration, and declared virtual actuator assumptions remain. No CAD hardware, independent balance, gait, or Hunt-curriculum promotion claim.\n'+end
 existing=index.read_text(encoding='utf-8')
 if begin in existing:
  left,tail=existing.split(begin,1);_,right=tail.split(end,1);existing=left+addition+right
 else:existing+='\n\n'+addition+'\n'
 index.write_text(existing,encoding='utf-8')

if __name__=='__main__':main()
