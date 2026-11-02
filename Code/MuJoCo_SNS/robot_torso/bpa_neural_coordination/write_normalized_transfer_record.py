"""Small durable experiment index; raw state arrays remain outside the repository."""
import json,hashlib,shutil
from pathlib import Path
from normalized_prepared_transfer import ROOT,OUT,CHECKPOINT,configuration,F

def main():
 cfg=configuration();rows=[]
 for path in sorted(OUT.glob('*.comparison.json')):
  c=json.loads(path.read_text());r=json.loads(Path(c['result']).read_text());a=c['assessment'];rows.append(dict(case=c['case'],accepted=a.get('accepted',False),failure=r['failure'],free_s=r['completed_samples']*.002,tilt_deg=a.get('tilt_deg'),sway_mm=a.get('sway_m',0)*1000,footload_bw=r.get('mean_ground_load_bw'),upward_bw=r.get('mean_upward_rig_support_bw'),downward_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),knee_max_deg=r.get('max_knee_angle_deg'),result=c['result']))
 audit=json.loads((OUT/'COMMON_STATE_AND_REPEAT_AUDIT.json').read_text());payload=dict(machine='easteregg2',configuration=str(CHECKPOINT),configuration_sha256=hashlib.sha256(CHECKPOINT.read_bytes()).hexdigest(),cases=rows,complete_state_audit=audit)
 (OUT/'TRANSFER_SUMMARY.json').write_text(json.dumps(payload,indent=2),encoding='utf-8')
 text='# Normalized w2lvar prepared transfer — easteregg2\n\n'
 text+='Frozen standing-only configuration: PF-to-MN gain 0, aggregate-Ib motor gain 0, tonic bias 5/6. All ten positive IIIN target-group mappings remain at the recorded inverse pool counts, with bilateral tying and every IIIN edge retained. This is an explicit engineering mapping variant; production graphs and weights are unchanged.\n\n'
 text+='All nominal studies use 0.6 s pinned preparation plus 7.4 s free mechanics under the recorded apparatus. Biological and virtual-BPA runs share the same corrected pose, static activation targets, force-law reference, and baseline apparatus. Each actuator/encoder treatment develops its own internal state. Matched fresh/replay comparisons require exact equality only within the same treatment.\n\n'
 text+='Sensory-spiking replaces 280 pure-input afferent/contact neurons with 16 actual LIF cells each (4,480 cells); 608 central/interneuron/MN neurons remain graded. Complete release audits include encoder/filter/refractory/count/history state, explicit no-RNG state, graded calcium/voltage state, sparse work buffers, runner counters, actual mechanics, native activation, and pressure/history/directional-memory/series-energy state.\n\n'
 text+='| Case | Pass | Moving time (s) | Tilt (deg) | Sway (mm) | Foot load (%BW) | Up/down apparatus (%BW) |\n| --- | --- | ---: | ---: | ---: | ---: | ---: |\n'
 for row in rows:
  num=lambda x,n=3:'—' if x is None else f'{x:.{n}f}'
  text+='| '+row['case']+' | '+str(row['accepted'])+' | '+num(row['free_s'],1)+' | '+num(row['tilt_deg'])+' | '+num(row['sway_mm'])+' | '+num(None if row['footload_bw'] is None else 100*row['footload_bw'])+' | '+num(None if row['upward_bw'] is None else 100*row['upward_bw'])+'/'+num(None if row['downward_bw'] is None else 100*row['downward_bw'])+' |\n'
 failures=[row for row in rows if row['failure']]
 if failures:
  text+='\n## Preserved failures\n\n'
  for row in failures:text+='- '+row['case']+': '+row['failure']+'\n'
  text+='\nNo source/guard relaxation or postfailure resizing is applied. A failed reduction in guidance is not promoted to a checkpoint.\n'
  text+='\nThe failing pectineus remains neurally recruited, not commanded off. Its requested tension becomes slightly lower than the geometry permits inside the declared sleeve domain. The static inverse reports flag 3 and returns a zero pressure value, which the current bank consumes as a vent command. Actual pressure then crosses the forward-domain limit. `INVERSE_FEASIBILITY_AUDIT.json` retains actual pressure, excitation, tension, and adjacent-time feasibility bounds. No new saturation policy was introduced. The long-run release state and entire successful short-run prefixes are byte-identical, including actual activation, pressure, neural inputs, apparatus, and spike counts.\n'
 text+='\nThe harsh held-out actuator uses 80/120 ms pressure time constants, 16 ms delay, a 0.8 force-map factor, and directional memory amplitude 0.12 (80 ms memory, 0.03 m/s velocity scale). Common nominal pressure preparation is preserved before the held-out switch. The perturbation is a 0.05 BW forward torso pulse lasting 150 ms at physical time 3 s. Replay supplies nominal recorded neural inputs, while the local pressure inverse remains state-dependent; it is not removal of every actuator-local feedback path.\n\n'
 harsh=[r for r in rows if r['case'].startswith('harsh')]
 if len(harsh)==4 and all(r['accepted'] for r in harsh):text+='All four harsh live/replay cases pass. These tests therefore show bounded robustness but do not demonstrate that newly measured neural feedback was necessary. This differs from the separately preserved syn6 sensory-spiking harsh test, whose matched replay failed the unchanged sleeve-domain guard.\n\n'
 text+='Pressure sleeves and dynamics are declared virtual design assumptions based on static characterization, not mounted CAD identification or measured dynamic hardware parameters. Dimensions are unchanged after any failure. Native `act` bookkeeping does not produce BPA force. The legacy runner field named `act` is excitation; actual biological activation is logged separately or reconstructed through the validated native activation kernel from measured release states and all 92 commands.\n\n'
 text+='All acceptance remains supported standing only: substantial lateral/rotational apparatus guidance persists, signed apparatus forces and actual visual mesh penetration are retained, and knee extension must remain at or below zero. No gait/curriculum, independent standing, CAD packaging, or full-neural-spiking claim. The summary lists every raw result path and full physical clearance/knee values.\n\n'
 text+='## Reproduction\n\nMachine: easteregg2. Interpreter: `D:\\Anaconda\\envs\\myo\\python.exe`. Run `'+str(Path(__file__).with_name('normalized_prepared_transfer.py'))+'` for the four nominal cases, or `'+str(Path(__file__).with_name('normalized_transfer_checks.py'))+' --phase fresh|heldout|push|guidance|hold20|audit` for the bounded checks. Existing completed cases are reused. Checkpoint/dependency hashes are validated before simulation. Outputs stay in `'+str(OUT)+'`.\n'
 (OUT/'EXPERIMENT_RECORD.md').write_text(text,encoding='utf-8')
 dest=OUT/'source_snapshot';dest.mkdir(exist_ok=True);manifest={}
 sources=[Path(__file__),Path(__file__).with_name('normalized_prepared_transfer.py'),Path(__file__).with_name('normalized_transfer_checks.py'),Path(__file__).with_name('family_prepared_comparison.py'),Path(__file__).with_name('native_activation_replay.py'),Path(F.__file__),F.SPINAL/'virtual_walker_20261031/spiking_sns_inputs.py',CHECKPOINT,ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json',ROOT/'virtual_bpa_design/HELDOUT_SLOW_WEAK_MEMORY.json',ROOT/'pinned_preparation/PREPARATION_REFERENCE.json']
 sources += [Path(__file__).with_name(name) for name in ('native_bpa_force.py','actuator_interface.py','cad_interface.py')]
 sources += [Path(__file__).with_name('audit_transfer_inverse_feasibility.py'),Path(F.__file__).parents[2]/'bpa_muscle.py']
 for p in sources:
  digest=hashlib.sha256(p.read_bytes()).hexdigest();target=dest/p.name
  if target.exists() and hashlib.sha256(target.read_bytes()).hexdigest()!=digest:target=dest/(p.stem+'_'+digest[:8]+p.suffix)
  if not target.exists():shutil.copyfile(p,target)
  manifest[str(p)]=dict(sha256=digest,snapshot=str(target))
 (dest/'MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');print(json.dumps(dict(cases=len(rows),record=str(OUT/'EXPERIMENT_RECORD.md'))))
 index=ROOT/'EXPERIMENT_RECORD.md';begin='<!-- normalized-prepared-transfer:start -->';end='<!-- normalized-prepared-transfer:end -->'
 addition=begin+'\n## Frozen normalized-w2lvar prepared transfer\n\n'+str(len(rows))+' recorded cases; '+str(sum(r['accepted'] for r in rows))+' accepted short studies/repeats/controls. The four graded/sensory biological/BPA nominal 8 s studies pass, with exact fresh repeats; harsh-actuator and small forward-pulse live/replay controls also pass. These controls do not show that fresh neural feedback is necessary. Both reduced-guidance 0.8 attempts and both nominal 20 s attempts fail the unchanged pectineus sleeve-domain guard. The 20 s attempts reproduce the entire successful short-run actual-state prefix exactly and fail only after 16.278/16.566 s moving. No geometry or guard change.\n\nFull record: `'+str(OUT/'EXPERIMENT_RECORD.md')+'`. Ten-case scoped family table: `'+str(ROOT/'family_comparison/CURRENT_PREPARED_FAMILY_COMPARISON.md')+'`. Legacy act is excitation; biological activation is measured or verified native-kernel reconstruction, and BPA force is driven by actual pressure. `INVERSE_FEASIBILITY_AUDIT.json` records that nonzero pectineus recruitment requests slightly infeasible tension, the inverse returns zero pressure with flag 3, and the current bank vents toward that sentinel. No new handling policy was introduced.\n'+end
 existing=index.read_text(encoding='utf-8')
 if begin in existing:
  left,tail=existing.split(begin,1);_,right=tail.split(end,1);existing=left+addition+right
 else:existing+='\n\n'+addition+'\n'
 index.write_text(existing,encoding='utf-8')

if __name__=='__main__':main()
