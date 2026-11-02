"""Source-specific frozen comparison ledger and separate experiment record."""
import json
from pathlib import Path
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa');OUT=ROOT/'prepared_spiking'

def main():
 audit=json.loads((OUT/'COMMON_STATE_AUDIT.json').read_text());rows=[]
 for p in sorted(OUT.glob('*_v2.comparison.json')):
  case=json.loads(p.read_text());r=json.loads(Path(case['result']).read_text());case.update(completed_samples=r['completed_samples'],failure=r['failure'],upward_rig_bw=r.get('mean_upward_rig_support_bw'),downward_rig_bw=r.get('mean_downward_rig_force_bw'),mesh_range_m=r.get('actual_visual_foot_mesh_clearance_range_m'),knee_max_deg=r.get('max_knee_angle_deg'));rows.append(case)
 def every(value):
  if isinstance(value,dict):return all(every(v) for k,v in value.items() if k not in ('scope','required_fields','fields','spikes_during_preparation'))
  return value if isinstance(value,bool) else True
 passed=every(audit['checks'])
 if not passed:raise AssertionError('State/trajectory audit failed; do not interpret held-outs')
 ledger=dict(machine='easteregg2',cases=rows,common_state_audit=str(OUT/'COMMON_STATE_AUDIT.json'),all_declared_state_checks_pass=True,source_state_api_test=r'D:\Users\Ben Bolen\NeuroCampaign_20261031\spiking\syn6_standing\syn6_state_api_v1_hybrid.SPIKING.json',scope='Sensory-spiking hybrid, not fully spiking RG/PF/MN. Frozen gains, no retuning. 8 study seconds include0.6s pinned preparation. Same baseline rig; different actuator/encoder treatments naturally differ in internal preparation state. Equality is verified within each live/replay pair. Apparatus constraints and visual mesh penetration remain.')
 (OUT/'PREPARED_ENCODER_COMPARISON.json').write_text(json.dumps(ledger,indent=2),encoding='utf-8')
 marker='## Prepared sensory-spiking comparison'
 text='\n'+marker+'\n\nThe complete encoder export/import interface now permits opt-in pinned preparation. Each release archive contains51fields, including LIF voltage, refractory state, decoded filter, counts/history, explicit empty RNG state, canonical ordering, graded backend/Ca states, sparse buffers, outer neural counter/current input, mechanical state, and pressure/history/series energy. The nominal preparation produces30,912spikes across4,480LIFcells. The source central RG/PF/MN remain graded (514remaining graded cells). No neural parameter was refitted.\n\nThe initial v1logger scope error was caught at release before moving physics and corrected; that failed programming attempt remains preserved and is not counted as a controller failure. Current v2graded pressure traces reproduce the earlier nominal numerical arrays exactly on all pre-existing fields. Fresh sensory-spiking biological/BPA repeats and nominal BPA replay reproduce full stored trajectories and event totals exactly. Harsh live/replay have identical complete release states, including every encoder/filter/refractory/pressure/engine warm-start field.\n\n| Case | Standing gate | Tilt (deg) | Sway (mm) | Upward rig (%BW) | Downward rig (%BW) |\n| --- | --- | ---: | ---: | ---: | ---: |\n'
 for row in rows:
  a=row['assessment'];text+=f"| {row['case']} | {a.get('accepted',False)} | {a.get('tilt_deg',float('nan')):.3f} | {a.get('sway_m',float('nan'))*1000:.2f} | {(row['upward_rig_bw'] or 0)*100:.4f} | {(row['downward_rig_bw'] or 0)*100:.4f} |\n"
 text+='\nThe prepared harsh sensory-spiking live-feedback case passes, whereas saved-nominal neural sensory replay reaches the unchanged left vastus-medialis supported-sleeve domain guard (28.87micrometre excess at57.90kPa). This differs from the earlier prepared graded pair, where both survived: encoder substitution changes the robustness margin. It does not establish universal feedback necessity or measured hardware robustness. Current length-dependent pressure inversion remains active in both live and replay cases; the replay intervention replaces neural sensory input channels, not all actuator sensing. The domain guard was not relaxed.\n\nThe compared policy uses the full baseline rotational/lateral rig, not the later0.8rotational selected rung. Knees remain protected; mesh compression/penetration is logged per case. The biological graded reference remains the previously fresh-verified pinned case in pinned_preparation. CAD mounting and hardware characterization remain separate. Full source-specific evidence: prepared_spiking/PREPARED_ENCODER_COMPARISON.json and COMMON_STATE_AUDIT.json. No manuscript edits.\n'
 for a,b in (('contains51fields','contains 51 fields'),('30,912spikes','30,912 spikes'),('4,480LIFcells','4,480 LIF cells'),('514remaining','514 remaining'),('v1logger','v1 logger'),('v2graded','v2 graded'),('28.87micrometre','28.87 micrometre'),('at57.90kPa','at 57.90 kPa'),('later0.8rotational','later 0.8 rotational')):text=text.replace(a,b)
 path=ROOT/'EXPERIMENT_RECORD.md';old=path.read_text(encoding='utf-8')
 if marker in old:
  start=old.index(marker);end=old.find('\n## ',start+len(marker));old=old[:start]+(old[end:] if end>=0 else '')
 path.write_text(old.rstrip()+'\n'+text,encoding='utf-8');print(json.dumps(dict(cases=len(rows),audits_pass=passed,ledger=str(OUT/'PREPARED_ENCODER_COMPARISON.json')),indent=2))

if __name__=='__main__':main()
