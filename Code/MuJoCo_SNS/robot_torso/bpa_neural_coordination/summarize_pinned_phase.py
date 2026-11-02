"""Refresh the ongoing prepared-state section from completed experiment files."""
import json
from pathlib import Path
ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')

def main():
 folder=ROOT/'pinned_preparation';ledger=json.loads((folder/'SUMMARY.json').read_text());audit=json.loads((folder/'STATE_MATCH_AUDIT.json').read_text())
 text=['## Pinned actuator-state preparation and bounded guidance','',
 'This is a new opt-in initialization protocol. The historical main-runner warmup and prior raw results remain unchanged. Pose stays pinned while the existing neurons and pressure states advance for 0.6 s. A separate activation-only update matches native MuJoCo muscle activation integration bit-exactly over 250 test steps; full-model biological activation preparation is a separate condition.', '',
 'The pressure experiment uses a common nominal preparation for all control pairs. Harsh actuator response changes begin at release or after 2 s of free physics. The controller weights and three fitted gains stay fixed. Nominal sensory replay starts only after preparation; it does not alter the shared preparation.', '',
 f"At release, mean pressure is {audit['release_summary']['pressure_mean_kpa']:.3f} kPa, maximum pressure is {audit['release_summary']['pressure_max_kpa']:.3f} kPa, and modeled series-spring energy totals {audit['release_summary']['series_spring_energy_J']:.6f} J. Pneumatic energy is not quantified without a volume model. The nominal fresh-process gate is **{audit['nominal_exact_fresh_gate']}**.", '',
 '| Case | Standing gate | Tilt | Sway | Upward / downward rig |',
 '| --- | --- | ---: | ---: | ---: |']
 for case,row in ledger['cases'].items():
  a=row['assessment']
  if a['accepted']:text.append(f"| {case} | pass | {a['tilt_deg']:.3f}° | {a['sway_m']*1000:.2f} mm | {row['upward_rig_bw']*100:.3f}% / {row['downward_rig_bw']*100:.3f}% BW |")
  else:text.append(f"| {case} | terminated at a characterized equilibrium limit | — | — | partial trace |")
 text+=['',
 'The harsh composite passes with both live feedback and nominal sensory replay after the common preparation, at release and after standing. This does not support a claim that feedback is necessary in every prepared condition. The earlier failed composite was dominated by initialization/loading sensitivity rather than an established steady actuator-capacity limit.', '',
 'The causal controls retain the identical prepared neural state, then reset pressure alone, or pressure plus queued commands. Both reset variants fail in both live/replay controls. Before the reset, their complete prepared states match the nominal release state exactly. Within every live/replay pair, initialized qpos/qvel, neural gates, pressure/history, represented series energy, and engine warm start match exactly. Thus pressure reserve at release is a consequential variable in this virtual setup. These are deliberate initial-state controls, not calibrated venting procedures.', '',
 'Rotational guidance is scaled for the three pelvis orientation joints, three lumbar joints, and two hip-rotation stabilizers; damping scales with the square root to preserve its nominal ratio. A 20% stiffness reduction passes and repeats bit-exactly. A 50% reduction reaches an erector-spinae pressure-supported sleeve-domain limit, so the ladder stops. No stroke range, pressure bound, standing threshold, or contact geometry was relaxed. Lateral and vertical guidance remain unchanged in the retained 0.8 rung.', '',
 'Prepared starts improve the tilt/sway envelope but increase the transient visual-foot penetration to approximately 6.9 / 5.8 mm. That contact limitation is retained explicitly. Knee extension remains limited to 0° with the existing 2° protective margin. These simulations retain multiple apparatus and passive-joint contributions; they are not independent standing or locomotion.', '',
 'The common preparation is a controlled pinned initialization: pose is held rather than its pin reactions being modeled as a measured physical harness. Force and apparatus comparisons refer to the subsequent moving-physics phase. The source study clock includes preparation; the MuJoCo physics clock starts at release.', '',
 'Current raw ledger: `'+str(folder/'SUMMARY.json')+'`. Exact state/trajectory audit: `'+str(folder/'STATE_MATCH_AUDIT.json')+'`. Runner: `D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\robot_torso\bpa_neural_coordination\run_pinned_preparation.py`. The latest selected 20 s hold is listed only once its completed result is present. No tuning was performed on these preparation or guidance probes.']
 record=ROOT/'EXPERIMENT_RECORD.md';old=record.read_text(encoding='utf-8');marker=text[0]
 tail=''
 if marker in old:
  old,section=old.split(marker,1);boundary=section.find('\n## ')
  if boundary>=0:tail=section[boundary:]
 record.write_text(old.rstrip()+'\n\n'+'\n'.join(text)+'\n'+tail,encoding='utf-8');print(json.dumps(dict(record=str(record),cases=len(ledger['cases']),selected_fresh=audit.get('selected_rung_fresh',{}).get('bit_exact')),indent=2))

if __name__=='__main__':main()
