"""Append the current adopted checkpoint while preserving the earlier chronology."""
import json
from pathlib import Path

ROOT=Path(r'D:\Users\Ben Bolen\NeuroCampaign_20261031\gait_bpa')
HERE=Path(__file__).resolve().parent

def main():
 report=json.loads((ROOT/'BPA_COMPARISON_ADOPTED.json').read_text());fresh=report['fresh_heldout_checks']
 lines=['## Adopted checkpoint after force-interface consistency checks','',
 'Session 1, easteregg2. The earlier probe descriptions above are retained as history. This section and `BPA_COMPARISON_ADOPTED.json` identify the adopted results.', '',
 'The isolated biological standing calibration uses the syn6 family with shared bilateral weights, shared static initialization, and three fitted scalars. It is distinct from the locked Hunt controller. The biological candidate passes the standing thresholds and all eight ±1% fitted-parameter corners. Its exact-JSON fresh repeat is bit-identical.', '',
 'The virtual BPA substitution keeps the same fitted neural configuration. Every moving-physics actuator is pressure driven: 86 configured routes, plus the same six inactive short rotators. Each route uses two 20 mm actuators; both sides share dimensions and force-reference scale. Rest lengths are derived from the standing envelope. Fixed fittings (50 mm), maximum contraction (25%), total series stiffness (400 kN/m), pressure response (40/60 ms), and delay (4 ms) are declared assumptions. These are not identified CAD hardware dimensions or dynamics, and full walking range and packaging remain unverified.', '',
 '| Adopted condition | Standing gate | Maximum tilt | Sway | Upward / downward rig force |',
 '| --- | --- | ---: | ---: | ---: |']
 for row in report['adopted_results']:
  a=row['assessment']
  if a['accepted']:lines.append(f"| {row['condition']} | pass | {a['tilt_deg']:.3f}° | {a['sway_m']*1000:.2f} mm | {row['upward_rig_bw']*100:.3f}% / {row['downward_rig_bw']*100:.3f}% BW |")
  else:lines.append(f"| {row['condition']} | aborted at equilibrium-domain boundary after {row['completed_physics_samples']*.002:.3f} s of free physics | — | — | partial trace; not comparable to a complete-run mean |")
 lines+=['',
 'The nominal 8 s study repeats bit-exactly in a fresh process. The exact-parameter 20 s hold also passes; 0.6 s is pinned neural initialization, so it contains 19.4 s of moving physics. The preallocated motor-output path reproduces the entire original nominal trajectory, neural states, pressure states, and force logs bit-exactly. Its observed total trial time was 19.89 s, versus 32.04 s for the earlier scalar-output run; these are observed run timings, not a controlled end-to-end benchmark. The separate neural-only reference/pairwise speed ratios are approximately 26×, 43×, and 49× for s3k, w2lvar, and syn6.', '',
 'The isolated 20% forward force-map loss and 16 ms pressure-delay cases pass with live SNS feedback. Both successes repeat bit-exactly in separate fresh processes. Replaying nominal sensory inputs preserves command biases and the local length-aware pressure inverse, but neither replay completes a valid-domain trial. The replay neural prefixes are exactly nominal. This supports a narrow conclusion: existing myostatic/load feedback kept these virtual models inside the characterized equilibrium range and standing envelope under these two untrained changes. It does not establish hardware robustness, independent balance, or recovery after an observed fall.', '',
 'The combined doubled pressure time constants, longer delay, weaker force curve, and loading-direction memory still fail. The modest 2% BW, 150 ms torso push passes with both live feedback and replay. Their score difference is too small to claim a useful neural advantage from that pulse.', '',
 'Two implementation corrections are explicit. First, the domain test now applies to sleeve length after force-dependent series stretch, not to unloaded total route length. Thirty-six valid equilibria match the scalar BPA model within 1.59e−8 N; 24 genuinely unsupported pressure/stretch cases still reject. Second, synthetic force-map uncertainty changes the forward curve before series equilibrium is solved, while pressure inversion keeps the nominal curve. Earlier post-equilibrium force-scaling probes are superseded only where affected. All earlier raw outputs remain intact.', '',
 'The apparatus remains active: pelvis translational and rotational springs, lumbar and hip-rotation stabilizers, passive joint structures, and the protective knee extension stops with a 2° margin. The model has no knee hyperextension, but it is not harness free. Full signed forces and generalized stabilizer loads are logged in `.apparatus_extended.npz`; limited vertical assistance does not imply that other apparatus loads are negligible. The loaded visual foot meshes reach about 3.4 mm below the numerical floor. No zero-penetration or walking claim is made.', '',
 'The reflex-standing path disables the BAL_* task channels. The correction demonstrated here uses the existing myostatic/load pathways, not a new Peterka/Ting task-level balance pathway. Neural output is normalized excitation; the actuator interface maps it to a force reference and then pressure. The actual realized tension and pressure are logged separately. There is no instantaneous MN synaptic cancellation, dynamic torque allocator, or teacher-force trajectory.', '',
 'The first 0.6 s uses the original pinned initialization. `mj_forward` does not integrate biological activation; the pressure bank starts at zero at release. No biological force actuator remains active during moving BPA physics. The source study clock includes this initialization, whereas apparatus and pressure clocks report MuJoCo free-physics time.', '',
 'Metadata-only cleanup removed four W2L logging aliases from the canonical neuron-name array: the model has 888 neurons, not 892. Numerical parity traces were unchanged. The corrected mirror audit counts zero unequal pairs in stock s3k and w2lvar, and 189 in stock syn6; the experimental syn6 variant explicitly pools those pairs. The earlier s3k count of eight was an audit-mirroring error, not a real controller asymmetry.', '',
 'Full CAD pressure-actuated standing remains unestablished because the 54 mounted BPA assignments are incomplete. The CAD route screen is a zero-series-stretch geometric screen, not a force-dependent compliance or packaging validation. The full-Gait virtual pressure result is not being relabelled as CAD robot standing.', '',
 '### Files and reproduction', '',
 '- Adopted result ledger: `'+str(ROOT/'BPA_COMPARISON_ADOPTED.json')+'`.',
 '- Exact neural parameters: `'+str(ROOT/'rustuna_syn6_v1/BEST.json')+'`.',
 '- Virtual actuator specification: `'+str(ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json')+'`.',
 '- Fresh held-out confirmations: `'+str(ROOT/'BPA_HELDOUT_FRESH_CHECKS.json')+'`.',
 '- Plot provenance: `'+str(ROOT/'figures/FIGURE_NOTES.md')+'`.',
 '- Runner: `'+str(HERE/'gait_standing_family.py')+'`.', '',
 'Example on easteregg2; use a new tag to preserve completed cases:', '',
 '```powershell',
 '& "D:\\Anaconda\\envs\\myo\\python.exe" "'+str(HERE/'gait_standing_family.py')+'" --family syn6 --out "'+str(ROOT/'pressure_standing')+'" --tag NEW_UNIQUE_TAG --parameters-json "'+str(ROOT/'rustuna_syn6_v1/BEST.json')+'" --bpa-configuration "'+str(ROOT/'virtual_bpa_design/VIRTUAL_GAIT_BPA_CONFIGURATION.json')+'"',
 '```', '',
 'Read the source-specific adopted ledger before reusing an older raw result. No manuscript numbers, original production winners, or mounted hardware values were changed.']
 record=ROOT/'EXPERIMENT_RECORD.md';text=record.read_text(encoding='utf-8');marker=lines[0]
 tail=''
 if marker in text:
  text,section=text.split(marker,1);boundary=section.find('\n## ')
  if boundary>=0:tail=section[boundary:]
 record.write_text(text.rstrip()+'\n\n'+'\n'.join(lines)+'\n'+tail,encoding='utf-8')
 pointer=HERE/'EXPERIMENT_RECORD.md';old=pointer.read_text(encoding='utf-8');heading='## Full-Gait virtual pressure comparison checkpoint'
 if heading not in old:
  old+='\n'+heading+'\n\nThe separate ongoing campaign record is `'+str(record)+'` on easteregg2. It contains the constrained syn6 standing calibration, exact fresh repeats, pressure-state substitution, held-out feedback/replay comparisons, and the explicit force-interface corrections. The adopted source-specific ledger is `'+str(ROOT/'BPA_COMPARISON_ADOPTED.json')+'`. These are virtual Gait2392 actuator designs with retained support; they are not identified mounted CAD robot results. The separate CAD record remains `'+str(ROOT.parent/'cad_bpa/EXPERIMENT_RECORD.md')+'`.\n'
  pointer.write_text(old,encoding='utf-8')
 print(json.dumps(dict(record=str(record),repo_pointer=str(pointer),fresh_heldout_cases=[dict(case=x['case'],accepted=x['assessment']['accepted'],bit_exact=x['fresh_bit_exact']) for x in fresh['results']]),indent=2))

if __name__=='__main__':main()
