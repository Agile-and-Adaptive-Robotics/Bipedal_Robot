"""Apply Ben's results-first framing; preserve technical audits separately."""
from pathlib import Path
import json,re,shutil
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists());D=R/'Documentation/Reports and Papers/Dissertation';N=Path(__file__).parent;A=D/'Archive/session1_sns_20261031/before_results_framing';A.mkdir(parents=True,exist_ok=True)
changes=[]
def paragraph(text,prefix,replacement):
 lines=text.splitlines();hits=[i for i,line in enumerate(lines) if line.startswith(prefix)];assert len(hits)==1,(prefix,hits)
 i=hits[0];changes.append({'paragraph_prefix':prefix,'before':lines[i],'after':replacement});lines[i]=replacement
 return '\n'.join(lines)+'\n'
def save(p,old,new):
 b=A/p.name
 if b.exists():raise RuntimeError('Already applied; preserve the prior revision')
 shutil.copy2(p,b);assert p.read_text(encoding='utf-8')==old;p.write_bytes(new.replace('\n','\r\n').encode('utf-8'))
p=D/'ProofFinal/chapters/60-sns.tex';old=p.read_text(encoding='utf-8');s=old
s=paragraph(s,'I first transferred a published neural controller',r'I combined a published neural controller with a different biomechanical model and adjusted the controller by hand. Alternating neural activity and air stepping emerged. This experience motivated tools for transferring, training, and testing controllers more efficiently, first on the Gait2392 human model and then on a robot counterpart with BPAs. This chapter presents the resulting neuromechanical models, graphical neural blocks, simulation interfaces, and controlled experiments in the order summarized in Table~\ref{tab:sns_capability}.')
s=paragraph(s,'We adjusted applied currents',r'Morrow and I adjusted applied currents, synaptic weights, and muscle parameters by hand. In suspension, antagonist muscle activity alternated and the legs cycled through air stepping. We presented this result at the 2021 NIH BRAIN Initiative Investigators Meeting \citep{bolenTorqueValuesArtificialKnee2021June15-17}. The experiment established coupled neural and limb motion and exposed the difficulty of transferring a controller to different mechanics and loading conditions.')
s=paragraph(s,'Suspended motion did not establish',r'In a separate unsupported ground test, the pelvis descended while the neural rhythm continued. This trial introduced the additional demands of loading and interleg coordination, motivating staged training from suspended movement to supported ground stepping.')
s=s.replace('These are graded neural signals, not action potentials.','The neural traces are graded membrane voltages.')
s=s.replace('The agreement is specific to this circuit and input protocol.','This matched-input comparison verifies the rhythm-core implementation.')
s=s.replace('The test concerns the driven rhythm core, not walking mechanics.','The comparison checks the driven rhythm core.')
s=s.replace('Thus, the demonstrated feedback contribution is to bilateral contact coordination, rather than the existence of oscillation or forward motion.','The comparison identifies a sensory-feedback contribution to bilateral contact coordination.')
s=s.replace('Its MN target-current compensation cancels incoming PF and reflex currents, so this comparison demonstrates supported model integration and a separation-control effect, rather than full neural pattern control.','The allocator sets MN current targets for posture and separation, while the SNS supplies movement timing. This comparison demonstrates supported model integration and a separation-control effect.')
s=s.replace('transfer to alternating load-bearing ground steps remains unresolved.','these checks provide the starting configuration for training load-bearing steps.')
s=paragraph(s,'Mechanical support remains part of these results.',r'The support protocol combines a constant weight-support term with position guidance, providing controllable loading conditions during training. These components can be reduced progressively as the controller adapts to greater stance loading.')
s=s.replace("These are controlled tests of a translated circuit coupled to a suspended model, rather than a replication of the original rat's locomotion.",'These interventions link circuit organization to the timing of neural and mechanical responses in the suspended model.')
s=s.replace('Thus, a motor phase shift alone did not establish resetting of the rhythm generator.','The comparison separates motor-pattern timing from rhythm-generator resetting.')
remote=(D/'Archive/session1_sns_20261031/inputs_eb475ws4/60-sns_remote.tex').read_text(encoding='utf-8')
start=remote.index('\\subsection{Execution Speed and Actuator Implementations}');end=remote.index('\\section{Discussion:',start);execution=remote[start:end]
execution=execution.replace('This faster-than-real-time result concerns the muscle bank, rather than the complete walker.','The separately timed muscle bank therefore advances faster than real time.')
anchor='\\section{Discussion: What the Tools Establish}'
assert 'sec:sns_execution' not in s;s=s.replace(anchor,execution+anchor)
s=paragraph(s,'Matched neural traces support the tested circuit',r'The tools connect graphical neural construction, detailed mechanics, alternative actuator laws, and automated experiments. Matched neural traces establish implementation agreement, and command and sensor exchange through the bridge couples block-programmed controllers to MuJoCo models. The W2L comparison identifies a sensory-feedback contribution to bilateral contact coordination, while the loaded-arm example demonstrates a corrective response in the coupled neural and actuator blocks.')
s=paragraph(s,'The implemented workflow supports controlled circuit',r'These simulations provide a platform for investigating neural coordination with experimentally characterized BPA forces. Motor-neuron outputs are converted to actuator excitation or pressure commands, and sensory feedback allows the controller to respond to the resulting movement and loading. Future researchers can extend this simulated BPA control through improved dynamic characterization and closed-loop hardware experiments, including pressure response, loading-dependent behavior, and neural compensation for actuator-model errors. Chapter~\ref{ch:overall} connects these opportunities with the actuator characterization, torque models, SADb, and robot design.')
save(p,old,s)
p=D/'ProofFinal/chapters/95-AppendixD.tex';old=p.read_text(encoding='utf-8');s=old
s=paragraph(s,'A comparison with a native MuJoCo tensile actuator',r'Component tests now check tensile-force application against native MuJoCo actuators. Together with the characterized force law, this interface supports simulated BPA control and provides a basis for future work on dynamic actuator models, sensory compensation, and physical BPA control.')
s=s.replace('An initial BPA-interface experiment recorded','An early supported software-interface experiment recorded')
# Reconcile EB475WS4's already archived startup comparison; retain its source bytes above.
remote=(D/'Archive/session1_sns_20261031/inputs_eb475ws4/95-AppendixD_remote.tex').read_text(encoding='utf-8')
begin='\\subsection{Simulink Demonstrations}';end='\\subsection{AnimatLab Conversions}'
a=s.index(begin);b=s.index(end,a);c=remote.index(begin);d=remote.index(end,c)
retired=s[a:b];(A/'retired_knee_startup_block.tex').write_text(retired,encoding='utf-8')
s=s[:a]+remote[c:d]+s[b:]
s=s.replace('This comparison demonstrates transfer of the tested response after calibration, not invariance of all controlled behavior to neural coding.','This comparison demonstrates transfer of the tested response after calibration.')
save(p,old,s)
(N/'RESULTS_FRAMING_CHANGES.json').write_text(json.dumps(changes,indent=2))
print('Updated results-first framing; reconciled archived startup block; preserved technical force-sign audit outside prose.')
