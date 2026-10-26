from pathlib import Path
import shutil, re, json, difflib, hashlib, subprocess, os, time
R=Path(r'D:\Github\Bipedal_Robot')
P=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023'
D=P/'Revision_D'; E=P/'Revision_E'
assert not E.exists(), 'Preserve existing Revision E; inspect before updating'
E.mkdir(); (E/'review').mkdir()
for p in D.iterdir():
    if p.name in ('chapters','figs'): shutil.copytree(p,E/p.name)
    elif p.suffix in ('.tex','.cls','.bib','.bst'): shutil.copy2(p,E/p.name)
s=(D/'chapters/60-sns.tex').read_text(encoding='utf-8')
(E/'review/60-sns_revision_D.tex').write_text(s,encoding='utf-8')
s=s.replace('with a separate pattern-formation layer. I added afferent feedback and heel- and toe-contact signals during the subsequent development.', 'with separate pattern-formation (PF) groups for each leg: one for the hip and one for the knee and ankle. The shared RG layer contained four neurons. We added afferent feedback at the motoneuron (MN) level.')
s=s.replace('This suspended stepping result preceded the later feedback and bilateral-RG extensions; it does not demonstrate walking with those additions.', 'The later bilateral-RG extension was a separate development; the suspended result does not establish its weight-bearing performance.')
start=s.index('I developed tools for connecting neural controllers')
end=s.index('Existing software supplied separate parts',start)
s=s[:start]+r'''The suspended walker raised a question that its observed motion could not answer: which parts of the neural controller, sensory feedback, and mechanical model were responsible for the behavior? Changing several of them together could produce movement without showing why it occurred. Before using the simulation to study control of the robot, I needed to check the translated circuits and their mechanical connections separately, then examine the coupled model with controlled changes to its inputs or mechanics.

One gap concerned the mechanical model. The Li walker used simplified biological muscles, whereas this project also required a human musculoskeletal reference and a robot model with BPA routing and force production. A controller that produced motion in one of these models could not be assumed to work in the others. The OpenSim-to-MuJoCo conversion and BPA-force calculation made those comparisons possible. Their contribution was an inspectable connection between neural output, actuator force, and sensory return. Geometry and force-direction checks establish parts of that connection; they do not establish equivalent locomotion across the models.

A second gap concerned the controller itself. Transferring a circuit required preserving its connectivity, synaptic strengths, neuron dynamics, and sensory inputs. Agreement in burst timing could otherwise conceal a changed mechanism, such as a rhythm imposed by alternating contact inputs. Matched circuit tests therefore distinguish intrinsic rhythmic activity from sensory-driven alternation. Coupled perturbation tests then ask whether a controller continues to produce supported progression when the mechanics change while its conductances remain fixed.

A third gap was practical: the mechanical and neural models used different construction environments. Scripted MuJoCo and SNS-Toolbox models supported automated evaluations, while the MATLAB mechanical models and CAD-derived mechanisms needed a graphical neural interface. The SNS-Simscape library supplied editable neural and sensory blocks for that workflow. Its verification must address both circuit behavior and exchange of forces and measurements with a mechanical model. Importing a mechanism alone answers neither question. The following sections describe the tools and the tests that establish their present scope.

'''+s[end:]
start=s.index('The tools connect musculoskeletal geometry')
end=s.index(r'\subsection{Model Conversion',start)
removed_table=s[start:end]
(E/'review/replaced_capability_table.tex').write_text(removed_table,encoding='utf-8')
s=s[:start]+r'''The tools connect musculoskeletal geometry, actuator forces, and neural circuits. Figure~\ref{fig:pipeline} shows these connections and their relationship to the proposed hardware system.

'''+s[end:]
start=s.index('MyoConverter supplies the initial')
end=s.index('For the BPA-routed robot model',start)
s=s[:start]+r'''My first conversion attempt used the robot model, which combined BPA routing on the Gait2392 skeleton with unchanged human muscle parameters. MyoConverter's muscle-equilibrium checks rejected infeasible combinations of routed length and muscle parameters. I therefore separated conversion of the human reference from reconstruction of the robot routing. This distinction matters because an error caused by an artificial-actuator route should not be interpreted as a property of the original human model.

The human-model conversion also required specific repairs. Three small hip rotators per side, the quadratus femoris, gemelli, and piriformis, were disabled at model load. The converted model retains 92 actuator entries, but these six contribute no force. The model retains the OpenSim and International Society of Biomechanics coordinate convention; frame transformations and individual-actuator tests check that joint motion and torque directions have the intended signs.

The source represents patella motion through muscle insertion locations that vary with knee angle. In the converted model, limits on the joints representing these moving path points conflicted with their coupling constraints. Removing those conflicting limits allowed the coupled path points to follow their prescribed trajectories. The rectus femoris path also required a patella-tracking via point to restore its knee-extension action. Explicit collision pairs introduced unwanted contacts, so contact was restricted to the foot--ground pairs. Passive joint torques supplied joint-stop behavior absent from the conversion. These repairs explain why conversion needed checks beyond successful loading: a model could load while producing an incorrect torque direction or unintended contact force.

Comparisons with the OpenSim attachment paths and single-actuator tests checked the transferred geometry and force directions. They do not establish that every muscle reproduces the source force throughout its range or that the controller produces walking. The neural and mechanical models exchange signals at a \qty{2}{\milli\second} solver time step.

'''+s[end:]
# Add new implementation evidence without treating period matching as waveform matching.
mark=r'\subsection{Closed-Loop Model Coupling}'
i=s.index(mark)
s=s[:i]+r'''A separate Simulink implementation of the Deng rhythm-generating circuit produced ten bursts over \qty{20}{\second} following a single \qty{10}{\nano\ampere}, \qty{1}{\milli\second} pulse. Its period was \qty{1.938}{\second}, compared with approximately \qty{1.94}{\second} in the NumPy reference. This test demonstrates self-sustained activity after a brief input, unlike the contact-driven comparison above. It does not demonstrate pointwise voltage agreement: the logged maximum voltage differences were approximately \qty{63}{\milli\volt}.

Additional tests checked the graphical implementation at different scales. In the five-neuron SNS-Toolbox tutorial circuit, final membrane voltages agreed with the NumPy implementation to within \qty{1.1e-5}{\milli\volt}. The persistent-sodium tutorial pair reproduced its winner and loser steady states, but retained transient differences. For the 410-neuron export, the maximum difference across all final membrane voltages at \qty{0.3}{\second} was \qty{4.201e-6}{\milli\volt} with a \qty{2}{\milli\second} Euler time step. This short-window check tests translation of the network equations and wiring; it does not establish long-duration agreement or transfer of a tuned walking controller.

A Li-derived contact-driven circuit was also assembled with 58 neurons and 98 connections. This implementation substitutes graded stance cells with slow adaptation for the original spiking cells. Under scripted alternating heel inputs, it produced rhythmic activity with an inter-leg correlation of $-0.271$. It is an adaptation tested with imposed contacts, rather than a reproduction of the original circuit's closed-loop walking experiment.

'''+s[i:]
start=s.index('The MuJoCo--Simulink bridge,')
end=s.index(r'\subsection{Perturbation Robustness',start)
s=s[:start]+r'''The MuJoCo--Simulink bridge was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} open-loop run without afferent return. A separate reflex experiment exercised commands and sensory measurements through the interface. These tests establish selected force and signal paths across the interface, rather than a complete transferred locomotor controller.

The Simulink air-stepping demonstration connects eight BPA-force elements to two planar legs with a fixed pelvis and no ground contact. It establishes that neural commands can drive the actuator and joint calculations together. However, the approximately \qty{0.31}{\second} swing period differed substantially from the RG period of \qty{1.938}{\second}. The observed limb oscillation therefore does not establish that the mechanics follow the neural rhythm. Mechanical imports from OpenSim and SolidWorks remain another stage of integration; the standalone reflex and air-stepping demonstrations do not establish control of those imported mechanisms.

'''+s[end:]
start=s.index('I tested the coupled W2L walker')
end=s.index(r'\section{Capabilities and Remaining Limits}',start)
(E/'review/excluded_perturbation_figure.tex').write_text(s[start:end],encoding='utf-8')
s=s[:start]+r'''The coupled W2L evaluations address a different question from the circuit-translation tests: whether a fixed controller tolerates changes in the body. The mass-change experiment held the conductances of the selected W2L controller fixed while scaling segment masses and inertias. Thigh or distal-leg properties were changed by $\pm10$ and $\pm20$ percent, and the whole leg by $\pm10$ percent. All ten conditions met the campaign's supported-progression criterion during \qty{30}{\second} runs: finite simulation, no detected fall, at least six heel-contact cycles per leg, median inter-leg phase between 0.35 and 0.65, and forward displacement greater than \qty{0.05}{\metre}. The observed forward displacements ranged from \qty{0.781}{\metre} to \qty{1.630}{\metre} over the scored windows.

These results concern the W2L mechanical model with a retained virtual harness, not the BPA-routed Gait2392 robot model. The harness supplied vertical feedback and constrained lateral motion, roll, and yaw. Its nominal weight-support setting was \qty{25}{\percent}; the realized vertical force also depended on the height-control terms. Thus the test demonstrates tolerance to the specified mechanical changes within the supported simulation, rather than unsupported balance or a constant \qty{25}{\percent} unloading force.

Separate searches reduced the nominal weight-support setting from 70 to 10 percent, with controller parameters allowed to change between settings. Each selected configuration passed a fresh-process \qty{30}{\second} evaluation. This is evidence that supported progression could be recovered at each tested setting; it is not a demonstration that one frozen controller tolerated the entire reduction. The local conductance-perturbation tests remained less robust: the selected configurations at nominal 25 and 10 percent support passed 9 and 10 of 12 probes, respectively, and both were recorded as failures of their acceptance criterion. Tolerance to segment-mass changes therefore should not be generalized to tolerance to arbitrary controller changes.

Exploratory network-deletion tests recorded 13 left-heel cycles after several pathway removals, two after commissural removal, and zero after Ib-feedback removal. However, the corresponding summary counted the left heel again for its right-side column. These records suggest tests of commissural and load-feedback contributions, but do not establish bilateral pathway necessity. Correctly scored bilateral contact and joint traces are needed before assigning that interpretation. Likewise, continued motion after a stimulus or push requires a defined recovery measure before it can establish recovery to the original gait.

'''+s[end:]
start=s.index('The AnimatLab attempt motivated a repeatable')
s=s[:start]+r'''The initial AnimatLab result showed suspended stepping-like motion but left its dependence on circuit organization, feedback, and mechanics unclear. The tools developed here make those dependencies testable in stages. Conversion repairs address whether forces act along the intended routes. Matched circuit tests address whether a translated implementation preserves selected neural responses. Coupled evaluations address how a controller behaves when its sensory input comes from the mechanical model. Keeping those questions separate prevents a successful component test from being treated as evidence for complete locomotion.

The new results strengthen specific parts of this chain. A brief-pulse test produced self-sustained neural bursts, short-window comparisons checked network translation, and fixed-controller mass changes preserved supported progression in the W2L model. They also identify remaining problems: the Li transcription changes the neuron dynamics, the air-stepping mechanics do not yet follow the RG period, and the deletion summary requires bilateral rescoring. None of these results establishes that the BPA-routed robot model reproduces the W2L result or that either model walks without a harness.

The resulting contribution is a set of neural, actuator, and mechanical interfaces with tests that distinguish implemented capabilities from unresolved model behavior. The next use of these tools is to compare controllers and feedback pathways under matched mechanical conditions, then evaluate whether the resulting behavior transfers to the robot. Chapter~\ref{ch:futurework} discusses hardware deployment and further extensions.
'''
s=s.replace('one discrete solver step','one discrete solver time step')
s=s.replace('step.', 'time step.') if False else s
(E/'chapters/60-sns.tex').write_text(s,encoding='utf-8')
before=(D/'chapters/60-sns.tex').read_text(encoding='utf-8')
(E/'review/changes.diff').write_text(''.join(difflib.unified_diff(before.splitlines(True),s.splitlines(True),fromfile='Revision_D/chapters/60-sns.tex',tofile='Revision_E/chapters/60-sns.tex')),encoding='utf-8')
def labels(x): return re.findall(r'\\label\{([^}]*)\}',x)
from collections import Counter
audit={'removed_labels':sorted(set(labels(before))-set(labels(s))), 'added_labels':sorted(set(labels(s))-set(labels(before))), 'removed_numeric_tokens':list((Counter(re.findall(r'\d+(?:\.\d+)?',before))-Counter(re.findall(r'\d+(?:\.\d+)?',s))).elements()), 'added_numeric_tokens':list((Counter(re.findall(r'\d+(?:\.\d+)?',s))-Counter(re.findall(r'\d+(?:\.\d+)?',before))).elements()),'bib_unchanged':hashlib.sha256((D/'thesis.bib').read_bytes()).hexdigest()==hashlib.sha256((E/'thesis.bib').read_bytes()).hexdigest()}
(E/'review/conservation.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print(E); print(json.dumps(audit))
