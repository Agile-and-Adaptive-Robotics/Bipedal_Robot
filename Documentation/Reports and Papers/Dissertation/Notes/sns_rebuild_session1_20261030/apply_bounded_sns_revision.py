"""Apply explicit SNS narrative/figure directives; preserve pre-edit source bytes."""
from pathlib import Path
import hashlib,json,re,shutil,zipfile,xml.etree.ElementTree as ET
R=next(p for p in Path(__file__).resolve().parents if (p/'AGENTS.md').exists())
D=R/'Documentation/Reports and Papers/Dissertation';N=Path(__file__).parent;A=D/'Archive/session1_sns_20261031';A.mkdir(parents=True,exist_ok=True)
changes=[]
def read(p):return p.read_bytes().decode('utf-8').replace('\r\n','\n')
def backup(p):
 dest=A/p.relative_to(D);dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists():raise RuntimeError('Existing backup prevents accidental repeat: '+str(dest))
 shutil.copy2(p,dest)
def save(p,new,old):
 backup(p);raw=p.read_bytes();assert raw.decode('utf-8').replace('\r\n','\n')==old
 p.write_bytes(new.replace('\n','\r\n').encode('utf-8'))
 changes.append({'file':str(p),'before_sha256':hashlib.sha256(raw).hexdigest(),'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'removed_labels':sorted(set(re.findall(r'\\label\{([^}]+)\}',old))-set(re.findall(r'\\label\{([^}]+)\}',new)))})
p=D/'ProofFinal/chapters/60-sns.tex';old=read(p);s=old
start=s.index('\\section{Attempting W2L Integration in MATLAB}');end=s.index('\\section{Knowledge Gaps and the New Tools}',start)
matlab=s[start:end].replace('\\section{Attempting W2L Integration in MATLAB}','\\paragraph{MATLAB implementation tests.}',1)
replacement=r'''\section{AnimatLab Results and Lessons}\label{sec:animatlab_results}

The suspended AnimatLab model provided an initial test of the connection between neural activity and joint motion. In a saved graded-neuron trial, alternating extensor and flexor activity propagated from the rhythm generator through the hip pattern-formation and motor-neuron subsystems. Hip, knee, and ankle angles varied periodically over the same interval, at approximately \qty{2.2}{\hertz}. Figure~\ref{fig:animatlab_short_response} shows a few cycles of these neural and mechanical signals on a common time axis.

\begin{figure}[p]
\centering
\includegraphics[width=\textwidth,height=0.77\textheight,keepaspectratio]{animatlab_short_neural_mechanical.png}
\caption[Short-window AnimatLab neural and mechanical responses]{Recorded graded activity in the rhythm generator, hip pattern-formation subsystem, and hip motor neurons, followed by left hip, knee, and ankle angles from the same suspended AnimatLab trial. Extensor and flexor pathways are blue solid and orange dashed, respectively. The \qtyrange{1}{2.5}{\second} interval shows approximately three cycles. Joint angles retain the original AnimatLab conventions, including positive knee flexion. These are graded neural signals, not action potentials.}
\label{fig:animatlab_short_response}
\end{figure}

Suspended motion did not establish that the controller could support the model on the ground. In a separate unsupported ground test, the pelvis fell while the neural rhythm continued. That ground test used a different configuration from the suspended recording, rather than a release during the same trial. Together, these observations distinguish coordinated joint motion from control of the body during loading. They also show why neural oscillation alone is an insufficient criterion for successful locomotor control.

I needed a workflow that retained identifiable neural connections while testing their interaction with different mechanical models and sensory conditions. Graphical editing remained useful for inspecting circuits and applying stimuli, while repeated controller fitting and mechanical tests required an efficient programmatic interface.

'''
s=s[:start]+replacement+s[end:]
previous='Transferring a controller between biomechanical models changes its sensory and actuator interactions. A circuit that produces air stepping may require different feedback and coordination under load. I needed to vary neural organization and actuator mechanics while retaining identifiable model, input, and support conditions. Table~\\ref{tab:sns_capability} summarizes the models and tools developed for these comparisons.'
new=r'''Transferring a controller between mechanical models changes its interaction with muscle geometry, joint motion, inertia, and sensory feedback. A circuit that coordinates a simplified walker must be evaluated again when coupled to a detailed human musculoskeletal model such as Gait2392 or to a robot actuated by BPAs. I needed to vary neural organization and actuator mechanics while retaining identifiable model, input, and support conditions.

The workflow therefore required both graphical circuit development and repeated numerical experiments. Graphical blocks support inspection of neural subsystems and their connections to mechanical components. Headless execution runs an experiment without opening its graphical interface, allowing repeated tests to be automated. An application programming interface (API) exposes commands for constructing, controlling, and recording experiments; a software development kit (SDK) supplies supporting libraries and tools. These interfaces complement graphical editing when controller fitting requires many trials.

Table~\ref{tab:sns_capability} organizes the models and tools developed for these requirements. The Li--Deng/W2L models provide initial controller tests, and Gait2392 adds detailed musculoskeletal mechanics. SNS-Simscape supplies graphical neural blocks and mechanical coupling, including a CAD-derived robot assembly. The MuJoCo--Simulink interface exchanges control and sensory signals, while experiment drivers support staged training, perturbation tests, and alternative neuronal implementations.'''
assert s.count(previous)==1;s=s.replace(previous,new)
anchor=s.index('\\subsection{SNS-Simscape Block Library and Editor}');insert=s.index('\n\n',anchor)+2;s=s[:insert]+matlab+s[insert:]
s=s.replace('The residual \\qty{2}{\\milli\\second} lag is one discrete solver step. The test concerns the driven rhythm core, not walking mechanics.','Extensor and flexor traces are paired by leg over a common three-second interval; voltages are relative to the resting reference. The residual \\qty{2}{\\milli\\second} lag is one discrete solver step. The test concerns the driven rhythm core, not walking mechanics.')
s=s.replace('\\citet{hunt_development_2017}: intrinsic rhythm generation','\\citet{hunt_inrc_2019}: intrinsic rhythm generation')
stage='4 & Ground walking & Ground contact and forward-progression task & Joint kinematics against a human reference cycle'
assert stage in s;s=s.replace(stage,'4 & Harness-supported ground stepping & Ground contact with measured external support & Alternating stance load, swing-foot unloading and clearance, forward progression, and joint kinematics')
lit=s.index('\\subsection{Literature Circuits and Spiking Variants}');insert=s.index('\n\n',lit)+2
s=s[:insert]+r'''Biarticular mechanics require explicit consideration when organizing the pattern-formation layer. \citet{deng_biarticular_2023} tested a rat-hindlimb model without ground contact and compared active, passive-only, and disabled biarticular muscles during perturbations. Their results motivate testing multi-joint coordination separately from monoarticular control, while retaining ground-contact tests as a distinct requirement for walking.

'''+s[insert:]
save(p,s,old)
p=D/'ProofFinal/chapters/95-AppendixD.tex';old=read(p);s=old.replace('Archived AnimatLab motifs: rhythm generation (A), pattern formation (B), and sensorimotor organization (C).',r'Archived AnimatLab neural subsystems: \textbf{(A)} rhythm generation, \textbf{(B)} pattern formation, and \textbf{(C)} sensorimotor organization.')
s=s.replace('{spiking_mirror_rhythm.pdf}','{spiking_mirror_rhythm.pdf}')
s=s.replace('The panels use different neural readouts: graded population difference versus individual spiking cells.','The panels show a common \\qtyrange{2}{4.5}{\\second} interval using different neural readouts: graded population difference versus spike events reconstructed from cumulative-count increments. Periods were estimated from the full recordings. Mechanical responses are reported separately in Figure~\\ref{fig:spike_ground}.')
assert s.count('\\begin{table}[h]')==1;s=s.replace('\\begin{table}[h]','\\begin{table}[htbp]')
save(p,s,old)
p=D/'ProofFinal/thesis.bib';backup(p);raw=p.read_bytes();entries=b''
if b'{hunt_inrc_2019,' not in raw:entries+=b'\r\n@misc{hunt_inrc_2019,\r\n  author = {Hunt, Alexander J.},\r\n  title = {Real-time control of dynamic walking in a dog-like robot},\r\n  year = {2019},\r\n  month = feb,\r\n  howpublished = {Presentation at the 2nd Intel Neuromorphic Research Community (INRC2) Meeting},\r\n  note = {February 15, 2019}\r\n}\r\n'
if b'{deng_biarticular_2023,' not in raw:entries+=b'\r\n@incollection{deng_biarticular_2023,\r\n  author = {Deng, Kaiyu and Hunt, Alexander J. and Chiel, Hillel J. and Quinn, Roger D.},\r\n  title = {Biarticular Muscles Improve the Stability of a Neuromechanical Model of the Rat Hindlimb},\r\n  booktitle = {Biomimetic and Biohybrid Systems},\r\n  year = {2023},\r\n  pages = {20--37},\r\n  publisher = {Springer},\r\n  doi = {10.1007/978-3-031-39504-8_2}\r\n}\r\n'
p.write_bytes(raw+entries);assert p.read_bytes()[:len(raw)]==raw
changes.append({'file':str(p),'operation':'ASCII-only append; original bytes preserved','added_keys':['hunt_inrc_2019','deng_biarticular_2023']})
for src,dst in [('animatlab_short_neural_mechanical_preview.png','60-sns/animatlab_short_neural_mechanical.png'),('fig6_3_compact_engine_overlay.png','60-sns/sns_w2l_traces.png'),('spiking_mirror_rhythm_short.pdf','95-AppendixD/spiking_mirror_rhythm.pdf')]:
 target=D/'Figures'/dst
 if target.exists():backup(target)
 shutil.copy2(N/'previews'/src,target);changes.append({'file':str(target),'source':str(N/'previews'/src)})
(N/'APPLIED_CHANGES_20261031.json').write_text(json.dumps(changes,indent=2))
# Extract the source slide text for subsequent method-wording verification.
ppt=Path('C:/Users/Ben Bolen/Downloads/19_02_15 - INRC.pptx')
if ppt.exists():
 with zipfile.ZipFile(ppt) as z:
  slides=[]
  for name in sorted((n for n in z.namelist() if re.fullmatch(r'ppt/slides/slide\d+\.xml',n)),key=lambda n:int(re.search(r'slide(\d+)',n).group(1))):
   texts=[t.text or '' for t in ET.fromstring(z.read(name)).iter() if t.tag.endswith('}t')];slides.append(name+'\n'+'\n'.join(texts))
 (N/'HUNT_2019_SOURCE_SLIDES.txt').write_text('\n\n'.join(slides),encoding='utf-8')
print(json.dumps(changes,indent=2))
