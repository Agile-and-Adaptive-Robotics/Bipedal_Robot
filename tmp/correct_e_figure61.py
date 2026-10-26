from pathlib import Path
import shutil,json,hashlib,difflib
R=Path(r'D:\Github\Bipedal_Robot');E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_E';V=E/'review/figure61_correction';V.mkdir(exist_ok=False)
source=Path(r'C:\Users\Ben Bolen\AppData\Local\Temp\codex-clipboard-630e834f-731e-4eac-b7a1-e268792ba1a2.png')
asset=E/'figs/60-sns/connor_ben_shared_rg_original.png';shutil.copy2(source,asset)
assert source.read_bytes()==asset.read_bytes()
(asset.with_suffix('.alt.txt')).write_text('Original Connor Morrow and Ben Bolen AnimatLab neural controller. A single shared four-neuron rhythm-generator layer controls bilateral pattern-formation groups for the hips and for the knees/ankles. The diagram is the original shared-RG organization, not Ben\'s later per-side RG extension. Image supplied by Ben and preserved byte-for-byte; panel lettering is added in LaTeX.\n',encoding='utf-8')
for n in ['60-sns.tex','95-AppendixD.tex']:shutil.copy2(E/'chapters'/n,V/('before_'+n))
p=E/'chapters/60-sns.tex';s=p.read_text(encoding='utf-8')
s=s.replace('Figure~\\ref{fig:animatlab_model_architecture} shows the model and controller organization.', 'Figure~\\ref{fig:animatlab_model_architecture}(A) shows the mechanical model, and Figure~\\ref{fig:animatlab_model_architecture}(B) shows the original shared-RG controller that Morrow and I built.')
s=s.replace('The later bilateral-RG extension was a separate development; the suspended result does not establish its weight-bearing performance.', 'This suspended result concerns the original shared-RG controller in Figure~\\ref{fig:animatlab_model_architecture}(B).')
old='I subsequently expanded the neural model into left- and right-side circuits with interconnected rhythm generators (RGs). I constructed the right-side subsystems from the left-side organization and linked the two RGs through commissural pathways. This later model is referred to as the Walker Two-Layer (W2L) model.'
new='Years later, I split the shared RG into one RG for each side and connected them through commissural pathways. I also extended sensory feedback to the RG and PF layers and added foot-contact feedback. This later bilateral-RG extension is referred to here as the Walker Two-Layer (W2L) model; it is distinct from the original controller shown in Figure~\\ref{fig:animatlab_model_architecture}(B).'
assert old in s;s=s.replace(old,new)
start=s.index(r'\begin{figure}[p]');end=s.index(r'\end{figure}',start)+len(r'\end{figure}')
s=s[:start]+r'''\begin{figure}[p]
\centering
\makebox[\textwidth][l]{\textbf{(A)}}\par\smallskip
\includegraphics[width=0.64\textwidth,height=0.21\textheight,keepaspectratio]{OurWalker.png}
\par\medskip
\makebox[\textwidth][l]{\textbf{(B)}}\par\smallskip
\includegraphics[width=\textwidth,height=0.50\textheight,keepaspectratio]{connor_ben_shared_rg_original.png}
\caption[Original AnimatLab walker and shared-RG controller]{Original AnimatLab model developed by Morrow and Bolen: (A) mechanical walker and (B) neural controller. The shared four-neuron RG drives bilateral PF groups, with a hip group and a knee/ankle group for each leg. This is the original controller, before the later split into one RG per side and the addition of feedback to the RG and PF layers and foot-contact feedback.}
\label{fig:animatlab_model_architecture}
\end{figure}'''+s[end:]
s=s.replace(r'\subsection{Graphical Circuit Construction}',r'\subsection{Block Programming Neural Pathways}')
s=s.replace('Passive joint torques supplied joint-stop behavior absent from the conversion.', 'The anatomical joint limits also had to be re-enabled; these driver joints were treated separately from the coupled path-point joints. Passive joint torques supplemented the joint behavior.')
old='The neural and mechanical models exchange signals at a \\qty{2}{\\milli\\second} solver time step.'
assert old in s;s=s.replace(old,old+' Appendix~\\ref{sec:app_conversion_repairs} records the individual repairs and their scope.')
p.write_text(s,encoding='utf-8')
p=E/'chapters/95-AppendixD.tex';s=p.read_text(encoding='utf-8');mark=r'\section{Model Construction and Recording}';i=s.index(mark)
addition=r'''\section{OpenSim-to-MuJoCo Conversion Repairs}\label{sec:app_conversion_repairs}

The conversion work in Section~\ref{sec:sns_existing_tools} involved two distinct problems: importing the robot's artificial-muscle routing with human muscle parameters, and repairing the converted human reference used for controller development. The first combination failed muscle-equilibrium checks because the routed lengths were incompatible with the retained human parameters. The subsequent human-model repairs described below belong to that reference-model loading path; they are not evidence that the BPA-routed model reproduces every human muscle force.

\paragraph{Joint directions and coordinate conventions.}
Single-actuator tests exposed reversed hip-flexion and hip-adduction actions. The corresponding hinge-axis directions were negated so that the muscle actions agreed with the adopted OpenSim conventions. The knee required a different treatment: its converted flexion-negative convention was retained. Rotating the knee without updating its coupled translational path-point coordinates gives a misleading geometry comparison, so those coordinates must follow the prescribed coupling during a pose or actuator test. These are specific joint and path checks, rather than a change from a left-handed to a right-handed model.

\paragraph{Small muscles and force scales.}
The quadratus femoris, gemelli, and piriformis actuators were disabled on both sides. Their six entries remain in the 92-actuator model, with no force contribution. A separate check found eight muscle-actuator entries with an active maximum-force scale of only \qty{1}{\newton}. Their active-force scales were restored from the source-model values. This correction concerns the human muscle actuators, not the BPA force law.

\paragraph{Patella-related paths and coupled joints.}
The converted moving muscle-path points use joints coupled to knee angle. Their sampled limits conflicted with the coupling constraints away from the sampled pose, producing large constraint forces and restricted motion. The coupled path-point joints were therefore left without those conflicting limits while retaining their prescribed trajectories. The anatomical driver joints were treated separately: limits were re-enabled on 14 hinges using the converted ranges. The rectus femoris route also lacked the patella tracking needed for knee-extension action; its intermediate path point was replaced by the vastus medialis patella-tracking point. The resulting force comparison and residual routing differences are discussed in Appendix~\ref{app:neuron_equations}.

\paragraph{Contact and virtual geometry.}
The converted contact section contained 19 explicit ground--body pairs. These predefined pairs are independent of the individual geometry contact flags, so changing those flags alone did not remove ground contact in a suspended test. Suspended tests remove the explicit contact pairs; ground-contact tests retain only the six foot--ground pairs. Virtual muscle-path marker geometry was also prevented from colliding. This separates intended foot loading from unintended contact involving path markers or other body segments.

The repairs establish specific geometry, sign, constraint, and contact corrections. They do not remove the rigid-tendon approximation or establish complete equivalence with the OpenSim model. The appropriate result remains the corresponding test: correct muscle action, a valid coupled trajectory, intended ground contact, or agreement in the measured force comparison. None alone establishes walking performance.

'''
# Appendix A carries the force comparison under its existing chapter label.
addition=addition.replace(r'Appendix~\ref{app:neuron_equations}',r'Appendix~\ref{app:neuron_equations}')
s=s[:i]+addition+s[i:];p.write_text(s,encoding='utf-8')
diff=''
for n in ['60-sns.tex','95-AppendixD.tex']:
 b=(V/('before_'+n)).read_text(encoding='utf-8');a=(E/'chapters'/n).read_text(encoding='utf-8');diff+=''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='before/'+n,tofile='after/'+n))
(V/'changes.diff').write_text(diff,encoding='utf-8')
(V/'asset_provenance.json').write_text(json.dumps({'source':'Ben supplied attachment in this chat','asset':str(asset.relative_to(E)),'sha256':hashlib.sha256(asset.read_bytes()).hexdigest(),'bytes_unchanged':True,'replaced_reference':'2layerImprove.png','old_asset_preserved':True,'panel_labels':'LaTeX (A)/(B), not raster modification'},indent=2),encoding='utf-8')
print('Corrected Figure 6.1, controller chronology, exact section heading, and added conversion-repair appendix.')
