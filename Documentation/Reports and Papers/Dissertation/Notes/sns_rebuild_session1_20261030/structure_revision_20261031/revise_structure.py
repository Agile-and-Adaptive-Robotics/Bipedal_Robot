from pathlib import Path
import re, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle,FancyArrowPatch
root=Path(r'D:\GitHub\Bipedal_Robot')
dis=root/'Documentation/Reports and Papers/Dissertation'
out=dis/'Notes/sns_rebuild_session1_20261030/structure_revision_20261031'
out.mkdir(exist_ok=True)
s=(dis/'Archive/session1_submission_20261031/source/60-sns.tex').read_text(encoding='utf-8')
a=(dis/'Archive/session1_submission_20261031/source/95-AppendixD.tex').read_text(encoding='utf-8')
def block(label,text=s):
 p=text.index('\\label{'+label+'}')
 b=max(text.rfind('\\begin{figure}',0,p),text.rfind('\\begin{sidewaysfigure}',0,p),text.rfind('\\begin{table}',0,p))
 typ='sidewaysfigure' if text.startswith('\\begin{sidewaysfigure}',b) else ('figure' if text.startswith('\\begin{figure}',b) else 'table')
 e=text.index('\\end{'+typ+'}',p)+len('\\end{'+typ+'}')
 return text[b:e]
old_intro=s[:s.index('\\section{The AnimatLab')]
walker=block('fig:animatlab_model_architecture')
old_controller=block('fig:animatlab_shared_rg')
air_plot=block('fig:animatlab_short_response')
treadmill=block('fig:walker_fall_to_supported')
cap_table=block('tab:sns_capability')
model_table=block('tab:sns_model_configurations')
wf=block('fig:simulation_workflow')
orig_story=s[s.index('\\section{The AnimatLab'):s.index('\\section{Knowledge Gaps')]
gap=s[s.index('\\section{Knowledge Gaps'):s.index('\\section{Methods and Verification')]
methods=s[s.index('\\section{Methods and Verification'):s.index('\\subsection{Li--Deng/W2L')]
w2l=s[s.index('\\subsection{Li--Deng/W2L'):s.index('\\subsection{Gait2392-Based')]
gait=s[s.index('\\subsection{Gait2392-Based'):s.index('\\subsection{SNS-Simscape Block')]
rest=s[s.index('\\subsection{SNS-Simscape Block'):]
gait=gait.replace('\\subsection{Gait2392-Based Models}','\\section{Gait2392 and Controller Transfer}')
gait=gait.replace(model_table,'')
# Keep the numerical diagnostics intact in Appendix D and replace their duplicated main-text discussion.
diag_start=gait.index('The Gait2392-based evaluations')
diag_end=gait.index('\\FloatBarrier',diag_start)
old_diagnostics=gait[diag_start:diag_end].strip()
gait=gait[:diag_start]+r'''These controller comparisons use a common mechanical reference and the same sensory and recording interfaces. Pathway removal tests how a controller uses feedback, while inverse-kinematics overlays show the resulting joint timing and amplitude relative to human motion. Appendix~\ref{sec:gait_library_results} retains the controller comparisons, and Appendix~\ref{sec:app_ik_overlay} retains the direct trajectory comparison. The next requirement was to make these experiments practical both through graphical circuit editing and through automated simulation.

'''+gait[diag_end:]
# Compress the motivating story, retain its measured result, and keep the knowledge gap early.
intro=r'''\chapter{Tools for Neuromechanical Simulation}\label{ch:sns}
\graphicspath{{../Figures/60-sns/}{figs/60-sns/}}

This chapter follows the development of tools for transferring neural controllers between mechanical models. An initial AnimatLab walker motivated the move to Gait2392, graphical neural construction, and accelerated simulation. These capabilities were then combined with staged controller training and BPA-actuated robot models.

'''
story=r'''\section{Motivation: The AnimatLab Story}\label{sec:prior_walker}\label{sec:animatlab_results}

For a course project in ME 610: Neuromechanical Modeling, Morrow and I coupled the bipedal mechanics of \citet{li_bipedal_2017} with the two-level central pattern generator (CPG) of \citet{deng_neuromechanical_2019} in AnimatLab \citep{cofer_animatlab_2010}. The model had bilateral thigh, shank, foot, and toe segments, sagittal-plane joints, and antagonist Hill-type muscles at the hip, knee, and ankle. I added a passive toe joint. A shared rhythm generator (RG) supplied bilateral pattern-formation (PF) and motor-neuron (MN) pathways, with afferent feedback at the MN level. We adjusted currents, synaptic weights, and muscle parameters by hand and presented the suspended stepping result at the 2021 NIH BRAIN Initiative Investigators Meeting \citep{bolenTorqueValuesArtificialKnee2021June15-17}.

In a saved graded-neuron trial, alternating activity propagated through the RG, PF, and MN layers while the hip, knee, and ankle cycled at approximately \qty{2.2}{\hertz}. The aligned recording is retained in Appendix~\ref{app:simulation_verification} (Figure~\ref{fig:animatlab_short_response}). When released toward the ground, the model lost support while the neural rhythm continued (Figure~\ref{fig:animatlab_model_architecture}). The controller schematic is shown separately in Figure~\ref{fig:animatlab_shared_rg}. This contrast between suspended motion and ground loading motivated a workflow for retuning a circuit when its mechanics and sensory conditions change.

'''+walker+'\n\n'
framing=r'''\section{Knowledge Gap and Framing}\label{sec:sns_gaps}

The central problem was how to transfer an identifiable neural controller to a different mechanical model and train it efficiently. Changing muscle geometry, joint motion, inertia, or actuator type changes both the force produced by a neural command and the sensory signals returned to the circuit. The target was a detailed Gait2392 human model, followed by a robot counterpart actuated by BPAs. Comparing controllers therefore required consistent command and sensory interfaces, explicit mechanical checks, and repeatable training conditions.

I developed a workflow that combines graphical circuit construction with automated experiments. OpenSim and MyoConverter supply musculoskeletal mechanics, SNS-Toolbox supplies neural formulations, and Simulink/Simscape supplies graphical neural and mechanical modeling \citep{delp_opensim_2007,caggiano_myosuite_2022,todorov_mujoco_2012,nourse_sns_2023}. My contributions are the model implementations, custom blocks, actuator interfaces, and experiment drivers that connect these tools. The following sections build these capabilities around the Gait2392 reference before combining them in the supported robot example.

'''
# Put the standalone controller float after the gap so it cannot split the two opening sections.
framing+=old_controller+'\n\n'
tools=r'''\section{Building the Neuromechanical Simulation Workflow}\label{sec:preliminary_sim_methods}\label{sec:preliminary_sim_results}

The Gait2392 reference defines the detailed mechanical target. The workflow also retains the simpler W2L model as a test case for checking circuit calculations, sensory pathways, and controlled perturbations before interpreting full-body results. Table~\ref{tab:sns_capability} summarizes the capabilities, and Table~\ref{tab:sns_model_configurations} distinguishes controller organization from mechanical body and actuator type.

Graphical blocks support circuit construction, subnetwork inspection, and interactive stimulus application. Headless execution runs the same experiment without opening its graphical interface, allowing parameter searches and repeated perturbation tests to be automated. An application programming interface (API) exposes commands for constructing, controlling, and recording experiments; a software development kit (SDK) supplies supporting libraries and tools. Saved configurations and traces connect these automated runs to graphical review.

'''+cap_table+'\n\n'+model_table+'\n\n'+r'''Figure~\ref{fig:simulation_workflow} shows how the model-construction routes and controller interfaces connect. Chapter~\ref{ch:overall} places them within the broader simulation-to-hardware workflow (Figure~\ref{fig:pipeline}). Detailed configurations and evaluation sets are retained in Appendix~\ref{app:simulation_verification}.

'''+wf+'\n\n'
evolution=r'''I replaced the original shared RG with one RG on each side, linked through commissural pathways. The bilateral PF and MN organization was retained, while afferent feedback was extended to the RG and PF layers and heel- and toe-contact channels were added. This later controller is referred to as Walker Two-Layer (W2L).

'''
w2l=w2l.replace('The newer bilateral W2L implementation',evolution+'The newer bilateral W2L implementation',1)
w2l=w2l.replace('\\FloatBarrier',r'''Anatomical action checks, initial-pose alignment, and bilateral recruitment adjustments produced nine left and eight right qualifying swings over ten seconds in the newer supported MuJoCo implementation (Figure~\ref{fig:walker_fall_to_supported}). Passive rails and guides retained controllable support while the legs alternated.

'''+treadmill+'\n\n'+r'''These controlled comparisons establish a reusable test case. To construct and alter the same kinds of neural pathways graphically, I next developed the SNS-Simscape blocks and editor.

\FloatBarrier''',1)
# Replace pseudo-paragraph headers with continuous prose; preserve labels.
rest=rest.replace('\\paragraph{MATLAB implementation tests.}\\label{sec:w2l_matlab_attempt}','\\label{sec:w2l_matlab_attempt}')
rest=rest.replace('I later explored MATLAB as one route to those controller experiments,','I used MATLAB for graphical controller experiments,',1)
rest=rest.replace('sns_library_implementation_current.pdf','sns_library_clean.pdf')
rest=rest.replace('\\caption[Simulink SNS library blocks and neuron implementation]{\\textbf{(A)} Current SNS-Simscape library and \\textbf{(B)} native non-spiking neuron implementation. The twelve library blocks include biological-muscle and diameter-specific BPA elements. Summed currents drive membrane integration; the threshold and saturation bound neural output.}',r'''\caption[Editable SNS blocks and graded-neuron signal flow]{\textbf{(A)} The twelve implemented SNS-Simscape library elements, redrawn as a functional block catalog. \textbf{(B)} Signal flow of the implemented graded neuron: synaptic conductance and reversal-potential products form the synaptic current, which is combined with applied and leak currents before membrane integration. Membrane voltage and normalized drive are separate outputs.}''')
# Add visible, local row spacing to both the requested curriculum table and the model overview.
for label in ['tab:sns_curriculum','tab:sns_model_configurations']:
 text=rest if label=='tab:sns_curriculum' else tools
 t=block(label,text)
 lines=t.splitlines()
 lines=[line+('\n\\addlinespace[6pt]' if line.endswith('\\\\') and not line.startswith('\\textbf') else '') for line in lines]
 newer='\n'.join(lines)
 if label=='tab:sns_curriculum':
  newer=newer.replace('\\small','\\small\n\\begingroup\n\\renewcommand{\\arraystretch}{1.22}',1).replace('\\end{tabularx}','\\end{tabularx}\n\\endgroup',1)
  rest=rest.replace(t,newer)
 else: tools=tools.replace(t,newer)
combined=intro+story+framing+gait+tools+w2l+rest
combined=combined.replace('\\noindent\\textbf{Verification.} ','')
# De-duplicate the cross-engine methods/results wording without changing its measurements.
combined=combined.replace('The contact-driven W2L rhythm core was implemented in NumPy and native Simulink blocks from those same neuron definitions, synaptic connections, and parameter values. Both implementations received the same input protocol. Burst periods matched','Under these matched inputs, burst periods matched')
combined=combined.replace('The MathWorks MuJoCo--Simulink bridge exchanges activation commands and mechanical measurements between a Simulink controller and an existing MuJoCo mechanical model. I paired this bridge with SNS-Simscape and tested','To test this connection, I paired the bridge with SNS-Simscape and compared')
# All pre-existing labels survive, including the short AnimatLab plot moved to the appendix.
insert=r'''The suspended AnimatLab recording provides the neural and mechanical context for the motivating ground test in Section~\ref{sec:prior_walker}. Alternating activity propagates through RG, PF, and MN layers over the same short interval as the joint motion.

'''+air_plot+'\n\n'
a=a.replace('\\section{Controller Ports and Reference Checks}',insert+'\\section{Controller Ports and Reference Checks}',1)
a=a.replace('\\section{Tuning and Supported Evaluations}','\\section{Tuning and Supported Evaluations}\n\n'+old_diagnostics+'\n',1)
old_labels=set(re.findall(r'\\label\{([^}]+)\}',s))
new_labels=set(re.findall(r'\\label\{([^}]+)\}',combined+a))
assert not old_labels-new_labels,old_labels-new_labels
assert combined.count('\\section{')==5
assert '\\textbf{Verification.}' not in combined
for name,text in [('60-sns.tex',combined),('95-AppendixD.tex',a)]: (out/name).write_text(text,encoding='utf-8')
(out/'SOURCE_CHECKS.json').write_text(json.dumps({'labels_preserved':True,'sections':re.findall(r'\\section\{([^}]+)\}',combined),'old_source_chars':len(s),'new_source_chars':len(combined),'appendix_chars':len(a)},indent=2),encoding='utf-8')
# Clean vector library diagram, following sns_build_library.m, with no changed model equations.
plt.rcParams.update({'font.family':'Arial','font.size':10.5,'pdf.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'dejavusans'})
fig=plt.figure(figsize=(6.5,6.1)); ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,6.5),ylim=(0,6.1));ax.axis('off')
ax.text(.12,5.93,'(A)  Implemented SNS-Simscape blocks',weight='bold',va='top',fontsize=11)
labels=['Graded neuron','Chemical synapse','Synaptic sum','Spiking LIF neuron','Ia spindle','Ib tendon organ','Activation','Biological muscle','BPA force','BPA 10 mm','BPA 20 mm','BPA 40 mm']
colors=['#eef5fb']*4+['#edf6ed']*2+['#fff1df']*2+['#fff1df']*4
for i,(label,color) in enumerate(zip(labels,colors)):
 x=.13+(i%4)*1.60;y=5.00-(i//4)*.66
 box=FancyBboxPatch((x,y),1.43,.49,boxstyle='round,pad=0.025,rounding_size=.055',facecolor=color,edgecolor='#34404b',linewidth=.9);ax.add_patch(box)
 ax.text(x+.715,y+.245,label,ha='center',va='center',fontsize=10.5)
ax.text(.12,3.39,'(B)  Graded-neuron signal flow',weight='bold',va='top',fontsize=11)
def rect(x,y,w,h,label,fontsize=10.5):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015,rounding_size=.03',facecolor='white',edgecolor='#34404b',lw=.9))
 ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=fontsize)
def arrow(points):
 for p,q in zip(points[:-2],points[1:-1]):ax.plot([p[0],q[0]],[p[1],q[1]],color='#34404b',lw=.95)
 ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle='-|>',mutation_scale=8,color='#34404b',lw=.95))
rect(.15,1.90,1.75,.76,'Synaptic current\n'+r'$I_{syn}=\sum g_jE_j-V\sum g_j$')
ax.text(.15,2.91,r'Synaptic inputs: $(g_j,\;g_jE_j)$',ha='left',fontsize=10.5)
arrow([(.95,2.81),(.95,2.68)])
ax.add_patch(Circle((2.28,2.28),.16,facecolor='white',edgecolor='#34404b',lw=.95));ax.text(2.28,2.28,'+',ha='center',va='center',fontsize=12)
arrow([(1.92,2.28),(2.10,2.28)])
ax.text(2.28,2.95,r'$I_{app}$',ha='center',fontsize=11);arrow([(2.28,2.85),(2.28,2.46)])
rect(2.75,2.04,.57,.48,r'$1/C_m$');arrow([(2.46,2.28),(2.72,2.28)])
rect(3.66,2.04,.54,.48,r'$\int dt$');arrow([(3.35,2.28),(3.63,2.28)])
arrow([(4.23,2.28),(4.60,2.28)])
ax.add_patch(Circle((4.40,2.28),.025,facecolor='#34404b',edgecolor='none'))
arrow([(4.40,2.28),(4.40,2.91),(5.04,2.91)]);ax.text(5.10,2.91,r'$V$ (mV)',ha='left',va='center',fontsize=11)
rect(4.63,1.96,1.01,.64,'Threshold\nand saturation',fontsize=10.5)
arrow([(5.67,2.28),(6.15,2.28)]);ax.text(6.16,2.47,r'$S_{drive}$',ha='center',fontsize=10.5)
rect(2.69,.94,1.83,.56,r'$I_{leak}=G_m(V_{rest}-V)$')
arrow([(4.40,2.28),(4.40,1.52)])
arrow([(2.66,1.22),(2.28,1.22),(2.28,2.10)])
arrow([(4.40,1.68),(5.96,1.68),(5.96,.55),(1.02,.55),(1.02,1.87)])
ax.text(4.37,.27,r'$C_m\dot V=I_{app}+I_{syn}+I_{leak}$',ha='center',fontsize=11)
fig.savefig(out/'sns_library_clean.pdf');fig.savefig(out/'sns_library_clean.png',dpi=200);fig.savefig(out/'sns_library_clean.svg');plt.close(fig)
(out/'sns_library_clean.alt.txt').write_text('A: twelve SNS-Simscape functional blocks in a clean catalog, including graded and spiking neurons, synapse and sum, Ia and Ib encoders, activation, biological muscle, and four BPA force blocks. B: exact implemented graded-neuron current sum, membrane integration, voltage output, normalized drive, and voltage-dependent leak and synaptic current return paths.',encoding='utf-8')
print((out/'SOURCE_CHECKS.json').read_text());print('Figure generated:',out/'sns_library_clean.pdf')
