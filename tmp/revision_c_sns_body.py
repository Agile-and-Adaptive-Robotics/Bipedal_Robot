from pathlib import Path
import re,json
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
original=(C/'review/before/60-sns.tex').read_text(encoding='utf-8-sig')
def block(kind,label):
    for m in re.finditer(r'\\begin\{'+kind+r'\}(?:\[[^\]]*\])?.*?\\end\{'+kind+r'\}',original,re.S):
        if r'\label{'+label+'}' in m.group():return m.group()
    raise ValueError(label)
def figure(label,caption):
    s=block('figure',label)
    start=s.index(r'\caption');brace=s.index('{',start);depth=1;i=brace+1
    while depth:
        depth += (s[i]=='{')-(s[i]=='}');i+=1
    short=re.match(r'\\caption\[([^]]+)\]',s[start:])
    optional='['+short.group(1)+']' if short else ''
    return s[:start]+r'\caption'+optional+'{'+caption+'}'+s[i:]+'\n\n'
head=(C/'chapters/60-sns.tex').read_text().split(r'\section{A Neuromechanical Simulation Toolchain}')[0]
walker=figure('fig:animatlab_model_architecture',r'AnimatLab walker (top) and bilateral controller hierarchy (bottom). The model combines the Li walker mechanics and Deng CPG organization. These diagrams show construction; the suspended stepping result is described in the text.')
head=head.replace(r'\section{Requirements for the New Tools}',walker+r'\section{Requirements for the New Tools}',1)
main=head+r'''\section{Software Tools}\label{sec:preliminary_sim_methods}

The tools connect musculoskeletal geometry, actuator forces, and neural circuits. Table~\ref{tab:sns_capability} separates the implemented interfaces from incomplete mechanical integration. Figure~\ref{fig:pipeline} shows their relationship to the proposed hardware system.

\begin{table}[htbp]
\caption[Neuromechanical tools and verification scope]{Implemented tools, representative checks, and limits. A successful interface test does not establish walking performance.}
\label{tab:sns_capability}
\small
\begin{tabularx}{\textwidth}{@{}>{\raggedright\arraybackslash}p{0.23\textwidth}>{\raggedright\arraybackslash}p{0.34\textwidth}>{\raggedright\arraybackslash}X@{}}
\toprule
\textbf{Tool} & \textbf{Implemented capability} & \textbf{Verification and limit} \\
\midrule
MuJoCo--SNS coupling & Converted musculoskeletal model, BPA force interface, and afferent feedback & Actuator geometry and dynamic signs checked; supported closed-loop evaluations, not unsupported walking \\
SNS-Simscape library & Editable neuron, synapse, sensory, and actuator blocks & Reflex demonstrations and controller comparisons; generic BPA demonstration parameters remain provisional \\
Controller ports and test scripts & Comparable neural implementations and repeatable evaluations & Rhythm-core agreement across NumPy and Simulink; agreement is specific to the tested network and inputs \\
Mechanical imports & OpenSim-derived and SolidWorks-derived mechanisms in Simscape Multibody & Mechanism import achieved; the SolidWorks knee's welded joint and incomplete muscle integration remain limitations \\
\bottomrule
\end{tabularx}
\end{table}

\subsection{Model Conversion and Actuator Coupling}\label{sec:sns_existing_tools}

MyoConverter supplies the initial OpenSim-to-MuJoCo translation. I added repair and audit steps for the custom muscle paths and joint conventions, then coupled the converted model to an SNS-Toolbox controller. The converted model carries 92 muscle-tendon actuators. Geometry comparisons against the source and individual-actuator tests check that the conversion preserves the intended attachment paths and torque directions. The neural and mechanical models exchange signals at a \qty{2}{\milli\second} step.

The BPA interface evaluates the quasi-static actuator force model from Chapter~\ref{ch:muscles} using the simulated muscle length and pressure. It provides an actuator model for controller development, rather than a simulation of valve and pressure dynamics. The controller receives length, velocity, force, and contact signals through the corresponding sensory channels. The correction terms of Chapter~\ref{ch:torque} do not by themselves verify the transferred routing or compliance of a dynamic model; these require separate checks.

\subsection{Graphical Circuit Construction}\label{sec:sns_new_tools}

The SNS-Simscape library implements neuron, synapse, afferent-encoder, activation, and muscle-force elements as editable Simulink blocks. It extends the laboratory's earlier modeling approach \citep{hunt_development_2017, young_analyzing_2019} to the BPA interfaces used here. The relevant distinction is the added actuator and mechanical interfaces, rather than graphical neural modeling itself. Figure~\ref{fig:simulink_sns_implementation} shows the core blocks and graded-neuron signal flow.

Two mechanical import routes were investigated. The OpenSim-derived body imports as a mechanism whose muscles must be connected through actuator blocks. The SolidWorks knee imports with its intended hinge welded because the exporter does not preserve the hinge mates. Consequently, the import is not evidence of a controlled biomimetic knee. These limits are separate from the standalone reflex demonstrations.

'''
main+=figure('fig:pipeline',r'Model conversion, closed-loop simulation, and proposed hardware deployment. Solid simulation interfaces connect neural activation to BPA force and return afferent signals. The embedded valve-control and network-evaluation system is future work (Chapter~\ref{ch:futurework}).')
main+=figure('fig:simulink_sns_implementation',r'Core SNS-Simscape blocks (top) and graded-neuron implementation (bottom). Summed currents drive membrane integration; threshold, inverse-slope gain, and saturation bound the output. Voltage feeds both current calculations. The generic BPA demonstration remains provisional.')
main+=r'''\section{Verification of the Implemented Interfaces}\label{sec:preliminary_sim_results}

The tests address three distinct questions: whether a circuit produces its expected response, whether implementations agree under matched inputs, and whether the neural and mechanical models exchange signals in a closed loop. The representative results below establish those capabilities. The larger evaluation sets and their limitations are collected in Appendix~\ref{app:simulation_verification}; neuron equations and calibration details remain in Appendix~\ref{app:neuron_equations}.

\subsection{Reflex Response}\label{sec:sns_simulink_demos}

The library assembles Ia stretch excitation, reciprocal inhibition, and Ib autogenic inhibition in a reduced-order antagonist knee model (Figure~\ref{fig:simulink_reflex_topology}). A second demonstration tests the stretch reflex in an arm holding a beverage container as it is loaded. With the reflex enabled, arm sag was \ang{2.00}, compared with \ang{7.85} when it was disabled. This comparison demonstrates the response of the implemented feedback pathway. It does not validate the provisional mechanical parameters against an experiment or establish quantitative control of the physical knee.

'''
main+=figure('fig:simulink_reflex_topology',r'Reduced-order antagonist knee-reflex circuit. Triangles mark excitation, and filled circles mark inhibition. Ia input excites the homonymous motor pool and inhibits its antagonist; Ib input inhibits the homonymous pool. Motor output drives activation and BPA-force blocks. Mechanical parameters and the generic BPA block remain provisional.')
main+=r'''\subsection{Matched Neural Implementations}\label{sec:crossengine_results}

The contact-driven W2L rhythm core was implemented in NumPy and native Simulink blocks from a common connectivity specification. Both implementations received the same input protocol. Burst periods matched at \qty{1.000}{\second}, all ten expected bursts were reproduced, and voltage-trace correlation was at least 0.99996 (Figure~\ref{fig:w2l_crossengine}). These results verify agreement for this circuit and protocol, rather than interchangeability of arbitrary controllers or simulators.

Removing the phasic heel-contact inputs caused the circuit to latch in both implementations. The agreement of that negative result matters: this network's alternating activity depends on the specified sensory drive. It should not be interpreted as an autonomous oscillator merely because it alternates under the full protocol.

'''
main+=figure('fig:w2l_crossengine',r'W2L rhythm-generator voltages under matched inputs in NumPy and Simulink. Burst periods and counts agree; full-trace correlation is at least 0.99996. The residual \qty{2}{\milli\second} lag is one discrete solver step. The test concerns the driven rhythm core, not walking mechanics.')
main+=r'''\subsection{Closed-Loop Model Coupling}

The MuJoCo--Simulink bridge was tested against an independent Python integration using single-muscle activation, different control rates, and sensor readback. Readback was bit-exact over the scored window. The exported neural network also drove all 92 actuators during a \qty{5}{\second} closed-loop run. These tests establish that commands and measurements pass through the interface; they do not establish successful gait control. In particular, under tonic drive alone the Simulink realization reached a fixed point while the NumPy reference oscillated. The network-level discrepancy remains a limit on transferring tuned controller parameters between implementations.

\section{Capabilities and Remaining Limits}\label{sec:sns_discussion}

The AnimatLab attempt motivated a repeatable way to connect neural circuits to mechanical models and test their individual interfaces. The resulting tools provide converted geometry, a BPA-force interface, graphical neural components, controller ports, and automated comparisons. The reflex and rhythm-core tests demonstrate specific circuit behaviors; the bridge tests demonstrate closed-loop signal exchange. Mechanical import, controller transfer, and locomotion performance remain separate validation problems.

The supported robot evaluations produced stepping-like motion and, in some configurations, forward progression. They did not demonstrate harness-free walking. Removing a tonic neural drive in a supported model does not remove the mechanical support. Spiking-neuron comparisons likewise preserved a continuous motor-activation interface but changed rhythm timing and, in some runs, limb coordination. The detailed comparisons are reported in Appendix~\ref{app:simulation_verification}; they do not support a general claim that changing neural coding leaves movement unchanged.

The next step is to use these interfaces and tests to improve weight-bearing control and then evaluate transfer to hardware. Chapter~\ref{ch:futurework} discusses those extensions. They are not prerequisites for using the demonstrated reflex, conversion, and cross-implementation tools.
'''
(C/'chapters/60-sns.tex').write_text(main,encoding='utf-8')

app=r'''\chapter{Appendix D: Supplementary Neuromechanical Verification}\label{app:simulation_verification}
\graphicspath{{figs/60-sns/}{figs/60-sns/AnimatLab_reference/}}

This appendix supplies supporting model records and evaluation sets for Chapter~\ref{ch:sns}. It separates implementation checks, gait comparisons, pathway-removal tests, and spiking conversions because they test different claims. Neuron equations and calibration tables are provided in Appendix~\ref{app:neuron_equations}.

\section{Model Construction and Recording}

The archived AnimatLab projects provide reference circuits for implementation comparisons. Figure~\ref{fig:animatlab_controller_motifs} shows their neural organization. A later reorganized bilateral project added commissural and biarticular pathways; it ran headless but latched under the tested drives. That assembly result is not evidence of stepping. The phase-1 recordings below concern an earlier instrumented configuration, not the reorganized network.

'''
app+=figure('fig:animatlab_controller_motifs',r'Archived AnimatLab motifs: rhythm generation (A), pattern formation (B), and sensorimotor organization (C). Ext and Flx denote antagonist pathways; IN denotes interneurons, MN motor neurons, and Ia and Ib muscle afferents. Triangles and circles distinguish excitatory and inhibitory connections.')
app+=r'''A \qty{10}{\second} instrumented run recorded joint, neural, sensory, and body-motion channels at a \qty{0.2}{\milli\second} interval. Contact-sensory recordings are membrane voltages, not contact forces. The initial external horizontal force and the limited window prevent the joint traces from establishing stable walking.

'''
app+=figure('fig:animatlab_phase1_preliminary',r'Bilateral hip (A), knee (B), and ankle (C) angles in the recorded joint conventions. All 50,001 samples in the \qty{10}{\second} window are retained. An external horizontal force acts during the first second. The recordings show motion, not validated walking or ground-contact gating.')
app+=figure('fig:mujoco_model_views',r'Converted Gait2392 robot at its reference pose (A) and independent synthetic antagonist-BPA knee at \ang{-35} (B). Cyan and gold identify flexor and extensor tendon routes. The images show construction, not controlled trajectories or validated torque predictions; scales differ.')
app+=r'''\section{Controller Ports and Reference Checks}\label{sec:sns_ports}\label{sec:crossengine_methods}

The matched W2L test compares implementations from a common circuit specification under identical drive and contact inputs, using burst period, bilateral phase, burst count, and voltage correlation. The Deng controller port was checked separately: a current-pulse test produced a \qty{1.938}{\second} period against a \qty{1.94}{\second} NumPy reference, with antagonist correlation $r=-0.86$. Agreement on these tests does not establish long-term phase agreement in every network.

The full exported network agreed with its NumPy reference to \qty{4.2e-06}{\milli\volt} over the early transient across 410 neurons. Its subsequent divergence and the tonic-drive fixed-point discrepancy reported in Chapter~\ref{ch:sns} limit parameter transfer. The OpenSim-derived Simscape mechanism contains 1,446 blocks and 149 joints; mechanism compilation is distinct from completing its neural and muscle coupling.

\subsection{Archived AnimatLab Models}\label{sec:animatlab_validation}

Headless reruns tested whether archived models reproduced their reference behavior in the original simulator. Eleven targetable rhythm and gait metrics passed their acceptance bands. The bilateral contact-driven model's rhythm generators had periods of \qty{0.450}{} and \qty{0.467}{\second}; the modern W2L model showed joint oscillations at \qty{2.25}{\Hz}. The Li-model rerun reached \qty{0.637}{\meter\per\second} with step periods near \qty{1.3}{\second} \citep{li_bipedal_2017}. These are reference-model checks, not results for the combined Deng--Li model or evidence of unsupported walking.

\section{Tuning and Supported Evaluations}\label{sec:sns_tuning_methods}\label{sec:sns_curriculum}

Free network gains were searched with a tree-structured Parzen estimator, with each candidate evaluated through simulation. The staged procedure separates suspended rhythm generation, sensory feedback, standing, and ground evaluations (Figure~\ref{fig:sns_curriculum}). Finalists were checked under \qty{1}{\percent} multiplicative parameter perturbations and changes in drive. Such checks probe the neighborhood of a selected configuration; they do not establish robustness to untested loads or removal of the harness.

The cycle-shape metric assigns zero to perfect agreement. A run without detectable cycles returns $-25$ as a failure sentinel. Because this value is numerically greater than the genuine-gait scores below, failed runs must be identified before ranking scores. The sentinel itself is not a valid measure of gait quality.

'''
app+=figure('fig:sns_curriculum',r'Staged tuning procedure. Air-stepping evaluations are suspended, and the reported ground evaluations retain mechanical support. Each stage introduces additional pathways or tasks. Scores require a separate rhythm/failure check because the no-cycle sentinel cannot be ranked as a gait score.')
app+=r'''\subsection{Reference-Library Comparison}\label{sec:gait_library_results}

The 43-cycle library contains 32 running trials from eight subjects at four speeds, ten walking references, and a predicted-walking cycle. The walking references span prescribed speeds of \qtyrange{0.5}{2.0}{\meter\per\second} and self-selected trials. The four configurations were scored without treating a good tuning-reference score as evidence of generalization. Library scores were substantially more negative (Table~\ref{tab:gait_sweep}), and the evaluations showed unilateral or asymmetric gait. This comparison bounds transfer to the available references; it does not validate normal human walking.

'''
app+=block('table','tab:gait_sweep')+'\n\n'
app+=r'''\subsection{Pathway-Removal Tests}\label{sec:walker_variant_results}

Individual pathway classes were removed from each selected configuration and the model re-evaluated on its tuning reference (Table~\ref{tab:walker_ablations}). Removing Ia feedback from the production configuration eliminated detected cycles. Other removals changed the score without necessarily eliminating movement. The different response of the six-synergy variant shows that these dependencies belong to the tested network configurations, rather than constituting a general claim about biological necessity. The systematic pattern-formation and motoneuron stimulation series was not completed.

'''
app+=block('table','tab:walker_ablations')+'\n\n'
app+=r'''\section{Spiking-Neuron Comparisons}\label{sec:spiking_methods}\label{sec:spiking_results}\label{sec:disc_spiking}

The spiking variants replace afferent, rhythm-generator, and interneuron populations while retaining graded motor neurons. Synapses onto motor neurons are calibrated by mean conductance, preserving the continuous activation interface. This construction does not require the network dynamics or resulting movement to remain unchanged. Calibration and full conversion tables are in Appendix~\ref{app:neuron_equations}.

\subsection{MuJoCo Network and Ground Evaluation}

Under constant drive, the non-spiking and spiking rhythm periods were \qty{0.399}{} and \qty{0.621}{\second}, respectively (Figure~\ref{fig:spike_rhythm}). The spiking period is \qty{55.6}{\percent} longer. Timing agreement therefore requires more than copying the selected parameters. The untuned spiking variant produced bilateral cycling in a ground evaluation where the non-spiking winner cycled predominantly on one side (Table~\ref{tab:spike_ground}; Figure~\ref{fig:spike_ground}). This individual comparison motivates further controlled testing; it does not establish superiority of spiking control.

'''
app+=figure('fig:spike_rhythm',r'Constant-drive rhythm in the non-spiking network (A) and spiking variant (B). The observed periods are \qty{0.399}{} and \qty{0.621}{\second}; the latter is \qty{55.6}{\percent} longer. The panels use different neural readouts: graded population difference versus individual spiking cells.')
app+=block('table','tab:spike_ground')+'\n\n'
app+=figure('fig:spike_ground',r'Ground-evaluation joint angles at the non-spiking winner parameters. The non-spiking model cycles predominantly on the right (A, C), while the untuned spiking variant cycles both legs (B, D). Solid indigo denotes left, and dashed orange right. Both use the supported evaluation protocol.')
app+=r'''\subsection{Simulink Demonstrations}

The graded and spiking knee-reflex demonstrations were compared over matched \qty{5}{\second} runs (Table~\ref{tab:spike_knee}; Figure~\ref{fig:spike_knee}). The beverage example retained approximately the same sag, \ang{2.02} versus \ang{2.00}, compared with \ang{7.85} without the reflex. These examples demonstrate transfer of particular reflex responses after calibration, not invariance of all controlled behavior to neural coding.

'''
app+=block('table','tab:spike_knee')+'\n\n'
app+=figure('fig:spike_knee',r'Knee-reflex demonstration at the matched operating point: graded baseline (solid indigo) and spiking variant (dashed orange). Full \qty{5}{\second} run (A) and final \qty{2}{\second} window (B). The spiking response has a slightly wider oscillation band.')
app+=r'''\subsection{AnimatLab Conversions}

The tested spiking conversions did not preserve the graded models' rhythm (Table~\ref{tab:spike_animatlab}). Converting thresholds while retaining graded synapses preserved near-baseline oscillations, whereas converting synaptic amplitudes and decay changed or eliminated alternation. These results identify the synapse conversion as a limitation of this mapping. They do not demonstrate that a suitably calibrated spiking implementation cannot reproduce the reference behavior.

'''
app+=block('table','tab:spike_animatlab')+'\n\n'
app+=r'''\section{Supported BPA-Robot Evaluations}\label{sec:sns_controller_comparison}

The BPA-actuated walker exercised the chain from neural activation through the actuator force model to multi-joint motion. In the supported air-stepping run, the harness carried \qty{70}{\percent} of the \qty{411}{\newton} body weight. The period was \qty{1.542}{\second}, knee excursion was \ang{63.5}, and peak forces were \qty{880}{} and \qty{1372}{\newton}, at a peak pressure of \qty{372}{\kPa}. This was a supported simulation, not a physical walking validation.

With split left/right drives and retuned feedback, the W2L robot progressed \qty{1.72}{\meter} over \qty{30}{\second} under the support protocol. It retained progression at \qty{60}{\percent} drive and passed 12 of 12 perturbed reruns, but did not satisfy the strict contact-duty walking criterion; below that drive the gait ceased. The mechanical support was not weaned.

In the separate human-model drive-weaning study, the selected score improved from $-161.57$ to $-145.07$, with ten contact cycles per leg and contact duties of 0.71 and 0.65. It passed eight of eight perturbed reruns and retained bilateral stepping down to zero hand drive. The motion remained in place under support. These two studies distinguish reduced neural assistance from reduced mechanical assistance; neither demonstrates harness-free walking.
'''
(C/'chapters/95-AppendixD.tex').write_text(app,encoding='utf-8')
p=C/'main.tex';s=p.read_text(encoding='utf-8');assert r'\include{chapters/95-AppendixD}' not in s;s=s.replace(r'\include{chapters/94-AppendixC}',r'\include{chapters/94-AppendixC}'+'\n'+r'\include{chapters/95-AppendixD}');p.write_text(s,encoding='utf-8')
(C/'review/02_SNS_body_and_appendix.md').write_text('Main chapter rebuilt around three interface checks: reflex response, matched rhythm implementations, and closed-loop bridge exchange. Detailed sweeps, ablations, spiking comparisons, tuning, reference reruns, and supported evaluations moved to new Appendix D, with original table values retained. GUI screenshots and the catalogue of unbuilt replication candidates omitted from the compiled text; existing files remain. Repetitive construction narrative, unsupported library novelty claims, tuning/cell-identity causation, global spiking-invariance claims, and promises of easy future tuning removed. Period comparison corrected to 55.6% longer, with the underlying periods unchanged. Figure captions shortened; no pixels or Chapter 3 assets edited. Numeric/label/figure differences are audited separately against the saved original.\n')
print('Saved short SNS chapter and supplementary Appendix D')
