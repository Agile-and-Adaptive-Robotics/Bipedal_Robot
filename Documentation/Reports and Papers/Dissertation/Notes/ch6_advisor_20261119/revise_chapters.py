"""Prepare revised chapters from the verified November 19 editor readbacks."""
from pathlib import Path
from zipfile import ZipFile
import re, sys, hashlib
from difflib import unified_diff
from collections import Counter

here = Path(__file__).resolve().parent
diss = here.parents[1]
chapter_dir = diss / "Final/chapters"
with ZipFile(here / "LIVE_SOURCES.zip") as z:
    originals = {n: z.read(n).decode("utf-8") for n in z.namelist()}
def lines_of(name):
    return originals[name].split("\n")
six = lines_of("60-sns.tex")
seven = lines_of("70-overall.tex")
app = lines_of("95-AppendixD.tex")
# Line numbers belong to the captured source. Leave every graphic, label, and
# reproducibility location intact unless its change is recorded below.
def put(lines, n, new):
    assert 0 < n <= len(lines)
    lines[n-1] = new
def cut(lines, first, last):
    for n in range(first,last+1):
        put(lines,n,None)
def finish(lines):
    return re.sub(r"\n{3,}", "\n\n", "\n".join(x for x in lines if x is not None))

new6 = {
7: r"""I developed the tools in this chapter to couple neural controllers to Gait2392, a human musculoskeletal model, and to mechanical models based on the computer-aided design (CAD) of our bipedal robot. Slow simulation and hand tuning in an earlier course project with Connor Morrow made it difficult to repeat the same experiment after changing the controller or the mechanical model. A scripted version of the walker in the MuJoCo physics engine made automatic tuning possible. The subsequent experiments use the human model, a robot counterpart with relocated muscle paths, and a model built from the lab robot's own geometry. I built a graphical editor to construct neural circuits for these experiments. A simulated knee reflex tests whether a sensory signal passes through a motor neuron and muscle to produce joint motion. The final experiments use scripts to tune controllers, test whether the feet alternate under load, and measure how much a faster neural calculation saves over a complete trial.""",
11: r"""I built my first neuromechanical model in AnimatLab, a graphical program in which a neural circuit can activate muscles and receive sensory signals from the moving body \citep{cofer_animatlab_2010}. For the course project, Morrow and I connected the circuit of \citet{deng_neuromechanical_2019} to the sagittal-plane biped of \citet{li_bipedal_2017} and adjusted neural currents, synaptic weights, and muscle parameters by hand. The circuit supplied a common rhythm to both legs, and opposing muscles pulled on each limb. I added a passive toe joint. The legs cycled in suspension, but the model tipped when released toward the ground (Figure~\ref{fig:animatlab_model_architecture}). Hunt and I reported the suspended model in a poster for the NIH BRAIN Initiative Investigators Meeting \citep{bolenTorqueValuesArtificialKnee2021June15-17}. Each parameter change required another slow run, making a systematic search impractical.""",
25: r"""\caption[Original controller-transfer experiment]{The Li--Deng model cycled its legs in suspension (A) and tipped when released for the ground test (B).}""",
29: r"""The fall could arise from a late extensor command, insufficient muscle force, or a small moment arm. Moving a muscle attachment can change its moment arm and, in a circuit that senses muscle length, the stretch signal returned to the controller. Distinguishing these explanations requires recording the neural command, muscle force, joint motion, and sensory feedback on the same clock while changing one part of the model at a time.""",
33: r"""MuJoCo (Multi-Joint dynamics with Contact) calculates joint motion and ground contact forces \citep{todorov_mujoco_2012}. SNS-Toolbox, a Python library, simulates the neurons and synapses of a synthetic nervous system (SNS) \citep{nourse_sns_2023}. I connected the two so that a script could reload a controller, change a pathway, and repeat a trial with the same body, external drive, and support.""",
35: r"""I developed the Walker Two-Layer (W2L) controller from the Li--Deng walker. Its central pattern generator (CPG) is a circuit that produces rhythmic motor commands. Each leg has its own rhythm generator (RG) to set the timing of its cycle. The two sides are linked, and heel and toe contact signals return to the circuit. This arrangement made it possible to switch off selected pathways and measure whether the feet still alternated under load.""",
37: r"""Two additional mechanical models came from OpenSim, software for simulating muscles and skeletons \citep{delp_opensim_2007}. Its Gait2392 human model provides a reference with 23 degrees of freedom and 92 musculotendon actuators (Chapter~\ref{ch:background}). In the robot counterpart, I relocated muscle paths to accommodate braided pneumatic actuators (BPAs), artificial muscles that shorten under pressure (Figure~\ref{fig:simulation_workflow}).""",
39: r"""A separate model came from the lab robot's SOLIDWORKS assembly. This CAD model preserves the joint mechanisms and muscle attachments of the design. The neural circuits and actuators can also be represented in Simulink, the block-diagram simulator used with MATLAB, using the SNS-Simscape library developed here. The MathWorks bridge sends commands from Simulink to a MuJoCo model and returns motion and force measurements. Table~\ref{tab:sns_model_configurations} identifies the bodies and actuators used in each experiment.""",
69: r"""MyoConverter translates OpenSim muscle paths and joint geometry into a model for MuJoCo \citep{caggiano_myosuite_2022,wang_myosim_2022}. I used it to convert the human reference, \path{gait2392_simbody.osim}, restored missing maximum-force scales from the OpenSim source, and checked the direction of each tested muscle's joint action. The model retains 92 actuator entries, including six disabled small hip rotators, and the OpenSim convention of negative knee flexion. Moving the joints through their ranges exposed changes in moment arm at moving points and wrap surfaces. A separate check verified that foot contact produced the ground forces (Appendix~\ref{sec:app_conversion_repairs}).""",
71: r"""To check whether the converted muscles pulled as strongly as their OpenSim counterparts, I compared their forces at 121 poses from a reference walk using the same joint positions and muscle activations. The reference forces and activations came from OpenSim static optimization, which estimates how muscles share the load at each pose \citep{opensim_static_optimization_docs}. For eight of twelve right-side muscles, MuJoCo reproduced the variation in force across poses and the mean force agreed within \qty{10}{\percent}. Two muscles illustrate why both checks matter. For rectus femoris, the mean forces were close, about \qty{48}{\newton} in OpenSim and \qty{52}{\newton} in MuJoCo. At one pose, however, OpenSim estimated \qty{96}{\newton} and MuJoCo estimated \qty{416}{\newton}. For the short head of biceps femoris, the forces rose and fell in a similar pattern, but MuJoCo's mean was about \qty{85}{\percent} higher. For both muscles, using the mean OpenSim force at every pose gave a smaller squared error than using the converted model's forces. These two muscles therefore require a separate check before using their forces to interpret a controller test. Appendix~\ref{sec:muscle-models-compare} gives the calculation and all twelve results.""",
81: r"""Moving muscle attachments on the OpenSim counterpart changes where an artificial muscle pulls. Testing our robot's four-bar knee and mass distribution also requires the geometry of the mechanism we designed (Chapters~\ref{ch:background} and~\ref{ch:torque}). I assembled the mechanical model in MuJoCo from the CAD links and attachment locations. The mechanism can also be imported into Simscape Multibody, the Simulink library for rigid bodies and joints. Experiments using the MathWorks bridge require an executable MuJoCo model.""",
85: r"""In the supported CAD stepping experiment, the SNS set the timing and an engineering controller distributed reference tensions among 54 actuators. Foot forces and clearances showed whether one leg accepted load while the other unloaded and swung. In a pressure-driven experiment, the separate BPA interface calculates the available tension from commanded pressure and routed length.""",
89: r"""I adapted W2L to the Gait2392 human model. In this variant, W2Lvar, pattern-formation (PF) channels time the activity of motor neurons (MNs), the neurons that recruit muscles. One pair of channels coordinates the hip muscles on each side, and another coordinates the knee and ankle muscles. The circuit contains 888 neurons and 7,868 connected neuron pairs. To change how a knee extensor is recruited, a researcher needs to find the group of motor neurons that drives it, trace the inputs to that group, and inspect the strength and sign of each connection.

AnimatLab provides a drawing surface for wiring neurons and synapses. In the SNS-Toolbox workflow, the researcher specifies those connections in Python and can render a diagram of the resulting network \citep{cofer_animatlab_2010,nourse_sns_2023}. A diagram helps locate a pathway. An editable diagram also lets the researcher change the connection where it is drawn. I built a graphical circuit editor and a library of Simulink blocks to make selected pathways in a connectome, the map of connections among neurons, accessible in this way. A nonspiking neuron continuously transmits changes in membrane voltage through its synapses over their operating range. A spiking neuron signals by emitting a discrete spike when its voltage reaches threshold. Scripts remain useful for generating repeated connections and testing many parameter settings. Table~\ref{tab:sns_capability} identifies the tools used for each task.""",
111: r"""I implemented SNS-Simscape as Simulink blocks for neural activity, synapses, sensory signals, muscle activation, and actuator force (Figure~\ref{fig:simulink_sns_implementation}). The library includes both nonspiking equations based on SNS-Toolbox \citep{nourse_sns_2023} and a leaky integrate-and-fire (LIF) block. Sensory blocks represent group Ia fibers, which respond to muscle stretch, and group Ib fibers, which respond to tension.

In the graphical editor, a researcher can place nonspiking neurons, draw synapses between them, inspect a synapse's strength and reversal potential, and save the circuit. The Build action generates a Simulink circuit from that saved description. Appendix~\ref{sec:app_signal_contracts} gives the editor command and build script. The knee and arm examples use separate builders to add their actuator, sensory, and mechanical blocks.""",
120: r"""The Sensory Afferent Database (SADb) records which sensory nerve fibers an experiment stimulated, which motor population responded, and under what conditions (Chapter~\ref{ch:database}). These observations guide the choice of a connection's sign and destination. A pulse to an input that projects to a motor neuron tests recruitment. A connection to a pattern-formation or rhythm-generating channel can also alter the timing of the next burst.""",
122: r"""\subsection{Verifying Neural and Mechanical Coupling}\label{sec:sns_simulink_demos}""",
124: r"""Two single-joint tests verify complementary parts of SNS-Simscape. The knee test asks whether a known sensory input reaches a motor neuron and moves the joint. The arm test starts on the mechanical side and asks whether the circuit senses and resists an increasing load.""",
126: r"""In the simulated knee, a brief pulse to the extensor Ia sensory input recruited the extensor motor neuron, increased extensor force, and extended the knee. The stimulated knee reached \ang{1.09} more extension than the control \qty{111}{\milli\second} after pulse onset (Figure~\ref{fig:sns_reflex_current}). The aligned traces show the order of the response. Appendix~\ref{sec:app_knee_transient} gives the wiring, initial conditions, and recovery test.""",
135: r"""As the simulated cup filled, the arm descended and stretched the biceps. With feedback enabled, the circuit recruited the biceps, which increased the force of its pneumatic actuator. The arm sagged \ang{2.00} with nonspiking feedback and \ang{7.85} with feedback disabled at the end of the same load ramp (Figure~\ref{fig:sns_beverage_current}). Together, the knee and arm examples verify that these blocks transmit sensory input to muscle force and return the mechanical response to the circuit. Appendix~\ref{sec:spiking_results} gives the matched spiking comparison.""",
147: r"""Figure~\ref{fig:sns_design_fit_workflow}A follows two connections in the saved left-leg W2L circuit. One PF neuron excites the extensor motor pool and the other excites the flexor pool. Moving a muscle attachment in CAD can alter its moment arm. For a controller that senses muscle length, the new route can also change the feedback it receives. To compare designs, a researcher can update the attachment and its sensor and actuator mapping, then rerun the saved circuit with the same external drive and load. The force and joint-angle records provide the comparison.""",
168: r"""With the circuit and mechanical model accessible from a script, a parameter search can run a series of trials without hand tuning between runs (Figure~\ref{fig:sns_design_fit_workflow}B). I organized the searches by function. Each stage specifies what to measure and which parameters to vary, then retests the selected settings. The targets progress from rhythm to recruitment to movement under load.""",
170: r"""For a PF channel, the target specifies when it should drive its motor neurons, for how long, and how much its activity should overlap that of neighboring channels. The test measures the signal transmitted across its synapse. Voltage can vary below the synapse's operating range without changing motor recruitment. Appendix~\ref{app:neuron_equations} describes the nonspiking and spiking formulations.""",
172: r"""I adapted the four-stage training method of \citet{hunt_development_2017} to the human model. The first stage sets the intrinsic period to approximately twice the target locomotor period. The second tunes the circuit to synchronize with the sensory input. The third fits feedforward motor recruitment, and the fourth refines it with sensory feedback. The recruitment targets came from inverse mechanics applied to reference motion. Agreement with those targets measures how closely the circuit reproduces the estimated muscle activation.""",
174: r"""The search scripts bound the parameters, tie corresponding left and right values when symmetry is intended, and use tree-structured Parzen estimation to select candidates from earlier scores. They save the parameters, components of the score, and reason each trial ended. A selected candidate is rerun in a fresh process before any intervention.""",
176: r"""The physical tasks change the loads and sensory information available to a saved controller (Table~\ref{tab:sns_curriculum}). During suspension, \emph{deafferented} means that sensory feedback to the controller is disabled and \emph{afferented} means that it is enabled. Suspension permits a test of recruitment and timing before the feet meet the ground. Standing and supported stepping then test whether the legs bear load and transfer it between the feet.""",
195: r"""Progress on one task does not establish performance on the next. Accurate motor recruitment can coexist with a timing error in the circuit, and cycling the legs in suspension requires no weight transfer. The tests therefore examine PF transmission separately from motor recruitment and foot loading separately from joint motion.""",
199: r"""In a pathway-removal experiment, the unaffected controller parameters, drive, support, and external force protocol stay fixed. In a stimulation experiment, the driver changes the injected current at one location and repeats the trial from the same state. Aligning RG, PF, and MN voltages with joint motion, actuator force, and support force identifies where the responses first diverge.""",
201: r"""The timing comparison uses the same neural and mechanical trajectories as its starting point. Timing neural updates and the complete simulation separately accounts for the time spent on initialization, mechanics, and recording.""",
203: r"""\subsection{Testing Bilateral Coordination in Suspension}\label{sec:preliminary_sim_results}\label{sec:sns_biw2l_results}""",
205: r"""I extended W2Lvar to coordinate muscles that span the hip and knee. The resulting nonspiking controller, BiW2Lvar+, was calibrated using the four-stage method above. Suspension provided a first test of whether the two legs could cycle in alternation before either foot had to support weight.""",
207: r"""With the same fitted weights, the Gait2392 model cycled both legs in each suspended trial (Figure~\ref{fig:sns_biw2l_air}). The period shortened from \qty{2.466}{\second} to \qty{1.300}{\second} in the second trial. Both sensory feedback and sustained neural drive changed between trials, so this comparison cannot attribute the shorter period to feedback alone. The feet carried no load. Appendix~\ref{sec:app_biw2l_air} gives the controller settings, support forces, and joint-angle measurements.""",
218: r"""The same recording tools also show whether a change in one neuron's voltage reaches the next neuron. A pulse delivered while the hip-flexor PF pathway was inactive changed its voltage without recruiting the iliacus motor neuron. The same pulse during its active phase reached the motor neuron. Appendix~\ref{sec:app_current_noise} gives the phase-specific stimuli and repeat trials under neural noise.""",
224: r"""A knee that flexes too little gives no direct indication of where the error arose. To check the neural command, I tested syn6, a controller that distributes muscle recruitment among six pattern-formation channels in each leg of the Gait2392 human model. In the original circuit, the fifth channel continued to excite its motor neurons during the part of the cycle when it should have been silent. Changing the connection strengths improved the fit to muscle recruitment while leaving that timing error in place. Recording the signal leaving the channel showed which part of the controller needed revision.""",
226: r"""I changed how the channels started and stopped their bursts. The revised circuit met its timing targets with prescribed sensory inputs calculated from the reference movement. When the moving, suspended legs supplied their own feedback, a subsequent repair kept the second channel to one burst per cycle, but the third still produced three or four bursts. The trial completed with the knee-extension limit respected, but its rhythm missed the acceptance criteria. The phase-averaged knee curves spanned \ang{19} on the left and \ang{34} on the right, compared with about \qtyrange{70}{72}{\degree} in the walking references (Figure~\ref{fig:sns_syn6_overlay}). Individual cycles varied, so the average curves conceal some of the knee motion. The extra bursts identify a remaining timing error in the controller. Neither foot supported weight in this trial. Appendix~\ref{sec:app_pf_timing} gives the timing criteria and trial conditions.""",
230: r"""\begin{figure}[!htbp]
\centering
\includegraphics[width=\textwidth]{syn6_suspended_kinematics.pdf}
\caption[Joint motion of the suspended six-synergy model]{Hip, knee, and ankle angles for seven suspended syn6 cycles (light blue) and their mean (solid blue). The reference curves are recorded walking kinematics from the OpenSim Subject01 example (dashed) and walking at \qty{1.25}{\meter\per\second} simulated by \citet{ong_predicting_2019} (dotted). The legs share phase zero at maximum right-knee flexion. Dark bars below each leg mark stance in each reference, and gray bars mark swing. The suspended model had no foot contact. Knee flexion is negative.}
\label{fig:sns_syn6_overlay}
\end{figure}
\FloatBarrier""",
234: r"""On supported level ground, the W2L experiment asked whether the feet still alternated when selected sensory pathways were removed. Drive, assistance, and the forward perturbation were held fixed. With the controller intact, the heels struck almost half a cycle apart, a pattern called antiphase coordination. Removing force feedback through the Ib pathways alone preserved that timing. Removing heel and toe contact, pooled extensor-force inputs, and their associated spinal connections together moved heel strikes outside the alternation criterion, even though the model continued forward (Figure~\ref{fig:sns_ablation_pair}). The force records exposed a loss of coordination that forward distance alone would have missed (Appendix~\ref{sec:app_w2l_perturbations}).""",
243: r"""A later controller, W2L108, also produced alternating swings on a moving treadmill with upward assistance from rails and guides. Each qualifying swing unloaded and cleared the whole foot. A second candidate used less assistance but missed its short-shuffle target (Appendix~\ref{sec:app_w2l108_ground}).""",
247: r"""On the separate CAD robot, the neural circuit supplied timing and an engineering controller allocated tension among 54 ideal actuators. The simulated feet exchanged load and cleared the ground with vertical assistance averaging \qty{77.5}{\percent} of body weight during the scoring interval (Figure~\ref{fig:ah001_supported_control}). These tests used prescribed tensions and guidance. The pressure-dependent BPA force law is available for experiments in which actuator length and pressure determine force. Appendix~\ref{sec:app_cad_supported_stance} gives the scoring window, assistance, and repeat trials.""",
308: r"""I developed the SNS-Simscape block library and graphical editor so that a researcher can draw a nonspiking circuit, inspect an individual connection, and generate a Simulink model from the saved circuit. The actuator and sensory blocks link muscle recruitment to force and joint motion in the single-joint demonstrations. For experiments with MuJoCo mechanics, the adapted interface sends actuator commands to a prepared model and returns joint and actuator measurements. The mechanical preparation tools make the Gait2392 human model, its rerouted robot counterpart, and a model assembled from the lab robot's CAD available for controller experiments.

The calibration and recording scripts make comparisons repeatable. A future researcher can load a saved circuit, change one sensory pathway, and repeat the trial with the body, drive, and support held fixed. The recordings show how that change affects neural activity, muscle force, joint motion, and foot contact. A change in muscle routing can be tested with the controller settings held fixed. The sparse neural implementation shortens repeated trials on the tested Gait2392 controller (Table~\ref{tab:sns_runtime}). These capabilities allow researchers to use human mechanics as a reference and simulate how a change in the controller or the robot's design affects coordination.""",
}
for n,t in new6.items():
    put(six,n,t)
for n in range(45, 68):
    six[n-1] = six[n-1].replace("Antagonist Hill-type muscles", "Opposing muscle models with contractile and elastic elements (Hill-type)")
for first,last in [(12,13),(216,217),(220,221),(228,229),(245,246),(256,257),(298,305)]:
    cut(six,first,last)

for i,line in enumerate(six):
    if line and r"\label{sec:sns_literature_results}" in line:
        six[i]=re.sub(r"\\subsection\{[^}]+\}",lambda _:r"\subsection{Tracing a Timing Error in the Six-Synergy Controller}",line,count=1)
for i,line in enumerate(app):
    if line and line.startswith("Feedforward and feedback fits in the s3k and syn6 families"):
        app[i]=line.replace("the s3k and syn6 families", "s3k (the Gait2392 reference controller with four pattern-formation cells per leg) and syn6 (the six-synergy controller)")
    elif line and line.startswith("An earlier syn6 air experiment"):
        app[i]=r"""Figure~\ref{fig:sns_syn6_overlay} uses the completed suspended trial in which a shutoff latch kept the second channel silent until the next extensor rhythm signal on the same side. The motor-neuron input weights were fixed, sensory feedback remained active, and the knee-extension guard remained enabled. The second channel had one onset in each of six left and five right rhythm cycles. The third channel had three or four onsets per cycle. Both feet recorded zero vertical load. Seven complete intervals between right-knee flexion peaks, from \qtyrange{8.004}{13.652}{\second}, were sampled at 200 phase points. All joints share the right-knee phase origin. The reported \ang{19} and \ang{34} excursions are the ranges of the phase-averaged left and right knee curves. Measuring each cycle's excursion before averaging gives \ang{48.6} and \ang{49.6}, respectively. Variation among cycles reduces the excursion of the mean curves. The primary record is \path{s2_latch/capture/aligned_trace.npz}. The script \path{stage2_timing_figures.py} calculates the reference overlays and the per-channel burst counts."""+"\n\n"+line
    elif line and r"\label{sec:app_execution_provenance}" in line:
        app[i]=line+r"""

Figure~\ref{fig:app_research_toolchain} follows a saved experiment from model preparation to inspection. The circuit editor supports graphical wiring, and each experiment can be run from its corresponding script.

\begin{figure}[htbp]
\centering
\includegraphics[width=\textwidth]{research_toolchain_workflow.pdf}
\caption[Neuromechanical experimentation workflow]{The source model determines the mechanical preparation required. The researcher connects a neural circuit to actuators and sensors, executes the simulation, and compares the recorded signals across trials. The solid boxes show steps supported by the tools listed in Table~\ref{tab:app_script_map}. The dashed band identifies a proposed interface for selecting and launching experiments across all six steps.}
\label{fig:app_research_toolchain}
\end{figure}"""
    elif line and line.startswith(r"Table~\ref{tab:app_script_map} gives the executable sources."):
        app[i]=line.replace("gives the executable sources.", "lists the scripts for preparing mechanics, assembling the neural and actuator model, executing a trial, and inspecting the saved results. The graphical editor covers circuit wiring. The preparation, parameter searches, and replay have model-specific script entry points.")

# Keep the experimental scope and the numerical details with the appropriate appendix.
put(app,45,r"""The reference was the muscle force estimated by OpenSim static optimization along the subject01 walking trajectory \citep{opensim_static_optimization_docs}. For each of 121 poses, the MuJoCo model received the corresponding joint positions and OpenSim activations. Its rigid-tendon actuators then supplied a second force estimate. This framewise comparison tests posture and activation. It does not reproduce the muscle velocities, compliant-tendon history, or forward dynamics of a walking trial.""")
put(app,47,r"""For each muscle, let $F_i$ be the OpenSim force at pose $i$, $\hat F_i$ the MuJoCo force, and $\bar F$ the mean of the OpenSim forces over the sampled poses. The force-agreement score in Table~\ref{tab:force-r2} is
\[
R^2=1-\frac{\sum_i(\hat F_i-F_i)^2}{\sum_i(F_i-\bar F)^2}.
\]
The denominator is the squared error from using one constant, $\bar F$, at every pose. A score of one indicates identical force traces. A score below zero means the converted model accumulated more squared error than that constant estimate. The separate squared Pearson correlation measures how the traces vary together and can remain high when their force magnitudes differ. The mean force ratio is MuJoCo divided by OpenSim. Lag is the offset, in analysis frames, of maximum cross-correlation. The calculation and the retained results are in \nolinkurl{Code/MuJoCo_SNS/spinal/_muscle_force_compare.py} and \nolinkurl{muscle_force_compare.csv}.""")
put(app,51,r"""  \caption[Muscle force at matched poses and activations]{Right-side muscle forces from MuJoCo and OpenSim across 121 poses of the subject01 walk, \qtyrange{0.5}{2.5}{\second}. The mean ratio compares force magnitude. Squared correlation compares the variation across poses. Each lag frame is approximately \qty{17}{\milli\second}.}""")
put(app,54,r"""  \begin{tabular}{lrrrr}""")
put(app,56,r"""    Muscle & Force $R^2$ & Shape $r^2$ & Mean ratio & Lag (frames) \\""")
import csv
repo = diss.parents[2]
with (repo/"Code/MuJoCo_SNS/spinal/muscle_force_compare.csv").open(encoding="utf-8") as f:
    rows = list(csv.reader(f))[1:]
assert len(rows) == 12
for n,row in enumerate(rows,58):
    cells = app[n-1].rstrip().removesuffix(r"\\").split("&")
    assert len(cells)==4
    cells[0] = cells[0].replace("biceps femoris short ", "biceps femoris short head ")
    put(app,n,"&".join(cells[:3]+[" "+format(float(row[4]),".3f")+" ",cells[3].rstrip()+" "])+r"\\")
put(app,74,r"""Eight muscles had force-agreement scores of 0.88--0.98 and mean forces within \qty{10}{\percent} of OpenSim. The rectus femoris mean ratio was 1.088, but the force estimates followed different patterns across poses ($r^2=0.384$). Its summed squared force error was about 2.15 times the mean-reference squared error. For biceps femoris short head, the shape correlation was higher ($r^2=0.868$), but the mean ratio of 1.846 indicates that the converted actuator produced about \qty{85}{\percent} more force on average. Its squared error was about 1.26 times the mean-reference error. Both scores compare errors in force magnitude at matched poses.""")
put(app,76,r"""The paired forces for every pose are archived in \path{force_arrays.npz} with the force-comparison figure components for this appendix. At \qty{1.550}{\second}, for example, rectus femoris produced \qty{96.38}{\newton} in OpenSim and \qty{415.57}{\newton} in MuJoCo. The short head of biceps femoris averaged \qty{112.80}{\newton} in OpenSim and \qty{208.27}{\newton} in MuJoCo. The rectus femoris path was altered during conversion, and the two models use different tendon mechanics. Distinguishing path, tendon, and activation effects would require a controlled comparison at matched motion and activation histories.""")
put(app,92,app[91]+r""" From \path{Code/Matlab/SNS_Simscape/gui}, the command \verb|app = SNSGui()| opens the editor. The Build action calls \path{sns_build_from_json.m}.""")
put(app,117,app[116]+r""" The reproduction driver is \path{run_transient.m}, using the archived \texttt{KneeReflexDemo} model.""")
put(app,200,app[199].replace("Hunt supplied the calibration method.", r"The four training stages follow \citet{hunt_development_2017}."))
# Do not change the meaning of continuously varying synapses in hybrid spiking networks.
def neuron_terms(text):
    text = text.replace("Non-Spiking (Graded) Chemical Synapses","Voltage-Dependent Chemical Synapses")
    text = re.sub(r"\bgraded (chemical )?synaps", lambda m:"voltage-dependent "+(m[1] or "")+"synaps", text, flags=re.I)
    text = text.replace("Voltages are graded membrane signals", "The traces show continuously varying membrane voltages")
    text = text.replace("graded membrane signals","continuously varying membrane voltages")
    text = text.replace("non-spiking","nonspiking").replace("Non-spiking","Nonspiking").replace("Non-Spiking","Nonspiking")
    return re.sub(r"\bgraded\b", lambda m:"Nonspiking" if m[0][0].isupper() else "nonspiking",text,flags=re.I)

put(seven,61,r"""For an earlier course project, I built a neuromechanical model in AnimatLab and worked with Morrow to couple the Deng circuit to the Li walker and tune it by hand. The legs cycled while suspended, and the model tipped when released. I developed the tools in Chapter~\ref{ch:sns} to trace the consequences of a controller or mechanical change from neural activity to muscle force and joint motion. For example, an extensor sensory pulse extended the simulated knee. In a supported walker, removing selected sensory inputs changed the timing of footfalls even though the body continued to move forward. Recording both contact and movement exposed the coordination deficit.""")
put(seven,67,r"""Future researchers can test whether particular sensory inputs contribute to the extra bursts in the third pattern-formation channel during movement. The suspended six-synergy trial provides a record of the extra bursts and the joint motion at their onset. Replaying the same recorded motion while changing one sensory pathway at a time can distinguish the effect of that pathway from the accompanying change in knee angle. The layer-specific stimuli already tested in a reduced circuit can then be applied to the six-synergy controller under load. Stimulating the rhythm generator, a pattern-formation channel, or a motor neuron in separate trials can distinguish a reset of the rhythm from a missing motor burst with preserved phase (Section~\ref{sec:app_layer_stimulus}). Published responses provide a reference for those interventions \citep{selionov_tonic_2009,rybak_modelling_2006,markin_afferent_2010,shevtsova_modeling_2016}.

The force comparison in Appendix~\ref{sec:muscle-models-compare} identifies two muscles for a separate mechanical test. Future researchers can prescribe the same motion and activation history to rigid- and compliant-tendon models, then vary activation with the tendon model fixed. On the robot model, a measured pressure history can drive a characterized BPA route. Repeating a controller or contact experiment with that force law tests whether the force available at the measured length changes coordination.""")
put(seven,84,r"""Figure~\ref{fig:pipeline} shows how future researchers can test the embedded interface before connecting the biped's valves. A real-time leg simulation can supply measurements to the physical controller and receive its actuator commands. Recorded signals can first exercise the sensory encoder, the command interface, and the watchdog that detects stale commands. After the compiled neural controller is connected to the hardware inputs and outputs, replaying the same perturbation through the software and embedded controllers would expose sampling delays, missed updates, or scaling errors. Figure~\ref{fig:simulation_workflow} identifies the existing simulation interfaces, and Appendix~\ref{sec:app_embedded_interfaces} describes the hardware components.""")
# Use researcher voice for proposed work, retaining authorship of completed work.
for n,pairs in {
29:[("I propose a geometrically nonlinear finite-element model for this question.","Future researchers can test this explanation with a geometrically nonlinear finite-element model."),("I would represent","The model can represent"),("I would embed","Researchers can embed")],
39:[("I would prescribe","Future researchers can prescribe"),("I would measure","they can measure")],
47:[("I would measure","researchers can measure")],
49:[("I would measure","future researchers can measure")],
71:[("For the biped I would first test","For the biped, future researchers can first test")],
93:[("I would reconstruct","Researchers can reconstruct")]
}.items():
    for old,new in pairs:
        assert old in seven[n-1],(n,old)
        seven[n-1]=seven[n-1].replace(old,new)
new = {
    "60-sns.tex":neuron_terms(finish(six)),
    "70-overall.tex":finish(seven),
    "95-AppendixD.tex":neuron_terms(finish(app)),
    "93-AppendixB.tex":neuron_terms(originals["93-AppendixB.tex"]),
    # Exact live copy. No editorial change to the Introduction.
    "10-introduction.tex":originals["10-introduction.tex"],
}
# Ben had already made two copyedits in Overleaf. Allow just those words to
# differ locally, and preserve the file's existing byte and line-ending format.
torque_file = chapter_dir / "30-torque.tex"
before_torque = torque_file.read_bytes()
after_torque = before_torque.replace(b"based on bones scans", b"based on bone scans").replace(
    b"intersection point of the internal crossed links define the",
    b"intersection point of the internal crossed links defines the")
assert after_torque.decode("utf-8").replace("\r\n", "\n") == originals["30-torque.tex"]
if before_torque != after_torque:
    stage = torque_file.with_suffix(".tex.nov19-tmp")
    stage.write_bytes(after_torque)
    stage.replace(torque_file)
for name, text in new.items():
    old = originals[name]
    for token in [r"\\label\{([^}]+)",r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)"]:
        expected=Counter(re.findall(token,old))
        if name=="60-sns.tex":
            added="fig:sns_syn6_overlay" if "label" in token else "syn6_suspended_kinematics.pdf"
            expected[added]+=1
        if name=="95-AppendixD.tex":
            added="fig:app_research_toolchain" if "label" in token else "research_toolchain_workflow.pdf"
            expected[added]+=1
        assert expected==Counter(re.findall(token,text)),(name,token)
    assert not re.search(r'(?<!\\)[\x00-\x08\x0b\x0c\x0e-\x1f]',text)
    filename = chapter_dir / name
    prior = filename.read_bytes()
    # Preserve the existing source's line-ending convention when writing the local files.
    eol = "\r\n" if prior.count(b"\r\n") > prior.count(b"\n")/2 else "\n"
    encoded = text.replace("\n",eol).encode("utf-8")
    if prior != encoded:
        staging = filename.with_suffix(".tex.nov19-tmp")
        staging.write_bytes(encoded)
        staging.replace(filename)
diff = "".join("".join(unified_diff(originals[name].splitlines(True),txt.splitlines(True),
                        fromfile="live/"+name,tofile="revised/"+name)) for name,txt in new.items())
(here/"MANUSCRIPT.diff").write_text(diff,encoding="utf-8")
