"""Revise only the chapter sources captured from the November 19 Overleaf editor."""
from difflib import unified_diff
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

here = Path(__file__).resolve().parent
diss = here.parents[1]
chapters = diss / "Final/chapters"
with ZipFile(here/"LIVE_CHAPTER_BASELINE.zip") as z:
    baseline = {name:z.read(name).decode("utf-8") for name in z.namelist()}
files = {name:baseline[name].split("\n") for name in
         ("50-database.tex","60-sns.tex","70-overall.tex","95-AppendixD.tex","96-AppendixE.tex")}
five, six, seven, detail, appendix = files.values()
def put(lines,n,value):
    assert 0<n<=len(lines)
    lines[n-1]=value
def replace(lines,n,old,new):
    assert lines[n-1].count(old)==1,(n,old)
    put(lines,n,lines[n-1].replace(old,new))

put(seven,41,r"""\subsection{Joint Design Factors}""")
replace(seven,61,"However, but the model","However, the model")
put(seven,53,r"""The SADb links sensory pathways to the experiments in which they were observed. For example, the two pathway labels in the McCrea record lead to the same nerve stimulus producing opposite responses under different conditions (Section~\ref{sec:sadb_worked_rule}). A researcher choosing feedback for an extensor controller can recover the stimulus and locomotor state associated with each sign before assigning that connection in a model.""")
put(seven,57,r"""The monthly literature watch supplies candidates for curation (Section~\ref{sec:sadb_methods}). Linking each reported response directly to its stimulus, preparation, and locomotor phase would make it easier to compare papers that use the same pathway name. The Bubble map already distinguishes citation communities by color. Grouping papers from each community spatially would give readers another way to find related studies.

Future researchers can revisit large language models for classification as the tools improve. In our earlier evaluation, finding related papers was more successful than assigning them to the right categories \citep{bolen_using_2024}. A model could propose a pathway label and identify the passage that supports it. A researcher would read the experiment and approve, correct, or reject the label before it enters SADb. \citet{chomsky_false_promise_2023} distinguish predicting plausible text from explaining why a conclusion follows. In the McCrea example, the words \emph{Ib disynaptic excitation} alone do not tell a curator which fibers were stimulated. Comparing suggested labels with independently curated records, and measuring the time required to check them, would show whether newer models help with that work.""")
replace(seven,29,r"Festo describes a pressure-sealing rubber wall reinforced by aramid fabric \citep[p.~2]{festo_2026}.",r"Festo describes a pressure-sealing rubber wall reinforced by aramid fabric \citep[p.~2]{festo_2026}. The exploded view in Figure~\ref{fig:bpa_exploded_wall}C illustrates how the wall, reinforcement, and end fitting can be represented in the model.")
replace(seven,34,"bpa_exploded_cutaway.pdf","bpa_loading_conditions.pdf")
replace(seven,35,"End deformation and proposed BPA experiments","BPA end deformation, proposed tests, and exploded wall model")
put(seven,39,seven[38]+r"""

\begin{figure}[H]
\ContinuedFloat
\centering
\includegraphics[trim=0 0 0 20bp,clip,width=\textwidth]{bpa_exploded_cutaway.pdf}
\caption[]{\textit{Continued.} (C) Exploded longitudinal half-section illustrating the proposed BPA wall model. The fixed fitting constrains the end. The separated surfaces identify the rubber wall, two oppositely wound aramid fiber families, and the outer covering. Sectioning a specimen would establish the dimensions, fiber angles, and interfaces for a finite-element model.}
\label{fig:bpa_exploded_wall}
\end{figure}""")

replace(six,11,"The legs cycled in suspension, but the model tipped when released toward the ground","The legs made stepping-like movements while the model was suspended in air. When we lowered it to the ground and removed the support, the model fell over")
# Ben’s later local export removes the personal attribution for the passive toe.
with ZipFile(diss/"Bolen_Dissertation.zip") as z:
    later_six=z.read("chapters/60-sns.tex").decode("utf8").replace("\r\n","\n")
expected=baseline["60-sns.tex"].replace(
    "I added a passive toe joint. The legs cycled in suspension, but the model tipped when released toward the ground",
    "There is a passive toe joint. The legs cycled in a stepping-like motion while in air, but the model fell when active support was removed")
assert later_six==expected,"Reconcile the newer chapter export before writing."
replace(six,11,"I added a passive toe joint.","There is a passive toe joint.")
put(six,23,r"""\caption[Original controller-transfer experiment]{The legs of the Li--Deng model made stepping-like movements while the model was suspended in air (A). The model fell over when it was lowered to the ground and its support was removed (B).}""")
replace(six,54,"suspended cycling","stepping-like motion in suspension")
replace(six,57,"suspended cycling example","suspended leg-motion example")
replace(six,73,"the robot model cycled its legs in suspension","the legs of the suspended robot model made repeated stepping-like movements")
put(six,112,six[111]+r"""

For larger controllers, I also developed a connectome block editor that groups neurons and synapses into nested diagrams. A researcher can open one motor pathway, inspect the cells and connections inside it, and return to the diagram of the whole circuit. A saved drawing of the Hunt controller exposes the rectus femoris and hamstring pathways on each side in this way. The drawing records how the cells are connected. Running that controller on a body also requires its neuron equations, actuators, and sensory inputs. Appendix~\ref{sec:app_execution_provenance} identifies the editor and saved example.""")
put(six,175,r"""For rhythm and sensory timing, the search scripts use tree-structured Parzen estimation to propose new parameter sets from earlier scores. For the s3k and syn6 controllers, a separate bounded least-squares fit adjusts the connections from pattern-formation channels to motor neurons to match the reference muscle recruitment. Corresponding weights are tied across sides when the fit calls for symmetry. The physical test then loads the saved settings, records the components of its score and the reason the trial ended, and repeats a selected trial in a fresh process before an intervention.""")
replace(six,196,"cycling the legs in suspension requires no weight transfer","moving the legs in suspension places no load on the feet")
put(six,200,r"""In a pathway-removal experiment, the unaffected controller parameters, drive, support, and external force stay fixed. For a stimulation experiment, the driver changes the injected current at one location and repeats the trial from the same state. For the Gait2392 trials, I developed a recorder that saves the timing of neural inputs and updates alongside joint states, muscle commands, foot loads, and harness forces. A researcher can align the signals at the stimulus or foot-contact event and locate the first response that differs. The saved joint states also allow the movement to be viewed again without running the controller a second time.""")
replace(six,206,"could cycle in alternation","could move in alternation")
replace(six,208,"the Gait2392 model cycled both legs in each suspended trial","both legs of the suspended Gait2392 model made repeated stepping-like movements in each trial")
replace(six,240,"Joint cycling continued","The joints continued to flex and extend")

replace(detail,185,"The legs cycled while the walker was held up and the model fell on release in the demonstration shown in", "The legs made stepping-like movements while the model was suspended in air. The model fell over when it was lowered to the ground and its support was removed, as shown in")
replace(detail,351,"continued joint cycling","continuing joint flexion and extension")
replace(detail,384,"the untuned spiking version cycled both legs in the supported ground trial, while the nonspiking winner cycled predominantly on the right","the untuned spiking version produced repeated movements in both legs during the supported ground trial, while movement in the nonspiking winner occurred predominantly on the right")
for i,line in enumerate(detail):
    if line.startswith(r"\caption[Matched neural and mechanical responses with nonspiking and spiking controllers]"):
        detail[i]=r"""\caption[Matched neural and mechanical responses with nonspiking and spiking controllers]{Both controllers used the same parameters and support. Columns show nonspiking (left) and untuned spiking (right) responses during the \qtyrange{10}{12.5}{\second} ground-trial window. Solid orange and dashed blue distinguish extensors and flexors in the right-side RG, PF, and muscle traces (vastus lateralis and semimembranosus). For hip, knee, and ankle angles they distinguish right and left. Spiking readouts are filtered. Muscle activation is the actuator state.}"""
        break
else:
    raise ValueError("Missing matched-ground figure caption")
for i,line in enumerate(detail):
    if line==r"\subsection{Simulink Demonstrations}":
        detail[i]=r"\FloatBarrier"+"\n"+line
        break
else:
    raise ValueError("Missing start of the Simulink demonstrations")
put(detail,457,detail[456]+r"""

To inspect the Hunt circuit by module, open \path{Hunt_fourcell_session1_editor.html} in the repository's \path{Neuromechanical_Models/Mujoco_SNS_models/} directory. The saved diagram groups enabled connections into 19 modules. Opening a rectus femoris or hamstring module shows its cells and local wiring. The drawing can be edited and saved. For a supported W2L108 project, the desktop editor can also write a synaptic edit to working copies of the native AnimatLab project and simulation files. Its launcher is \path{Code/MuJoCo_SNS/connectome_app/Start_Connectome_Editor.cmd}.

For a saved Hunt body trial, the aligned recorder writes neural and mechanical timestamps, joint states, muscle commands, foot loads, and the harness wrench, together with the compiled model used for replay. A separate mechanics file contains the individual actuator forces. The plotting tools read each signal on its own time base. To render a different view, the replay restores the saved positions and recomputes the geometry for each frame. The plotted signals and rendered movement therefore describe the same trial.""")
put(detail,477,detail[476]+r"""
\path{connectome_block_editor.html} (P) & Inspect or edit a nested circuit drawing and save the diagram. The Hunt example opens with its modules already connected. \\""")
put(detail,480,detail[479]+r"""
\path{controller_profiles_20261104/train.py}; \path{fit_waveform_edges.py} (P) & Search bounded rhythm and timing parameters, or fit bilateral pattern-formation weights to motor recruitment by least squares. \\
\path{hunt_training/verification_capture.py}; \path{verification_figures.py} (P) & Record signals and sampling times, plot selected intervals, and replay saved body positions with the corresponding model. \\""")

put(five,4,r"""Choosing a sensory connection for a locomotor controller requires more than the name of a receptor. A researcher needs to know which nerve was stimulated, which motor neurons responded, and whether the animal was at rest or moving. Those details are scattered across studies that used different species, preparations, and tasks. I developed the Sensory Afferent Database (SADb) to keep a proposed feedback pathway linked to the experiment that supports it. The collection also links to computational models so that a researcher can compare an experimental finding with the assumptions used to represent it. This work forms part of the AARL's NeuroNex Communication, Coordination, and Control in Nervous Systems (C3NS) effort.""")
put(five,8,"")
put(five,10,r"""SADb Explorer is the browser application for the collection. To use the distributed copy, open \path{SADb_audit/dist/sadb_app.html} from the project repository in a browser. The Table tab retrieves individual records, Pivot counts papers that share classifications, and Bubble map follows citation links (Figure~\ref{fig:sadb_clusters}A). These three views work offline. Search (online) queries the OpenAlex bibliographic index for prospective additions. The researcher can save candidate records for later curation.""")
put(five,19,r"""I curated the database in Airtable, an online database system. Linked tables connect papers, pathways, and models. Exporting the records with those links preserves the route from a pathway name to the relevant paper, its curation note, and the original study in the standalone Explorer.""")
put(five,32,r"""\usetikzlibrary{arrows.meta}
\begin{tikzpicture}[every node/.style={font=\small},curationstep/.style={draw,rounded corners,align=center,text width=3.3cm,minimum height=0.9cm},flow/.style={-{Latex[length=3.4mm,width=2.6mm]},line width=1.1pt}]""")
put(five,39,five[38].replace("[->,line width=0.6pt]","[flow]"))
put(five,41,r"""\caption{Curation of SADb records. The monthly search supplies candidate papers. A curator reads each source, decides whether it belongs in the collection, and records the classifications and supporting evidence.}""")
put(five,47,r"""The searchable Table contains 1,393 entries across five record types: 943 Papers, 247 Models, 157 Review Papers, 42 Feedback pathways, and 4 Recorders. It includes archived entries. Pivot and Bubble map draw from the 943 Papers. Appendix~\ref{app:sadb_records} explains which records were retained in the working database and how to interpret the coverage counts.""")
replace(five,49,"Within the Papers category, 730 papers are connected by 11,130 directed citation pairs.","Of these Papers, 730 cite at least one other Paper in the collection. Counting both citing and cited papers gives 816 distinct Papers connected by 11,130 directed citation pairs.")
put(five,51,r"""A title search in Table opens a paper's record. Figure~\ref{fig:sadb_explorer} shows the search for the study by \citet{deng_biomechanical_2022}. Opening its row reveals the article DOI, the linked model, and a note about the effect of biomechanical feedback on the simulated rhythm. A researcher looking for physiological inputs can instead filter the same table by pathway, afferent, or study animal.""")
put(five,55,r"""\includegraphics[width=\textwidth]{sadb_table_record_guide.pdf}""")
put(five,56,"")
put(five,57,"")
put(five,58,r"""\caption[Finding a paper and its model in SADb Explorer]{A title search retrieves the paper by \citet{deng_biomechanical_2022}. Selected portions of its Table view and record drawer show the search, the DOI that opens the paper, and the link to its model record. The curation note explains why the study matters for controller design.}""")
replace(five,64,r"selector displays papers up to one, two, or three citation links from the focus paper, following either its references or papers that cite it (Figure~\ref{fig:sadb_clusters}, callout~(E)).","selector in callout~(B) displays papers up to one, two, or three citation links from the focus paper. Callout~(E) marks the resulting neighborhood, which includes both references and citing papers.")
replace(five,64,"and 555 papers lie within two links.","and the two-link display contains 555 papers including the focus paper itself.")
put(five,68,r"""\subsection{From a Pathway Record to a Controller Rule}\label{sec:sadb_worked_rule}""")
put(five,69,r"""In Table, selecting \emph{Ib disynaptic excitation} in the \emph{Pathway} filter returns six Papers. The record for \citet{mccrea_disynaptic_1995} carries both Ia and Ib tags and links to two pathways, \emph{Ib disynaptic excitation} and \emph{Ib disynaptic inhibition}. The paper explains why. In decerebrate cats at rest, stimulating the plantaris nerve at group~I strength inhibited medial gastrocnemius motor neurons. With the same stimulus during extension in midbrain-evoked fictive locomotion, the recorded response was excitatory. In this preparation the neural rhythm was present while the animal was immobilized. Nerve stimulation recruited group~I fibers and did not isolate Ib fibers.""")
put(five,70,r"""Figure~\ref{fig:sadb_group_i_rule} shows one rule a modeler could take from this record. During simulated rest, a brief input from the plantaris group~I pathway would reduce activity in the synergist motor pool. During the active extensor phase, the same input would increase its activity. The measurements constrain the sign and timing of the net response. Choosing the interneurons, synaptic strengths, and delays that produce it is a separate modeling decision. Appendix~\ref{app:sadb_records} shows how the result changes in a different preparation.

\begin{figure}[!htbp]
\centering
\includegraphics[width=\textwidth]{sadb_group_i_rule.pdf}
\caption[An experimental constraint on a sensory feedback rule]{The plantaris nerve stimulus and recorded medial gastrocnemius motor pool were the same in both conditions. Stimulation inhibited the motor neurons at rest and excited them during the extensor phase of midbrain-evoked fictive locomotion \citep{mccrea_disynaptic_1995}. The two signed connections represent a candidate model's net response through spinal circuitry.}
\label{fig:sadb_group_i_rule}
\end{figure}""")
put(five,72,r"""\subsection{Comparing Groups of Papers}

Pivot answers a different question. With the Table filters cleared, choose \emph{Afferent types} and \emph{Animal}, then select the Ib--Cat cell (Figure~\ref{fig:sadb_pivot}). Table opens the 21 Papers counted in that cell. \emph{Back to pivot} restores the comparison. A paper can carry more than one afferent or animal tag and appear in several cells. A cell count consequently describes how much material a researcher can inspect for a particular combination.""")
replace(five,83,r"The controller of Chapter~\ref{ch:sns} implements the mammalian subset.",r"These records let a researcher identify candidate pathways across species. Chapter~\ref{ch:sns} tests selected mammalian feedback connections in specific controllers.")
put(five,94,r"""For a modeler, the value of a curated record lies in the distinctions it preserves. The plantaris experiment yielded a different sign at rest and during extension. A controller with a single positive gain for every group~I input would miss that difference. A modeler can return to the source to choose which preparation the controller should represent and what stimulus would test it.""")
put(five,96,r"""Table~\ref{tab:sadb_afferents} in Appendix~\ref{app:sadb_records} relates each receptor to a possible model input. Table~\ref{tab:sadb_pathways} lists the pathway records and the papers linked to them. For example, \emph{Ia monosynaptic excitation} concerns the motor pool of the muscle containing the spindle, while \emph{Ia reciprocal inhibition} concerns the antagonist motor pool and an intervening inhibitory neuron. Selecting the pathway as well as the afferent specifies which connection the researcher wants to examine.""")
put(five,104,r"""SADb supplies the source, preparation, stimulus, and measured response behind a candidate rule. The rule can then be tested by holding the mechanical model fixed, changing the sensory connection, and repeating the same stimulus. Chapter~\ref{ch:sns} develops the tools for that comparison.""")
put(appendix,4,r"""Chapter~\ref{ch:database} follows a researcher from one pathway search to a candidate rule. This appendix explains which records each view counts, why two labels can apply to the same paper, and which sources support the named pathways.""")
put(appendix,14,r"""For a single paper, use Table to check the record type and archive status. For a group of papers, start in Pivot and open a cell to read the underlying records. Tagged groups overlap. For example, a Paper tagged with both Ia and Ib can appear in both afferent rows for the same animal. Adding those rows would count that Paper twice. The citation neighborhood in Figure~\ref{fig:sadb_clusters} answers a separate question about which papers cite one another.""")
put(appendix,20,r"""Table~\ref{tab:sadb_afferents} relates each receptor to a candidate model input, and Table~\ref{tab:sadb_pathways} gives the searchable pathway names. The two pathway links in the McCrea record illustrate why the preparation belongs alongside the label (Section~\ref{sec:sadb_worked_rule}). In the acute spinal preparation treated with clonidine, the inhibition seen at rest diminished during fictive locomotion, but the short-latency excitation seen in the midbrain-evoked preparation did not appear \citep{mccrea_disynaptic_1995}. A model of the spinal preparation would therefore need a different test from the one drawn in Figure~\ref{fig:sadb_group_i_rule}. The Ib pathway label is a way to retrieve this paper. Its Ia and Ib afferent tags preserve the broader identity of the stimulated group~I fibers.""")

replace(six,7,"I built a graphical editor to construct neural circuits for these experiments.","I developed Simulink blocks and a graphical interface for building nonspiking circuits.")
replace(six,89,"A diagram helps locate a pathway. An editable diagram also lets the researcher change the connection where it is drawn. I built a graphical circuit editor and a library of Simulink blocks to make selected pathways in a connectome, the map of connections among neurons, accessible in this way.","To change one connection in a large network, a researcher first has to locate the two neurons and identify the synapse between them. I developed a graphical interface and a library of Simulink blocks so that a researcher can wire a nonspiking circuit on screen and build the corresponding Simulink model.")
put(six,108,r"\subsection{Building a Nonspiking Circuit in Simulink}")
put(six,112,r"""In the SNS-Simscape interface, \texttt{SNSGui}, a researcher can place nonspiking neurons, connect them with synapses, inspect a synapse's strength and reversal potential, and save the circuit. Choosing \emph{Build} generates the neural portion of a Simulink model. Appendix~\ref{sec:app_signal_contracts} gives the command and build script. Separate builders connect the knee and arm circuits to their actuators, sensory inputs, and mechanical models.

For larger controllers, I used a separate drawing tool to organize a connectome, the map of connections among neurons, into smaller diagrams. In the Hunt drawing, for example, a researcher can open the rectus femoris and hamstring pathways on either side and inspect their wiring. Changes made there are saved in the drawing. To change the controller used in a simulation, the researcher must change its model or configuration separately.""")
replace(six,297,"SNS-Simscape block library and graphical editor","SNS-Simscape block library and its Simulink circuit interface")
replace(detail,96,"The graphical editor's Build action","The SNS-Simscape interface's Build action")
replace(detail,442,"The circuit editor supports graphical wiring","The SNS-Simscape interface supports graphical wiring of a Simulink neural circuit")
put(detail,457,detail[456].split("\n\n")[0].replace("The graphical editor covers circuit wiring.","The SNS-Simscape interface builds the neural portion of a Simulink model.")+r"""

To inspect the Hunt connectome, open \path{Hunt_fourcell_session1_editor.html} from \path{Neuromechanical_Models/Mujoco_SNS_models/} in the repository. The drawing groups its wiring into 19 modules. Open a rectus femoris or hamstring module to see its cells and connections. Redrawing a connection changes the diagram. The simulation uses its own controller definition.

For a saved Hunt body trial, the recorder writes input and neural-update times, joint states, muscle commands, foot loads, and harness forces, together with the compiled model used for replay. A separate mechanics file contains the individual actuator forces. The plotting tools read each signal on its own time base. To render a different view of the same trial, the replay restores the saved positions and recomputes the geometry for each frame.""")
put(detail,477,detail[476].replace("Inspect or edit a nested circuit drawing and save the diagram. The Hunt example opens with its modules already connected.","Draw, inspect, and save a map of connections among neurons. The Hunt drawing groups connected cells by pathway."))
replace(five,47,"The searchable Table contains 1,393 entries across five record types: 943 Papers", "The searchable Table contains 1,393 entries. These comprise 943 Papers")
replace(five,49,"The Louvain algorithm identifies citation communities: groups of papers with comparatively dense citation connections.", "The Louvain algorithm groups papers that frequently cite one another into citation communities.")
replace(five,51,"A title search in Table opens a paper's record.","In Table, a researcher can search by title and select a row to open the paper's record.")
replace(five,51,"the linked model","the name of the associated model")
replace(five,58,"the link to its model record","the model named in the database")
replace(five,62,"citation count on the vertical axis","citation count on a logarithmic vertical axis")
put(five,76,r"\includegraphics[width=\textwidth]{sadb_pivot_guide.pdf}")
put(five,77,r"\caption[The SADb Explorer pivot view]{The Ib row and Cat column intersect at 21. Selecting that cell opens the 21 matching Papers in Table. The return control restores the pivot.}")
put(five,53,r"\begin{figure}[!htbp]")
put(five,62,"\n".join([five[61], "", five[63], "", five[65]]))
map_text=five[61]
for index in (62,64,66):
    put(five,index,"")
put(five,79,five[78]+"\n\n"+r"\subsection{Following Citations in Bubble Map}"+"\n\n"+map_text)
replace(appendix,23,r"References link to corpus records; Ia, II, and Ib signal summaries follow",r"The cited studies give examples of these pathways. Ia, II, and Ib signal summaries follow")
replace(appendix,28,"Recorded pathway roles (anchors)","Examples of pathway roles")
replace(appendix,33,"metabolic and nociceptive signals, converging with cutaneous and group II--III afferents as flexor-reflex afferents & phase-dependent resetting", "metabolic and nociceptive signals & convergence with cutaneous and group II--III flexor-reflex afferents and phase-dependent resetting")

replace(five,10,"SADb_audit/dist/sadb_app.html","SADb_audit/app/sadb_app.html")
put(five,51,r"""The Table tab retrieves both a paper and the evidence recorded for it. For example, choosing \emph{Ib disynaptic excitation} in the \emph{Pathway} filter returns six Papers (Figure~\ref{fig:sadb_explorer}). Select the row for \citet{mccrea_disynaptic_1995} to open its record. The drawer shows the source DOI, the Cat classification, both Ia and Ib afferent tags, and two pathway labels, \emph{Ib disynaptic excitation} and \emph{Ib disynaptic inhibition}. Reading the experiment explains why both labels apply.""")
put(five,58,r"""\caption[Finding the experimental evidence for a pathway]{Selected portions of the SADb Explorer show (A) the pathway filter, (B) the six-result search with the McCrea paper selected, and (C) its record, including the source DOI and afferent and pathway classifications. The DOI opens the experiment that supports the labels.}""".replace("with the McCrea paper selected","and the first author of the matching paper"))
put(five,69,five[68][five[68].index("In decerebrate cats"):])
put(five,70,"\n"+five[69])
replace(appendix,12,"426 have preparation tags","426 have an animal or human study tag")
replace(detail,457,"Open a rectus femoris or hamstring module to see its cells and connections.","Open the module labeled RF HAM four-cell PF on either the left or right side to see the rectus femoris and hamstring cells and their connections.")
replace(six,299,"A future researcher can load a saved circuit, change one sensory pathway, and repeat the trial", "A researcher can load a saved controller configuration in its simulation script, change one sensory pathway, and repeat the trial")

replace(five,49,"The Louvain algorithm groups papers that frequently cite one another into citation communities.","The Louvain algorithm identifies groups of papers with relatively dense citation links.")
put(five,58,r"""\caption[Finding the experimental evidence for a pathway]{The SADb Explorer pathway filter (A) returns six Papers. Selecting the McCrea row in the results (B) opens the record (C), which gives the source DOI, the afferent classes, and the two pathway labels. The DOI opens the experiment that supports the classifications.}""")
put(five,77,r"""\caption[The SADb Explorer pivot view]{With afferent type and animal selected (A), the Ib row and Cat column intersect at 21 (B). Selecting the cell opens the 21 matching Papers in Table. The return control restores the pivot (C).}""")
closing=next(i for i,line in enumerate(seven) if line.startswith("The four contributions make it possible"))
seven[closing]=r"""Figure~\ref{fig:discussion_neuromechanical_workflow} shows two ways to investigate a change in coordination. To test a sensory connection, a researcher can repeat a trial with the same body, drive, and support and compare the timing of muscle recruitment and foot contact. To test an actuator route, the researcher can keep the neural circuit and test conditions fixed, then compare the force and motion after the attachment is moved. In either experiment, aligning the recordings at the stimulus or foot contact makes it possible to identify which signal changed first.

\begin{figure}[!htbp]
\centering
\includegraphics[width=5.5in]{neuromechanical_experiment_workflow.pdf}
\caption[Comparing controller and mechanical changes in simulation]{The reference trial establishes the response under a specified stimulus, drive, and support. In each repeat, one change tests either the controller or the mechanical model. Neural activity, muscle force, joint motion, and foot load show where the responses first diverge. Appendix~\ref{sec:app_execution_provenance} identifies the model-specific tools used to prepare, run, and inspect a trial.}
\label{fig:discussion_neuromechanical_workflow}
\end{figure}

"""+seven[closing]
replace(six,85,"Designing and Testing Neural Controllers Graphically","SNS-Simscape Block Programming for Neural Controllers")
replace(six,89,"To change one connection in a large network, a researcher first has to locate the two neurons and identify the synapse between them. I developed a graphical interface and a library of Simulink blocks so that a researcher can wire a nonspiking circuit on screen and build the corresponding Simulink model.","I developed the SNS-Simscape block library and its graphical interface so a researcher can place nonspiking neurons, set synaptic parameters where the blocks connect, and build the corresponding circuit in Simulink.")
put(six,112,six[111].split("\n\n")[0])
put(detail,477,detail[476].split("\n")[0])
parts=detail[456].split("\n\n")
assert len(parts)==3
parts[1]=r"""The saved Hunt connectome drawing in \path{Neuromechanical_Models/Mujoco_SNS_models/Hunt_fourcell_session1_editor.html} shows the left and right rectus femoris and hamstring circuits in the modules labeled RF HAM four-cell PF. Changes to the drawing are saved separately from the controller used in a simulation."""
put(detail,457,"\n\n".join(parts))
put(seven,69,seven[68]+r"""

Future researchers can develop an application that links a connectome drawing to the controller used in a simulation. A researcher could select a sensory connection in the drawing, change its strength in the controller, and rerun a saved experiment with the same mechanical model, stimulus, and support. Reproducing the result of the same change made directly in the controller would test that the application updated the intended synapse.""")
def once(lines,old,new):
    matches=[i for i,line in enumerate(lines,1) if old in line]
    assert len(matches)==1,(old,matches)
    replace(lines,matches[0],old,new)

once(six,r"the 1,025-cell \texttt{hunt} source configuration","an earlier 1,025-cell configuration of BiW2Lvar+")
once(six,"In the second timing test, Hunt-PT+RC used a Peterka--Ting-inspired balance pathway (PT) that responds to motion of the model's center of mass. Its Renshaw cells (RC)", "The second timing test used a 1,033-cell configuration with a Peterka--Ting-inspired balance pathway that responds to motion of the model's center of mass. Its Renshaw cells")
once(six,r"The first two rows use the 1,025-cell source-setting \texttt{hunt} implementation, SNS-Toolbox NumPy reference versus sparse execution, in one \qty{16}{\second} stage-1 protocol pair. The last two use the separate 1,033-cell implementation labelled Hunt-PT+RC, sparse versus compact execution.","The first two rows compare reference and sparse execution of a 1,025-cell BiW2Lvar+ configuration in a single paired "+r"\qty{16}{\second} stage-1 protocol. The last two compare sparse and compact execution of the 1,033-cell configuration with balance feedback and Renshaw cells.")
once(six,"For the 1,033-cell Hunt-PT+RC controller","For the 1,033-cell configuration with balance feedback and Renshaw cells")
once(detail,r"The software selector \texttt{hunt} identifies this implementation in the records. ","")
once(detail,r"the 1,025-cell source configuration selected by \texttt{hunt}","the earlier 1,025-cell BiW2Lvar+ configuration")
once(detail,"the 1,033-cell Hunt-PT+RC circuit","the 1,033-cell circuit with the balance pathway and recurrent Renshaw inhibition")
parts=detail[456].split("\n\n")
parts[1]=r"""A saved connectome drawing in \path{Neuromechanical_Models/Mujoco_SNS_models/} shows the four-cell rectus femoris and hamstring circuit on each side of an earlier BiW2Lvar+ configuration. Changes to the drawing are saved separately from the controller used in a simulation."""
parts[2]=parts[2].replace("For a saved Hunt body trial","For the recorded Gait2392 body trial")
put(detail,457,"\n\n".join(parts))
put(five,19,r"\begin{samepage}"+"\n"+five[18]+r"\par"+"\n"+r"\end{samepage}")
put(five,53,r"\begin{figure}[H]")
put(five,70,five[69]+"\n"+r"\FloatBarrier")
replace(five,15,"SADb Explorer with Grillner (1985) selected.","Citation map of the 943 Papers in an earlier Explorer snapshot, with Grillner (1985) selected.")
put(five,10,r"""\begin{samepage}
Open SADb Explorer in a browser using \path{SADb_audit/app/sadb_app.html} in the project repository. Three tabs work offline (Figure~\ref{fig:sadb_clusters}A). Table retrieves individual records, Pivot counts papers that share classifications, and Bubble map follows citations. Search (online) queries the OpenAlex bibliographic index for candidates to save for later curation.\par
\end{samepage}""")
seven[closing]=r"""Figure~\ref{fig:discussion_neuromechanical_workflow} shows how to locate the first change in coordination. To test a sensory pathway, repeat the stimulus with the body, drive, and support held fixed. To test an actuator route, keep the controller fixed. A force shortfall at a known pressure points to the actuator or its route. If force is adequate but the foot accepts load late, compare when the muscle was recruited and when the joint moved.

\begin{figure}[!htbp]
\centering
\includegraphics[width=5.5in]{neuromechanical_experiment_workflow.pdf}
\caption[Comparing controller and mechanical changes in simulation]{Each repeat changes one sensory pathway or one actuator route. Aligning the responses at the stimulus or foot contact reveals which signal diverges first. Appendix~\ref{sec:app_execution_provenance} identifies the tools used to prepare, run, and inspect each trial.}
\label{fig:discussion_neuromechanical_workflow}
\end{figure}"""
# Use the dissertation class's normal break settings. Grouping a whole timing paragraph
# was unnecessary, and neither the layout nor the science requires a chapter-wide override.
for number in (3,4,5):
    put(six,number,"")
for i,paragraph in enumerate(six):
    if paragraph.startswith(r"{\interlinepenalty=10000 "):
        assert paragraph.endswith(r"\par}")
        six[i]=paragraph.removeprefix(r"{\interlinepenalty=10000 ").removesuffix(r"\par}")
    if paragraph.startswith((r"\widowpenalty=\snsSaved",r"\clubpenalty=\snsSaved")):
        six[i]=""
for i,paragraph in enumerate(five):
    if paragraph.strip()==r"\widowpenalty=10000":
        five[i]=""
previous=json.loads((here/"SOURCE_HASHES.json").read_text(encoding="utf-8")) if (here/"SOURCE_HASHES.json").is_file() else {}
new_hash={}
diff=[]
for name,lines in files.items():
    source=baseline[name]
    target=chapters/name
    before=target.read_bytes()
    eol="\r\n" if before.count(b"\r\n")>before.count(b"\n")/2 else "\n"
    new="\n".join(lines)
    while "\n\n\n" in new:
        new=new.replace("\n\n\n","\n\n")
    if name=="60-sns.tex":
        new=new.rstrip("\n")+"\n"
    local=before.decode("utf-8").replace("\r\n","\n")
    matching=(local==source or hashlib.sha256(before).hexdigest()==previous.get(name) or local==new)
    if name=="70-overall.tex":
        from capture_baseline import old,author
        matching|=local.replace(old,author)==source
    assert matching, "The working source has a concurrent edit: "+name
    encoded=new.replace("\n",eol).encode("utf-8")
    if encoded!=before:
        stage=target.with_suffix(".nov19-review-tmp")
        stage.write_bytes(encoded)
        stage.replace(target)
    new_hash[name]=hashlib.sha256(encoded).hexdigest()
    diff.extend(unified_diff(source.splitlines(keepends=True),new.splitlines(keepends=True),fromfile="live/"+name,tofile="revised/"+name))
(here/"SOURCE_HASHES.json").write_text(json.dumps(new_hash,indent=2)+"\n",encoding="utf-8")
(here/"MANUSCRIPT.diff").write_text("".join(diff),encoding="utf-8")
bib=diss/"Final/thesis.bib"
entry=(here/"Chomsky_2023.bib").read_bytes()
old_bib=bib.read_bytes()
if b"chomsky_false_promise_2023" not in old_bib:
    assert hashlib.sha256(old_bib).hexdigest()==json.loads((here/"LIVE_CAPTURE.json").read_text())["previous_local_bibliography_sha256"]
    if old_bib.endswith(b"\r\n"):
        entry=entry.replace(b"\r\n",b"\n").replace(b"\n",b"\r\n")
    stage=bib.with_suffix(".nov19-review-tmp")
    stage.write_bytes(old_bib+b"\n"+entry)
    stage.replace(bib)
assert bib.read_bytes().count(b"@misc{chomsky_false_promise_2023,")==1
