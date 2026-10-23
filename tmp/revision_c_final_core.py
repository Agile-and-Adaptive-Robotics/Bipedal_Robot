from pathlib import Path
import re,shutil,hashlib,json,fitz
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';S=C/'review/final_sweep_20261023'
def paragraph(s,start,new):
    a=s.index(start);b=s.find('\n\n',a)
    if b<0:b=len(s)
    return s[:a]+new+s[b:]
def getblock(s,kind,label):
    for m in re.finditer(r'\\begin\{'+kind+r'\}.*?\\end\{'+kind+r'\}',s,re.S):
        if r'\label{'+label+'}' in m.group():return m.group()
    raise ValueError(label)
def cap(block,text,short=None):
    a=block.index(r'\caption');b=block.index('{',a);depth=1;e=b+1
    while depth:
        if block[e]=='{':depth+=1
        if block[e]=='}':depth-=1
        e+=1
    return block[:a]+r'\caption'+('['+short+']' if short else '')+'{'+text+'}'+block[e:]
p=C/'chapters/60-sns.tex';s=p.read_text()
s=paragraph(s,'Air stepping did not establish',r'''The missing step was a tested link between spinal neural circuits and the human-scale and BPA-actuated mechanical models needed to study this robot. Air stepping in the combined AnimatLab model did not establish how the same controller would behave with human musculotendon mechanics, robot BPA routing, or ground-contact feedback. I therefore needed models and tests that could separate the effects of the circuit, actuator mechanics, and numerical implementation.''')
s=s.replace(r'\section{Requirements for the New Tools}',r'\section{Knowledge Gaps and the New Tools}')
s=paragraph(s,'The existing workflow lacked four connections',r'''Within this project, five model combinations were missing: a Gait2392 human model coupled to a CPG in MuJoCo--SNS-Toolbox; its BPA-routed robot counterpart coupled to a CPG; the combined Li walker and Deng controller with foot-contact feedback; graphical neural construction coupled to a maintained mechanical simulator; and the SolidWorks robot assembly coupled to MATLAB BPA models and a neural controller. These missing combinations prevented a controlled comparison of the same neural organization across biological and artificial actuation. I developed the conversion, coupling, graphical-library, and controller-port tools below to build and test those connections. The capability table identifies the connections demonstrated and those still incomplete.''')
s=paragraph(s,'The software reviewed in Chapter',r'''Existing software supplied separate parts of the solution. OpenSim supplied musculoskeletal geometry, MyoConverter translated it to MuJoCo, and SNS-Toolbox supplied neuron and synapse models \citep{delp_opensim_2007, caggiano_myosuite_2022, todorov_mujoco_2012, nourse_sns_2023}. The MuJoCo version used here supports scripted simulation and parameter searches but uses rigid tendons; its Python interface does not provide the block-based circuit assembly needed for this workflow. Simulink supplies that graphical construction, while Simscape Multibody imports mechanical geometry, including SolidWorks-derived mechanisms. AnimatLab retains the reference controllers but is no longer maintained. The tools connect these complementary environments and use matched tests to check what transfers between them.''')
s=s.replace('Converted musculoskeletal model, BPA force interface, and afferent feedback','Human and BPA-routed robot models, actuator interface, and afferent feedback')
s=s.replace('Comparable neural implementations and repeatable evaluations','Ports of reference circuits and matched-input evaluations')
s=s.replace('Rhythm-core agreement across NumPy and Simulink; agreement is specific to the tested network and inputs','Rhythm-core agreement across NumPy and Simulink; complete Li--Deng contact-feedback walking remains unverified')
s=s.replace("It extends the laboratory's earlier modeling approach \\citep{hunt_development_2017, young_analyzing_2019} to the BPA interfaces used here. The added components connect neural signals to the actuator and mechanical interfaces.",r"Compared with the earlier SNS models of \citet{hunt_development_2017} and \citet{young_analyzing_2019}, the added elements concern BPA actuation, afferent encoding, and mechanical interfaces. The library retains the conductance-based neuron and synapse formulation; it does not introduce a new neuron model. Its verification compares selected implementations and circuits, rather than benchmarking all of the earlier libraries.")
p.write_text(s,encoding='utf-8')
p=C/'chapters/50-database.tex';s=p.read_text()
first=s.index('\n\nThe primary deliverable')
s=s[:first]+r'''

The knowledge gap was the lack of a structured account, available to this project, of which measured afferent properties and feedback pathways could support the robot's controller. Relevant findings were dispersed across species, preparations, tasks, and primary papers. Reviews organized parts of this evidence \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}, but selecting a feedback rule still required recovering its experimental conditions and separating measurements from modeling assumptions. I built SADb to connect those findings to their sources, afferent classes, and pathway roles so that controller choices could be traced to the evidence and remaining unknowns could be identified.'''+s[first:]
e=C/'chapters/96-AppendixE.tex';a=e.read_text();bubble=getblock(a,'figure','fig:sadb_clusters');explorer=getblock(s,'figure','fig:sadb_explorer')
s=s.replace(explorer,bubble);a=a.replace(bubble,explorer)
s=s.replace('exploring citation neighborhoods (Figure~\\ref{fig:sadb_explorer})','exploring citation neighborhoods (Figure~\\ref{fig:sadb_clusters})')
s=s.replace('Paper matching uses DOI where available and normalized titles otherwise.','Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\\ref{app:sadb_records} (Figure~\\ref{fig:sadb_explorer}).')
s=s.replace('Papers hold bibliographic information and short mechanism summaries grounded in the source text.','Papers hold bibliographic information and curation notes: short summaries of the reported mechanism, experimental conditions, and relevance to feedback design, grounded in the source text.')
a=a.replace("Figure~\\ref{fig:sadb_clusters} shows the application's paper-level neighborhood view.",r"The citation-neighborhood view is shown in Chapter~\ref{ch:database} (Figure~\ref{fig:sadb_clusters}); Figure~\ref{fig:sadb_explorer} provides the complementary paper-table view.")
p.write_text(s,encoding='utf-8');e.write_text(a,encoding='utf-8')
p=C/'chapters/10-introduction.tex';s=p.read_text().replace("The author's experience",'My experience')
s=paragraph(s,'Neural control introduces two further requirements',r'''Neural control exposes two further gaps. First, the project lacked a structured, searchable resource connecting measured sensory-afferent properties to the species, preparations, and pathway roles on which controller design depends. Those findings were distributed across primary papers and reviews \citep{van_de_crommert_review_1998, windhorst_muscle_2007, windhorst_spinal_2022}. Second, the project lacked tested neuromechanical models linking a spinal CPG to both a human-scale musculoskeletal model and its BPA-routed robot counterpart. Published CPG architectures and simulation tools supplied components \citep{rybak_modelling_2006, mccrea_organization_2008, markin_afferent_2010, nourse_sns_2023}, but their availability did not establish controller behavior after model conversion, actuator substitution, or the addition of contact feedback. The SADb and simulation chapters address these evidence and implementation gaps, respectively.''')
p.write_text(s,encoding='utf-8')
p=C/'chapters/15-background.tex';s=p.read_text();cat=getblock(s,'figure','fig:decerebrate_cat');s=s.replace(cat,cap(cat,r'Restored and colorized frame from archival footage of a decerebrate-cat treadmill experiment \citep{decerebrate_cat_footage}.','Decerebrate-cat treadmill experiment'))
b1=getblock(s,'figure','fig:rybak_cpg_circuit');b2=getblock(s,'figure','fig:shevtsova_spinal_circuit')
b1=cap(b1,r'Two-level CPG: rhythm generation (RG), pattern formation (PF), motor output, and afferent feedback. Reproduced from \citet[Fig.~1]{rybak_modelling_2006}.','Rybak two-level locomotor CPG concept')
b2=cap(b2,r'Extended CPG with pattern-formation and motor populations for the two-joint posterior biceps--semitendinosus (PBSt) and rectus femoris (RF). Reproduced from \citet[Fig.~5.8]{shevtsova_modeling_2016}.','Shevtsova extended CPG model with two-joint-muscle circuits')
b1=b1.replace(r'\begin{figure}[htbp]',r'\begin{figure}[p]').replace(r'width=0.55\textwidth',r'width=0.80\textwidth,height=0.30\textheight,keepaspectratio').replace(r'\end{figure}','')
b2=b2.replace(r'\begin{figure}[htbp]','').replace(r'width=0.6\textwidth',r'width=0.88\textwidth,height=0.31\textheight,keepaspectratio')
old1=getblock(s,'figure','fig:rybak_cpg_circuit');old2=getblock(s,'figure','fig:shevtsova_spinal_circuit');s=s.replace(old1,b1+'\n'+r'\par\medskip'+'\n'+b2).replace(old2,'')
# Keep the computational-neuroscience section general; CPG-specific structures remain in their own sections.
st=s.index('Computational neuroscience expresses');en=s.index(r'\section{Synthetic Nervous Systems}',st)
s=s[:st]+r'''Computational neuroscience expresses hypotheses about neural circuitry as executable models. In conductance-based models, membrane voltages evolve through capacitance, leak, synaptic, and applied currents, and connections represent identifiable pathways \citep{nourse_sns_2023}. Model predictions can be compared with recordings or with responses to stimulation and pathway removal. This provides a way to test a proposed mechanism rather than infer it from a wiring diagram alone. Synthetic nervous systems use this approach to control physical and simulated bodies.

'''+s[en:]
p.write_text(s,encoding='utf-8')
p=C/'chapters/70-overall.tex';s=p.read_text()
s=s.replace("and the monthly skimmer--curator workflow sweeps the recent literature for candidate records and routes every one of them to a human ruling.","The monthly literature-scraping extension is under development; it is separate from the completed curation and application described here.")
s=paragraph(s,"The database's growth now has two engines",r'''A monthly literature-scraping extension is being implemented to identify candidate additions to the database. Its retrieval and classification must be evaluated against human decisions before it can support routine updates. Further curation can also target missing afferent classifications, pathway links, and preparation metadata in the existing records. The AI-assisted triage study establishes the relevant division of responsibility: automated tools draft and retrieve, while I retain the inclusion and interpretation decisions \citep{bolen_using_2024}.''')
p.write_text(s,encoding='utf-8')
# Reuse the script-generated, committed physical-test plots byte-for-byte.
manifest={}
for source,name in [(R/'Code/Matlab/Mesh_Optimization/Results/Flx20mm_compare.png','flexor_20mm_physical.png'),(R/'Code/Matlab/Mesh_Optimization/Results/Ext20mm_compare_20261022.png','extensor_20mm_physical.png')]:
    shutil.copy2(source,C/'figs/30-torque'/name);manifest[name]={'source':str(source),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
(S/'physical_figure_sources.json').write_text(json.dumps(manifest,indent=2))
for name in ['actuator_technology_comparison','rybak_cpg_circuit','shevtsova_spinal_circuit']:
    d=fitz.open(C/'figs/15-background'/f'{name}.pdf');d[0].get_pixmap(matrix=fitz.Matrix(1,1)).save(S/f'{name}.png')
print('Updated knowledge-gap narratives, bubble placement, cat caption, paired CPG figures, and active-work status')
