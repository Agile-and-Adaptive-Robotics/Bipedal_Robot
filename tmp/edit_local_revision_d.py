from pathlib import Path
import shutil,json,re,difflib,hashlib
R=Path(r'D:\Github\Bipedal_Robot')
E=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_review_20261023'
C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C'
D=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261023/Revision_D'
if D.exists(): raise RuntimeError('Refusing to overwrite an existing revision')
D.mkdir(parents=True); (D/'review').mkdir()
shutil.copytree(C/'chapters',D/'chapters')
shutil.copytree(C/'figs',D/'figs')
for p in C.iterdir():
    if p.is_file() and p.suffix.lower() in ['.tex','.bib','.bst','.cls']: shutil.copy2(p,D/p.name)
for p in (E/'matrix_source').glob('*.tex'): shutil.copy2(p,D/'chapters'/p.name)
shutil.copytree(D/'chapters',D/'review/before')
changes=[]
def replace(file,old,new,ids,reason):
    p=D/'chapters'/file; s=p.read_text(encoding='utf-8')
    assert s.count(old)==1,(file,old[:80],s.count(old))
    p.write_text(s.replace(old,new),encoding='utf-8')
    changes.append(dict(file=file,ids=ids,before=old,after=new,reason=reason))
def para(file,start):
    s=(D/'chapters'/file).read_text(encoding='utf-8')
    return next(p for p in s.split('\n\n') if p.startswith(start))

f='10-introduction.tex'
replace(f,'Attempts to implement neural control models exposes two further gaps.','Attempts to implement neural control models exposed two further gaps.',['INTRO-grammar'],'Correct subject/verb agreement and tense.')
replace(f,'and was difficult to sift through for new researchers to the field.','and were difficult for researchers entering the field to locate and compare.',['INTRO-grammar'],'Make the dispersed findings the explicit plural subject.')
replace(f,'they had not been combined into a complete into a complete human neuromechanical model that could be used for rigurous testing of different hypotheses.','this project still needed to combine and test them in the human and BPA-routed robot models.',['INTRO-scope'],'Remove duplicate text and misspelling; scope the missing combination to this project instead of asserting an unverified field-wide first.')

f='50-database.tex'
replace(f,'Archived records remain available alongside the retained collection.','Records removed from the working Airtable base remain available in the application as archived records. This separates the collection used for ongoing curation from the broader set of records considered during its assembly.',['M-064'],'Define archived by its actual app-only status; source sadb_app.html archive control/header.')
replace(f,'The pivot and citation views operate on papers; the table view searches across record types.','The pivot and citation views operate on papers; the table view searches across record types. For a selected paper, references are the works it cites, whereas citing papers are works that cite it. The citation view distinguishes these directions.',['M-063','M-069'],'Explain citation direction; verified build_citation_graph.py:8 and build_app.py:320-322,725-744.')
old=para(f,'The working base is contained')
new=r'''The working base is contained in Airtable, an online database program. Its linked tables organize Papers, Feedback Pathways, Models, and Review Papers. A paper record contains bibliographic information and curation notes: short summaries of its reported mechanism, experimental conditions, and relevance to feedback design. A feedback-pathway record describes a sensory influence on motor control, such as Ib feedback influencing the stance-to-swing transition, and links that description to supporting papers. A model record connects a computational implementation to its source literature; a review record summarizes literature relevant to the collection. Thus, a pathway identifies a mechanism, the linked papers supply its experimental support, and model records identify implementations to examine. These links organize evidence; they do not establish that every relationship has been entered or verified. Paper matching uses DOI where available and normalized titles otherwise. The paper-table view is shown in Appendix~\ref{app:sadb_records} (Figure~\ref{fig:sadb_explorer}).'''
replace(f,old,new,['M-065','DB-22'],'Explain record types and relationships in reading order; preserve existing example and explicitly avoid claiming complete links.')
replace(f,'Classification used controlled vocabularies for preparations and afferent classes.','Classification used controlled vocabularies to group papers by their experimental preparation and afferent class, rather than relying on variations in free-text descriptions.',['DB-methods'],'Explain the purpose of controlled vocabulary without inventing category definitions.')
replace(f,'selected animal preparations as columns.','selected experimental preparations as columns.',['M-087'],'Use the existing preparation category without implying that every category is an animal species.')
old=para(f,'The present corpus is uneven.')
new=old+' Further work should check and complete the links between feedback pathways and their supporting papers, and between reviews, papers, and pathways. Each link needs to identify the evidence for the relationship, so that a reader can distinguish a directly measured sensory effect from a model assumption or a review summary. The priority is to verify these relationships and fill gaps in the existing collection before treating it as a complete basis for controller design.'
replace(f,old,new,['M-091'],'Answer the explicit future-work request without claiming all requested links are absent.')

f='60-sns.tex'
old=para(f,'The missing step was a tested link')
new='The AnimatLab result established that the combined model could produce stepping-like motion under suspension, but it did not establish weight-bearing walking. It also left unanswered how the controller would behave with the human musculotendon mechanics or BPA routing used elsewhere in this project. Investigating those questions required a way to connect the circuit to each mechanical model, inspect its signals, and test the connections before interpreting whole-model behavior.'
replace(f,old,new,['M-105','M-106'],'Replace abstract missing-step language with observed outcome, limitation, and specific requirement; the existing 92-actuator value remains later in this chapter.')
old=para(f,'I developed several software tools')
new=r'''I developed tools for connecting neural controllers to the mechanical models used in this project. The work addressed five integration tasks:
\begin{enumerate}
    \item Convert the Gait2392 human model from OpenSim to MuJoCo and connect it to an SNS-Toolbox CPG controller.
    \item Couple a BPA-routed robot counterpart to a CPG controller so that artificial-actuator mechanics can be examined.
    \item Combine the Li walker mechanics, Deng neural organization, and foot-contact feedback, and test the resulting controller and mechanical interfaces.
    \item Construct neural circuits graphically in Simulink and connect them to mechanical simulations.
    \item Bring SolidWorks-derived mechanisms and MATLAB BPA models into the same development workflow as the neural controller.
\end{enumerate}
These tasks have different completion levels. The converted models, actuator-force calculation, graphical neural blocks, and selected interface tests are described below. Importing a mechanism or reproducing a rhythm does not establish that a complete coupled model walks. Table~\ref{tab:sns_capability} distinguishes the demonstrated interfaces from the combinations that still require integration and testing.'''
replace(f,old,new,['M-107','M-108','M-109','M-110','M-111','M-112','M-113','M-114','M-115','M-116','M-117','M-118','M-119','M-120','M-121','M-122'],'Keep the advisor five-part structure but distinguish tasks from completed systems, remove tested-it ellipsis, and define graphical construction.')
replace(f,'Mechanical imports & OpenSim-derived and SolidWorks-derived mechanisms in Simscape Multibody ...&','Mechanical imports & OpenSim-derived and SolidWorks-derived mechanisms in Simscape Multibody &',['SNS-table'],'Remove literal ellipsis from table cell; preserve all capabilities and limits.')
old=para(f,'MyoConverter supplies the initial')
new=r'''MyoConverter supplies the initial OpenSim-to-MuJoCo translation. I added procedures for repairing custom muscle paths and checking joint conventions, then coupled the converted model to an SNS-Toolbox controller. The converted model carries 92 muscle-tendon actuators. Comparing attachment paths with the OpenSim source checks the transferred geometry; activating individual actuators checks the direction of the resulting joint torque. These checks concern conversion fidelity, not the ability of the controller to produce walking. The neural and mechanical models exchange signals at a \qty{2}{\milli\second} step.'''
replace(f,old,new,['M-123','M-124','M-125','M-126'],'Define checks rather than implying unspecified audit discoveries; actual defects and repairs remain a ZCode evidence request.')
old=para(f,'The BPA interface evaluates')
new=r'''For the BPA-routed robot model, an actuator-force calculation connects the neural and mechanical simulations. It takes the simulated muscle length and pressure as inputs and evaluates the quasi-static force model from Chapter~\ref{ch:muscles}; the resulting force acts along the simulated actuator path. This calculation does not include valve or pressure dynamics. Length, velocity, force, and foot-contact signals return to the controller through sensory channels. The correction terms of Chapter~\ref{ch:torque} do not by themselves establish that transferred routing and compliance are correct in a dynamic model; those properties require separate checks.'''
replace(f,old,new,['M-127'],'Define BPA interface inputs, force output, feedback, and limits at first mention without changing model functionality.')

f='15-background.tex'
replace(f,r'SNS Toolbox \citep{nourse_sns_2023}',r'SNS-Toolbox supplies the neuron and synapse models used to simulate the controller \citep{nourse_sns_2023}.',['BG-fragment'],'Complete a sentence fragment with the existing scoped attribution.')
old=para(f,'I used these tools to develop')
new=r'''I used these tools to develop a robotic OpenSim leg model, convert it to MuJoCo, and couple it to a synthetic nervous system. Chapter~\ref{ch:sns} describes the assembled interfaces and their evaluation. The coupled simulation supports repeated parameter, stimulation, and pathway-removal tests under controlled conditions. Whether those tests support a biological hypothesis depends on the model's validation and the experimental comparison, as well as simulation speed.'''
replace(f,old,new,['BG-throughput'],'Remove unsupported overnight throughput/statistical-power assertions; explicitly logged, not silently harmonized. No numerical results or figures removed.')
old='A robot that implements a validated neuromechanical hypothesis is precisely such a replacement method: every deletion, lesion, and stimulation experiment that runs on the synthetic nervous system is an experiment that does not require a living subject, and the robot population can be duplicated exactly, so that ``subjects\'\' are statistically identical across conditions, an experimental repeatability no breeding colony can offer. The failure modes of these experiments damage hardware, not tissue.'
new='Robotic and computational models allow proposed neural interventions to be tested under controlled conditions before corresponding animal experiments are considered. The same model and initial conditions can be reused across interventions to isolate the effect of a changed pathway. Repeating a model does not create independent biological subjects, and its relevance to animal experiments depends on validation against biological observations.'
replace(f,old,new,['BG-ethics-repeatability'],'Clarify repeatability versus biological replication; retain Three Rs definition/citations and the existing figure.')

for c in changes:
    s=(D/'chapters'/c['file']).read_text(encoding='utf-8'); c['new_line']=s[:s.index(c['after'])].count('\n')+1
(D/'review/changes.json').write_text(json.dumps(changes,indent=2,ensure_ascii=False),encoding='utf-8')
diff=''; conservation={}
for file in ['10-introduction.tex','15-background.tex','50-database.tex','60-sns.tex']:
    b=(D/'review/before'/file).read_text(encoding='utf-8'); a=(D/'chapters'/file).read_text(encoding='utf-8')
    diff+=''.join(difflib.unified_diff(b.splitlines(True),a.splitlines(True),fromfile='before/'+file,tofile='after/'+file))
    patterns={'labels':r'\\label\{[^}]*\}','citations':r'\\cite\w*(?:\[[^\]]*\])*\{[^}]*\}','figures':r'\\begin\{figure\}.*?\\end\{figure\}'}
    conservation[file]={k:re.findall(p,b,re.S)==re.findall(p,a,re.S) for k,p in patterns.items()}
    assert all(conservation[file].values()),(file,conservation[file])
(D/'review/changes.diff').write_text(diff,encoding='utf-8')
(D/'review/conservation.json').write_text(json.dumps(conservation,indent=2),encoding='utf-8')
print(json.dumps({'working_copy':str(D),'edit_count':len(changes),'conservation':conservation},indent=2))
