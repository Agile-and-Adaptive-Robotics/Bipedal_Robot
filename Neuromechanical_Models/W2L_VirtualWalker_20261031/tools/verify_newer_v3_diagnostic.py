from pathlib import Path
from lxml import etree as E
from collections import Counter
import json,numpy as np,hashlib
R=Path(__file__).resolve().parents[3];W=R/'Neuromechanical_Models/W2L_VirtualWalker_NewerV3Guided_20261031';name='W2L_VirtualWalker_NewerV3Guided';p=E.parse(str(W/(name+'.aproj'))).getroot();s=E.parse(str(W/(name+'_Standalone.asim'))).getroot();meta=json.loads((W/'V3_DIAGNOSTIC_BUILD.json').read_text());S=Path(meta['source']);old=E.parse(str(S/'W2L_VirtualWalker_Newer_GUI1nA_Standalone.asim')).getroot();pn={n.findtext('ID'):n for n in p.findall('.//NervousSystem//Node')};aliases={i:n.findtext('LinkedNodeID') for i,n in pn.items() if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.OffPage'}
def resolve(i):
 while i in aliases:i=aliases[i]
 return i
def val(n,k):
 z=n.find(k);return float(z.get('Actual')) if z.get('Actual') is not None else float(z.text)
pl=[l for l in p.findall('.//NervousSystem//Links/Link') if l.findtext('ClassName')=='IntegrateFireGUI.DataObjects.Behavior.Synapse'];sl=list(s.iter('Connexion'));assert len(pl)==len(sl)==204
assert Counter((l.findtext('ID'),resolve(l.findtext('OriginID')),resolve(l.findtext('DestinationID')),l.findtext('SynapticTypeID')) for l in pl)==Counter((l.findtext('ID'),l.findtext('SourceID'),l.findtext('TargetID'),l.findtext('SynapseTypeID')) for l in sl)
gt={q.findtext('ID'):q for q in p.findall('.//SynapseTypes/Link')};st={q.findtext('ID'):q for q in s.iter('SynapseType')};edges={q.findtext('ID'):q for q in sl};oldedges={q.findtext('ID'):q for q in old.iter('Connexion')};oldtypes={q.findtext('ID'):q for q in old.iter('SynapseType')};zero=set(meta['changes'][0]['edges_disabled']);newzero=meta['changes'][0]['new_zero_type']
for i,t in st.items():
 for g,n,m in [('EquilibriumPotential','Equil',1000),('MaxSynapticConductance','SynAmp',1e6),('PreSynapticThreshold','ThreshV',1000),('PreSynapticSaturationLevel','SaturateV',1000)]:assert abs(val(gt[i],g)*m-float(t.findtext(n)))<1e-10
for i,l in edges.items():
 if i in zero:assert l.findtext('SynapseTypeID')==newzero and float(st[newzero].findtext('SynAmp'))==0
 else:
  t=l.findtext('SynapseTypeID');assert t==oldedges[i].findtext('SynapseTypeID');assert E.tostring(st[t],method='c14n')==E.tostring(oldtypes[t],method='c14n')
sn={n.findtext('ID'):n for n in s.iter('Neuron')};gn={i:n for i,n in pn.items() if '.Neurons.' in (n.findtext('ClassName') or '')};assert len(sn)==len(gn)==96
for i,n in gn.items():
 assert abs(val(n,'TimeConstant')*1000-float(sn[i].findtext('TimeConst')))<1e-10
 for gate in ('CaActivation','CaDeactivation'):assert abs(val(n,gate+'/TimeConstant')*1000-float(sn[i].findtext(gate+'/TimeConstant')))<1e-10
for label in ('L RG ext','R RG ext'):
 n=next(n for n in sn.values() if n.findtext('Name')==label);assert float(n.findtext('TonicStimulus'))==0;assert val(gn[n.findtext('ID')],'TonicStimulus')==0
contact=[q for q in s.findall('.//Adapters/Adapter') if q.findtext('SourceDataType')=='ContactCount'];assert len(contact)==4
for q in contact:assert float(q.findtext('Gain/C'))==1e-9 and float(q.findtext('Gain/D'))==0
robot=next(n for n in s.iter('RigidBody') if n.findtext('Name')=='Root');assert robot.find('Joint') is None and robot.findtext('Freeze')=='False'
arr=np.load(W/'Guide diagnostics.npz');a=arr['data'];labels=arr['labels'].tolist();g=lambda k:a[:,labels.index(k)];raw=np.load(W/'Virtual walker verification.npz');b=raw['data'];names=raw['labels'].tolist();knees=np.column_stack([10+np.degrees(b[:,names.index('knee_L_angle')]),55+np.degrees(b[:,names.index('knee_R_angle')])])
neur=np.loadtxt(W/'SI contact neural check.txt',skiprows=1);neur=neur[(neur[:,1]<=10)&np.any(neur[:,2:]!=0,axis=1)];headers=(W/'SI contact neural check.txt').open().readline().strip().split('\t');periods={}
for label in ('L RG flx','R RG flx'):
 j=next(i for i,h in enumerate(headers) if h.startswith(label+' ['));v=neur[:,j];t=neur[:,1];ix=np.flatnonzero((v[1:]>-.06)&(v[:-1]<=-.06))+1;onsets=t[ix];onsets=onsets[onsets>=1];periods[label]={'onsets_s':onsets.tolist(),'median_period_s':float(np.median(np.diff(onsets))) if len(onsets)>2 else None}
for n in p.findall('.//NervousSystem//Node'):
 if n.find('DiagramXml') is not None:E.fromstring(n.findtext('DiagramXml').encode())
report={'source_unchanged':all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in meta['source_hashes'].items()),'GUI_runtime_neural_graph_and_used_type_parameters_match':True,'only_two_effective_synapse_outputs_zeroed':sorted(zero),'other_edge_types_and_parameters_preserved':True,'c1_preserved':True,'RG_ext_tonic_nA':[0,0],'neurons':96,'connections_including_two_zero':204,'finite_physical':bool(np.isfinite(a).all()),'finite_all96_neural':bool(np.isfinite(neur).all()),'neural_range_mV':[float(neur[:,2:].min()*1000),float(neur[:,2:].max()*1000)],'root_y_min_final_m':[float(g('root_pos_Y').min()),float(g('root_pos_Y')[-1])],'root_x_change_m':float(g('root_pos_X')[-1]-g('root_pos_X')[0]),'knee_min_flexion_deg_LR':knees.min(axis=0).tolist(),'knee_hyperextension_samples_LR':np.sum(knees<0,axis=0).tolist(),'RG_periods':periods,'scope':'Native10s diagnostic transfer, not finalV3 architecture or crossengine walkingpass; no native qualified footclearance/forcegate evaluated.'}
(W/'V3_NATIVE_VERIFICATION.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
