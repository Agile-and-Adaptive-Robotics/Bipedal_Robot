from pathlib import Path
from lxml import etree as E
from collections import Counter
import numpy as np,json,hashlib,argparse
parser=argparse.ArgumentParser();parser.add_argument('--variant-id',default='SIContact');args=parser.parse_args();assert args.variant_id.replace('_','').isalnum()
R=Path(__file__).resolve().parents[3];name='W2L_VirtualWalker_Newer_'+args.variant_id;W=R/('Neuromechanical_Models/'+name+'_20261031');p=E.parse(str(W/(name+'.aproj'))).getroot();s=E.parse(str(W/(name+'_Standalone.asim'))).getroot();build=json.loads((W/'BUILD_SYNC.json').read_text());coeff=build['contact_corrections'][0];expected_c=coeff.get('candidate_C',coeff.get('corrected_C'));expected_d=coeff.get('candidate_D',coeff.get('corrected_D'))
pn={n.findtext('ID'):n for n in p.findall('.//NervousSystem//Node')};aliases={i:n.findtext('LinkedNodeID') for i,n in pn.items() if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.OffPage'}
def resolve(i):
 while i in aliases:i=aliases[i]
 return i
def number(n,k):
 z=n.find(k);return float(z.get('Actual')) if z.get('Actual') is not None else float(z.text)
gl=[l for l in p.findall('.//NervousSystem//Links/Link') if l.findtext('ClassName')=='IntegrateFireGUI.DataObjects.Behavior.Synapse'];sl=list(s.iter('Connexion'))
ge=Counter((l.findtext('ID'),resolve(l.findtext('OriginID')),resolve(l.findtext('DestinationID')),l.findtext('SynapticTypeID')) for l in gl);se=Counter((l.findtext('ID'),l.findtext('SourceID'),l.findtext('TargetID'),l.findtext('SynapseTypeID')) for l in sl);assert ge==se
gt={n.findtext('ID'):n for n in p.findall('.//SynapseTypes/Link')};st={n.findtext('ID'):n for n in s.iter('SynapseType')}
for i,t in st.items():
 if t.findtext('Type')=='NonSpikingChemical':
  for a,b,m in [('EquilibriumPotential','Equil',1000),('MaxSynapticConductance','SynAmp',1e6),('PreSynapticThreshold','ThreshV',1000),('PreSynapticSaturationLevel','SaturateV',1000)]:assert abs(number(gt[i],a)*m-float(t.findtext(b)))<1e-10
sa={n.findtext('ID'):n for n in s.iter('Adapter')};contact=[]
neuron_fields={'RestingPotential':('RestingPot',1000),'RelativeSize':('Size',1),'TimeConstant':('TimeConst',1000),'InitialThreshold':('InitialThresh',1000),'RelativeAccomodation':('RelativeAccom',1),'AccomodationTimeConstant':('AccomTimeConst',1000),'AHP_Conductance':('AHPAmp',1e6),'AHP_TimeConstant':('AHPTimeConst',1000),'MaxCaConductance':('GMaxCa',1e6),'TonicStimulus':('TonicStimulus',1e9),'TonicNoise':('Noise',1000)}
sn={n.findtext('ID'):n for n in s.iter('Neuron')}
for i,n in pn.items():
 if '.Neurons.' in (n.findtext('ClassName') or ''):
  for a,(b,m) in neuron_fields.items():
   if n.find(a) is not None:assert abs(number(n,a)*m-float(sn[i].findtext(b)))<1e-10
  for gate in ('CaActivation','CaDeactivation'):
   for a,m in [('MidPoint',1000),('Slope',1),('TimeConstant',1000)]:assert abs(number(n,gate+'/'+a)*m-float(sn[i].findtext(gate+'/'+a)))<1e-10
 if (n.findtext('ClassName') or '').endswith(('PhysicalToNodeAdapter','NodeToPhysicalAdapter')):
  for a in ('A','B','C','D'):
   if n.find('Gain/'+a) is not None:assert abs(number(n,'Gain/'+a)-float(sa[i].findtext('Gain/'+a)))<1e-14
for i,n in pn.items():
 if n.findtext('DataTypeID')=='ContactCount':
  assert float(sa[i].findtext('Gain/C'))==number(n,'Gain/C')==expected_c;assert float(sa[i].findtext('Gain/D'))==number(n,'Gain/D')==expected_d;contact.append(i)
ids=[n.text for n in s.iter('ID')];duplicates={i:k for i,k in Counter(ids).items() if k>1};assert not duplicates,duplicates
for n in p.findall('.//NervousSystem//Node'):
 if n.find('DiagramXml') is not None:E.fromstring(n.findtext('DiagramXml').encode())
a=np.load(W/'Virtual walker verification.npz')['data'];labels=json.loads((W/'diagnostic_columns.json').read_text());y=a[:,labels.index('root_y')];fall=np.flatnonzero(y<.8);j=int(fall[0]) if len(fall) else None
if j is not None:fall_info={k:float(a[j,labels.index(k)]) for k in ['Time','root_x','root_y','root_z','Left hip wing_x','Left hip wing_y','Left hip wing_z','Right hip wing_x','Right hip wing_y','Right hip wing_z']}
else:fall_info=None
evidence={'effective_edge_ID_endpoint_type_match':True,'used_synapse_type_parameters_match':True,'mapped_neuron_intrinsic_parameters_match':True,'exported_polynomial_adapter_coefficients_match':True,'contact_coefficients_match':True,'runtime_ID_duplicates':duplicates,'neural_chart_finite':bool(np.isfinite(np.load(W/'SI contact neural check.npz')['data']).all()),'gross_height_loss_detection':{'criterion':'root world Y < 0.8 m (descriptive fall flag, not a general balance gate)','first_sample':fall_info},'rail_world_x_interval_m':[-5.454,-1.454],'scope':'Native finite run is a failed supported-motion case; no gait or balance pass claimed. GUI opening not verified.','source_files_unchanged':{path:hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha for path,sha in json.loads((W/'BUILD_SYNC.json').read_text())['source_hashes'].items()}}
assert all(evidence['source_files_unchanged'].values())
(W/'EXPORT_VERIFICATION.json').write_text(json.dumps(evidence,indent=2));print(json.dumps(evidence,indent=2))
