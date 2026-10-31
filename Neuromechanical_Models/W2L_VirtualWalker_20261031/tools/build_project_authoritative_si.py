"""Export the editable NEWER project's neural graph into a separate SI-corrected companion."""
from pathlib import Path
from historical_build_guard import require_historical_regeneration
require_historical_regeneration()
from lxml import etree as E
import copy,json,hashlib,uuid,shutil,argparse
from collections import Counter
parser=argparse.ArgumentParser();parser.add_argument('--variant-id',default='SIContact');parser.add_argument('--contact-na',type=float,default=20.);parser.add_argument('--offset-na',type=float,default=-.4);args=parser.parse_args()
assert args.variant_id.replace('_','').isalnum()
R=Path(__file__).resolve().parents[3];S=R/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031';NAME='W2L_VirtualWalker_Newer_'+args.variant_id;W=R/('Neuromechanical_Models/'+NAME+'_20261031');W.mkdir(exist_ok=True)
P=E.XMLParser(strip_cdata=False);contact_c=args.contact_na*1e-9;contact_d=args.offset_na*1e-9
pt=E.parse(str(S/'W2L_VirtualWalker_Newer.aproj'),P);p=pt.getroot();st=E.parse(str(S/'W2L_VirtualWalker_Newer_Standalone.asim'),P);s=st.getroot()
def put(n,k,v):
 z=n.find(k)
 if z is None:z=E.SubElement(n,k)
 z.text=str(v);return z
def val(n,k,default=None):
 z=n.find(k)
 if z is None:return default
 if z.get('Actual') is not None:return float(z.get('Actual'))
 try:return float(z.text)
 except:return z.text
def num(n,k,v,scale='None',display=None):
 z=put(n,k,'');z.attrib.update({'Value':str(v if display is None else display),'Scale':scale,'Actual':str(v)});return z
def uid(label):return str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-w2l-project-authoritative-si-20261031/'+label))
pn={n.findtext('ID'):n for n in p.findall('.//NervousSystem//Node')};aliases={i:n.findtext('LinkedNodeID') for i,n in pn.items() if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.OffPage'}
def resolve(i):
 while i in aliases:i=aliases[i]
 return i
def real(i):
 i=resolve(i);n=pn[i];return n.findtext('LinkedBodyPartID') or i
contact=[]
for n in pn.values():
 if n.findtext('DataTypeID')=='ContactCount' and n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.PhysicalToNodeAdapter':
  g=n.find('Gain');contact.append({'project_adapter_id':n.findtext('ID'),'target_id':n.findtext('DestinationID'),'source_C':val(g,'C'),'source_D':val(g,'D'),'candidate_C':contact_c,'candidate_D':contact_d})
  num(g,'C',contact_c,'nano',args.contact_na);num(g,'D',contact_d,'nano',args.offset_na)
assert len(contact)==4
put(p,'ProjectName',NAME);put(p,'SimulationFile',NAME+'.asim');put(p.find('Simulation'),'ProjectPath',str(W)+'\\');put(s,'ProjectPath',str(W)+'\\')
for f in S.glob('*.aform'):shutil.copy2(f,W/f.name)
for f in S.glob('*.bmp'):shutil.copy2(f,W/f.name)
shutil.copy2(S/'diagnostic_columns.json',W/'diagnostic_columns.json')
# The GUI's actual neuron properties are the editable authority, including all intrinsic values.
neurons={n.findtext('ID'):n for n in s.iter('Neuron')};changes=[]
fields={'RestingPotential':('RestingPot',1000),'RelativeSize':('Size',1),'TimeConstant':('TimeConst',1000),'InitialThreshold':('InitialThresh',1000),'RelativeAccomodation':('RelativeAccom',1),'AccomodationTimeConstant':('AccomTimeConst',1000),'AHP_Conductance':('AHPAmp',1e6),'AHP_TimeConstant':('AHPTimeConst',1000),'MaxCaConductance':('GMaxCa',1e6),'TonicStimulus':('TonicStimulus',1e9),'TonicNoise':('Noise',1000)}
gui_neurons={i:n for i,n in pn.items() if '.Neurons.' in (n.findtext('ClassName') or '')};assert set(neurons)==set(gui_neurons)
for i,n in gui_neurons.items():
 a=neurons[i];put(a,'Name',n.findtext('Text'));put(a,'Enabled',n.findtext('Enabled','True'))
 for g,(r,m) in fields.items():
  v=val(n,g)
  if v is not None:
   v=float(v)*m;before=float(a.findtext(r,'0'))
   if abs(v-before)>1e-12:changes.append({'neuron':i,'field':r,'runtime_before':before,'project_value':v})
   put(a,r,format(v,'.17g'))
 put(a,'BurstInitAtBottom',n.findtext('InitAtBottom','True'))
 for gate in ('CaActivation','CaDeactivation'):
  src=n.find(gate);dst=a.find(gate)
  if src is not None:
   if dst is None:dst=E.SubElement(a,gate)
   put(dst,'ID',src.findtext('ID'))
   for key,m in [('MidPoint',1000),('Slope',1),('TimeConstant',1000)]:put(dst,key,format(float(val(src,key))*m,'.17g'))
   put(dst,'ActivationType',src.findtext('ActivationType'))
# Export all synapse types actually referenced by GUI edges. Source native-only types are removed.
gui_links=[l for l in p.findall('.//NervousSystem//Links/Link') if l.findtext('ClassName')=='IntegrateFireGUI.DataObjects.Behavior.Synapse'];needed={l.findtext('SynapticTypeID') for l in gui_links}
gui_types={n.findtext('ID'):n for n in p.findall('.//SynapseTypes/Link')};native_types={n.findtext('ID'):n for n in s.iter('SynapseType')};containers={n.findtext('Type'):n.getparent() for n in native_types.values()}
type_changes=[]
for container in set(containers.values()):
 for n in list(container):
  if n.tag=='SynapseType':container.remove(n)
for i in sorted(needed):
 n=gui_types[i];kind=n.findtext('ClassName').split('.')[-1];out=copy.deepcopy(native_types[i]) if i in native_types else E.Element('SynapseType');put(out,'ID',i);put(out,'Name',n.findtext('Name'));put(out,'Type',kind)
 if kind=='NonSpikingChemical':
  for g,r,m in [('EquilibriumPotential','Equil',1000),('MaxSynapticConductance','SynAmp',1e6),('PreSynapticThreshold','ThreshV',1000),('PreSynapticSaturationLevel','SaturateV',1000)]:
   v=float(val(n,g))*m;old=out.findtext(r)
   if old is None or abs(float(old)-v)>1e-12:type_changes.append({'type':i,'field':r,'runtime_before':old,'project_value':v})
   put(out,r,format(v,'.17g'))
 elif i not in native_types:raise ValueError('Need explicit exporter for new type '+kind)
 containers[kind].append(out)
# Rewrite connection objects from the GUI graph; OffPage references resolve to original neurons.
nc=s.find('.//Connexions');assert nc is not None
before_edges=Counter((n.findtext('SourceID'),n.findtext('TargetID')) for n in nc.findall('Connexion'));before_ids={n.findtext('ID'):n for n in nc}
for n in list(nc):nc.remove(n)
for l in gui_links:
 n=E.SubElement(nc,'Connexion');i=l.findtext('ID');a,b=resolve(l.findtext('OriginID')),resolve(l.findtext('DestinationID'));kind=gui_types[l.findtext('SynapticTypeID')].findtext('ClassName').split('.')[-1]
 for k,v in [('ID',i),('SourceID',a),('TargetID',b),('Type',{'NonSpikingChemical':1,'Electrical':2,'SpikingChemical':0}[kind]),('SynapseTypeID',l.findtext('SynapticTypeID')),('Delay',float(val(l,'ConductionDelay',0))*1000),('G',float(val(l,'SynapticConductance',0))*1e6)]:put(n,k,v)
after_edges=Counter((n.findtext('SourceID'),n.findtext('TargetID')) for n in nc)
# GUI adapters include the actual physical source/target bindings and complete gain triplets.
ac=s.find('.//Adapters');assert ac is not None
old_adapters=list(ac);adapter_template=next(n for n in old_adapters if n.findtext('Type')=='PhysicalToNode')
for n in list(ac):ac.remove(n)
gain_skip={'AssemblyFile','ClassName','LimitsReadOnly','LimitOutputsReadOnly','UseParentIncomingDataType'};gain_rename={'LowerLimitScale':'LowerLimit','UpperLimitScale':'UpperLimit','LowerOutputScale':'LowerOutput','UpperOutputScale':'UpperOutput'}
adapter_changes=[]
for i,n in pn.items():
 klass=n.findtext('ClassName') or ''
 if not klass.endswith(('PhysicalToNodeAdapter','NodeToPhysicalAdapter')):continue
 physical='PhysicalToNodeAdapter' in klass;kind='PhysicalToNode' if physical else 'NodeToPhysical';a,b=real(n.findtext('OriginID')),real(n.findtext('DestinationID'));source_type=n.findtext('DataTypeID');target_type=n.findtext('TargetDataTypeID')
 old=next((q for q in old_adapters if q.findtext('Type')==kind and q.findtext('SourceID')==a and q.findtext('TargetID')==b),None)
 out=copy.deepcopy(old if old is not None else adapter_template)
 for key,v in [('ID',i),('Name',n.findtext('Text')),('Type',kind),('SourceModule','AnimatLab' if physical else 'IntegrateFireSim'),('TargetModule','IntegrateFireSim' if physical else 'AnimatLab'),('SourceID',a),('TargetID',b),('SourceDataType',source_type),('TargetDataType',target_type)]:put(out,key,v)
 gain=out.find('Gain');gain.clear()
 for q in n.find('Gain'):
  if q.tag in gain_skip:continue
  # Gain parameters are scalar SI values, unlike the IntegrateFire module's mV/ms/nA convention.
  put(gain,gain_rename.get(q.tag,q.tag),q.get('Actual') if q.get('Actual') is not None else q.text or '')
 ac.append(out);adapter_changes.append({'gui_id':i,'runtime_old_id':None if old is None else old.findtext('ID'),'kind':kind,'source':a,'target':b,'gain_type':gain.findtext('Type')})
# Additional headless-only neural telemetry does not change the model's dynamics.
charts=s.find('DataCharts');base=next(c for c in charts if c.find('DataColumns') is not None);chart=copy.deepcopy(base);put(chart,'ID',uid('neuralchart'));put(chart,'Name','SI contact neural check');put(chart,'OutputFilename','SI contact neural check.txt')
columns=chart.find('DataColumns');columns.clear();ct=next(c for c in s.findall('.//DataColumn') if c.findtext('DataType')=='MembraneVoltage')
for i,n in gui_neurons.items():
 c=copy.deepcopy(ct);put(c,'ID',uid('neuralcolumn/'+i));put(c,'ColumnName',n.findtext('Text')+' ['+i[:8]+']');put(c,'TargetID',i);columns.append(c)
charts.append(chart)
pt.write(str(W/(NAME+'.aproj')),encoding='utf-8',xml_declaration=True);st.write(str(W/(NAME+'_Standalone.asim')),encoding='utf-8',xml_declaration=True);shutil.copy2(W/(NAME+'_Standalone.asim'),W/(NAME+'.asim'))
report={'machine':'easteregg2','authority':'requested editable .aproj neural graph and properties','variant_name':NAME,'project':str(W/(NAME+'.aproj')),'standalone':str(W/(NAME+'_Standalone.asim')),'neurons':len(gui_neurons),'neural_connections':len(gui_links),'contact_corrections':contact,'neuron_parameter_changes_from_previous_runtime':changes,'synapse_type_changes_from_previous_runtime':type_changes,'runtime_only_endpoints_removed':[{'source':a,'target':b,'source_name':pn[a].findtext('Text'),'target_name':pn[b].findtext('Text')} for a,b in (before_edges-after_edges).elements()],'project_only_endpoints_added':[{'source':a,'target':b,'source_name':pn[a].findtext('Text'),'target_name':pn[b].findtext('Text')} for a,b in (after_edges-before_edges).elements()],'adapter_bindings':adapter_changes,'source_hashes':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (S/'W2L_VirtualWalker_Newer.aproj',S/'W2L_VirtualWalker_Newer_Standalone.asim')},'physics_payloads':'source standalone retained; only ProjectPath, neural export, adapters, and diagnostic chart changed','GUI_open':'not verified: CUA access denied'}
report['contact_gain_rationale']={'C_nA_per_contact':args.contact_na,'D_nA':args.offset_na,'interpretation':'1 nA/contact, zero offset preserves the GUI scalar magnitudes under an explicit nA-unit assumption; 20 nA/contact with -0.4 nA offset is a distinct engineering candidate derived from the old runtime scalar values. Neither is asserted to be Ben\'s intended or experimentally identified gain.','is_pure_GUI_magnitude_preserving_reinterpretation':args.contact_na==1 and args.offset_na==0}
(W/'BUILD_SYNC.json').write_text(json.dumps(report,indent=2));(W/'.gitignore').write_text('*.aproj\n*.asim\n*.aform\n*.bmp\n*.npz\n*.txt\n*.log\n!README.txt\n')
print(json.dumps({k:v for k,v in report.items() if k not in ('adapter_bindings','synapse_type_changes_from_previous_runtime')},indent=2))
