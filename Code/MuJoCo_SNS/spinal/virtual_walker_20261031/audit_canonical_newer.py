"""Read-only effective GUI source audit; standalone is compared, not assumed."""
import json,hashlib,xml.etree.ElementTree as E
from pathlib import Path
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[3];folder=REPO/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031';project=folder/'W2L_VirtualWalker_Newer.aproj';native=folder/'W2L_VirtualWalker_Newer_Standalone.asim'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def actual(e):
 if e is None:return None
 return float(e.get('Actual',e.get('Value',e.text or '0')))
g=E.parse(project).getroot();nodes={n.findtext('ID'):n for n in g.iter('Node')};neurons={i:n for i,n in nodes.items() if '.Neurons.' in n.findtext('ClassName','')};types={n.findtext('ID'):n for n in g.iter('Link') if n.findtext('ClassName','').endswith('SynapseTypes.NonSpikingChemical')};bodies={b.findtext('ID'):b for b in g.iter('RigidBody')};relevant={i for i,n in neurons.items() if n.findtext('Text','').endswith(' II') or 'contact' in n.findtext('Text','').lower() or ' RG ' in n.findtext('Text','')};edges=[]
def resolve(i):
 seen=set()
 while i in nodes and nodes[i].findtext('LinkedNodeID'):
  assert i not in seen;seen.add(i);i=nodes[i].findtext('LinkedNodeID')
 return i
for link in g.iter('Link'):
 if link.findtext('ClassName','')!='IntegrateFireGUI.DataObjects.Behavior.Synapse':continue
 pre=resolve(link.findtext('OriginID'));post=resolve(link.findtext('DestinationID'));typ=types[link.findtext('SynapticTypeID')];equil=actual(typ.find('EquilibriumPotential'))*1000;edges.append(dict(id=link.findtext('ID'),source_id=pre,source=neurons[pre].findtext('Text'),target_id=post,target=neurons[post].findtext('Text'),type_id=typ.findtext('ID'),type_name=typ.findtext('Name'),effective_gmax_uS=actual(typ.find('MaxSynapticConductance'))*1e6,local_GUI_conductance_uS=actual(link.find('SynapticConductance'))*1e6,reversal_mV=equil,threshold_mV=actual(typ.find('PreSynapticThreshold'))*1000,saturation_mV=actual(typ.find('PreSynapticSaturationLevel'))*1000,sign_at_minus60mV='depolarizing' if equil>-60 else 'hyperpolarizing' if equil< -60 else 'zero driving force'))
adapters=[]
for i,n in nodes.items():
 if not n.findtext('ClassName','').endswith('PhysicalToNodeAdapter'):continue
 target=resolve(n.findtext('DestinationID'));source=resolve(n.findtext('OriginID'));physical=nodes.get(source);bodyid=physical.findtext('LinkedBodyPartID') if physical is not None else None;body=bodies.get(bodyid);gain=n.find('Gain');drawn=[]
 for link in g.iter('Link'):
  if link.findtext('OriginID')==i:drawn.append(dict(link_id=link.findtext('ID'),target_id=resolve(link.findtext('DestinationID')),target=nodes[resolve(link.findtext('DestinationID'))].findtext('Text')))
 adapters.append(dict(id=i,GUI_name=n.findtext('Text'),source_gui_node=source,source_gui_name=physical.findtext('Text') if physical is not None else None,physical_id=bodyid,physical_name=body.findtext('Name') if body is not None else None,source_data_type=n.findtext('DataTypeID'),target_id=target,target=neurons[target].findtext('Text'),polynomial_ABCD_A=[actual(gain.find(k)) for k in 'ABCD'],coefficient_prefixes={k:dict(gain.find(k).attrib) for k in 'ABCD'},drawn_adapter_outputs=drawn,visual_binding_agrees=all(x['target_id']==target for x in drawn)))
selected=[]
for i in sorted(relevant):
 n=neurons[i];selected.append(dict(id=i,name=n.findtext('Text'),rest_mV=actual(n.find('RestingPotential'))*1000,tonic_nA=actual(n.find('TonicStimulus'))*1e9,incoming=[e for e in edges if e['target_id']==i],outgoing=[e for e in edges if e['source_id']==i],physical_inputs=[a for a in adapters if a['target_id']==i]))
r=E.parse(native).getroot();native_neurons={n.findtext('ID'):n.findtext('Name') for n in r.iter('Neuron')};native_edges={x.findtext('ID'):(x.findtext('SourceID'),x.findtext('TargetID')) for x in r.iter('Connexion')};gui_edges={x['id']:(x['source_id'],x['target_id']) for x in edges};result=dict(machine='easteregg2',authority='canonical newer .aproj, never the original80 or generated sibling',project=str(project),project_sha256=sha(project),standalone=str(native),standalone_sha256=sha(native),GUI_neurons=len(neurons),GUI_edges=len(edges),native_neurons=len(native_neurons),native_edges=len(native_edges),GUI_native_endpoint_equal=gui_edges==native_edges,neural_channels=selected,physical_adapters=adapters,visual_adapter_binding_mismatches=[a for a in adapters if not a['visual_binding_agrees']],nonspiking_strength_rule='Actual nonspiking conductance is the referenced synapse type MaxSynapticConductance/SynAmp, scaled by presynaptic threshold/saturation. Link-local SynapticConductance/native G does not replace it; verified IntegrateFireModule.cpp lines1204–1221.')
out=HERE/'runs/canonical_newer_audit';out.mkdir(exist_ok=True);(out/'AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:result[k] for k in ('GUI_neurons','GUI_edges','native_neurons','native_edges','GUI_native_endpoint_equal','project_sha256','standalone_sha256')},indent=2))
