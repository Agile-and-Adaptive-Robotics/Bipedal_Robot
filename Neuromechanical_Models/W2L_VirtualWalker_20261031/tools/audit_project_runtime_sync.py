from pathlib import Path
from lxml import etree as E
import json
from collections import Counter
R=Path(__file__).resolve().parents[3];W=R/'Neuromechanical_Models/W2L_VirtualWalker_Newer_20261031'
p=E.parse(str(W/'W2L_VirtualWalker_Newer.aproj')).getroot();s=E.parse(str(W/'W2L_VirtualWalker_Newer_Standalone.asim')).getroot()
for n in p.findall('.//NervousSystem//Node'):
 if n.findtext('DataTypeID')=='ContactCount':
  print('GUI_ADAPTER',n.findtext('ID'),n.findtext('Text'),'target',n.findtext('DestinationID'),'C',E.tostring(n.find('Gain/C')).decode(),'D',E.tostring(n.find('Gain/D')).decode())
pl={l.findtext('ID'):l for l in p.findall('.//NervousSystem//Links/Link') if l.findtext('ClassName')=='IntegrateFireGUI.DataObjects.Behavior.Synapse'}
sl={l.findtext('ID'):l for l in s.iter('Connexion')}
print('COUNTS',len(pl),len(sl));print('RUNTIME_ONLY',list(sl.keys()-pl.keys()));print('PROJECT_ONLY',list(pl.keys()-sl.keys()))
for i in sl.keys()-pl.keys():print('RUNTIME_EXTRA',E.tostring(sl[i]).decode())
aliases={n.findtext('ID'):n.findtext('LinkedNodeID') for n in p.findall('.//NervousSystem//Node') if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.OffPage'}
def resolve(i):return aliases.get(i,i)
pe=Counter((resolve(l.findtext('OriginID')),resolve(l.findtext('DestinationID')),l.findtext('SynapticTypeID')) for l in pl.values());se=Counter((l.findtext('SourceID'),l.findtext('TargetID'),l.findtext('SynapseTypeID')) for l in sl.values())
print('RUNTIME_ONLY_EFFECTIVE',list((se-pe).elements()));print('PROJECT_ONLY_EFFECTIVE',list((pe-se).elements()))
for i in pl.keys()-sl.keys():print('PROJECT_EXTRA',E.tostring(pl[i]).decode())
pn={n.findtext('ID'):n.findtext('Text') for n in p.findall('.//NervousSystem//Node')};sen=Counter((a,b) for a,b,c in se.elements());pen=Counter((a,b) for a,b,c in pe.elements())
print('RUNTIME_ONLY_ENDPOINTS',[(pn.get(a,a),pn.get(b,b),a,b) for a,b in (sen-pen).elements()]);print('PROJECT_ONLY_ENDPOINTS',[(pn.get(a,a),pn.get(b,b),a,b) for a,b in (pen-sen).elements()])
print('RUNTIME_TYPE',E.tostring(next(n for n in s.iter('SynapseType') if n.findtext('ID')=='cafe0043-0000-4000-8000-900000000043')).decode())
print('GUI_TYPE',E.tostring(next(n for n in p.findall('.//SynapseTypes/Link') if n.findtext('ID')=='cafe0044-0000-4000-8000-900000000044')).decode())
print('RUNTIME_NEURON',E.tostring(next(n for n in s.iter('Neuron') if n.findtext('ID')=='cafe0201-0000-4000-8000-900000000201')).decode()[:6500])
print('GUI_NEURON',E.tostring(next(n for n in p.findall('.//NervousSystem//Node') if n.findtext('ID')=='cafe0201-0000-4000-8000-900000000201')).decode())
