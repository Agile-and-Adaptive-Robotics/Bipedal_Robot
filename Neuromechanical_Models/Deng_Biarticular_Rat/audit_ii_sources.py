"""Static source-wiring inspection only; no model edits or runs."""
from pathlib import Path
import json,xml.etree.ElementTree as E
HERE=Path(__file__).resolve().parent;REPO=HERE.parents[1]
def inspect(p):
 r=E.parse(p).getroot();nodes={n.findtext('ID'):n.findtext('Name') for n in r.iter('Neuron')};bodies={n.findtext('ID'):n.findtext('Name') for n in r.iter('RigidBody')};types={n.findtext('ID'):n for n in r.iter('SynapseType')};out=[]
 for a in r.iter('Adapter'):
  if a.findtext('SourceDataType')!='II':continue
  target=a.findtext('TargetID');g=a.find('Gain');nexts=[]
  for c in r.iter('Connexion'):
   if c.findtext('SourceID')!=target:continue
   s=types.get(c.findtext('SynapseTypeID'));nexts.append(dict(target=nodes.get(c.findtext('TargetID')),synapse=s.findtext('Name') if s is not None else None,equilibrium=s.findtext('Equil') if s is not None else None,gmax=s.findtext('SynAmp') if s is not None else None))
  out.append(dict(body=bodies.get(a.findtext('SourceID')),target=nodes.get(target),enabled=a.findtext('Enabled'),gain={k:g.findtext(k) for k in ('Type','A','B','C','D')} if g is not None else {},outgoing=nexts))
 return dict(source=str(p),II_adapters=out)
result=[inspect(HERE/'Biarticular_net_Standalone.asim')]
for p in (REPO/'Neuromechanical_Models/W2L_VirtualWalker_20261031').glob('*.asim'):
 if 'Standalone' in p.name:result.append(inspect(p));break
(HERE/'runs/II_SOURCE_WIRING_AUDIT.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
