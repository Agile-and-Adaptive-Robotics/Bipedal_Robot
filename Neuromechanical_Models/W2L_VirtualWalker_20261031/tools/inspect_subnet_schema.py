from pathlib import Path
from lxml import etree as E
from collections import Counter
W=Path(__file__).resolve().parents[1];O=W/'build/subnetwork_rebuild_20261031';old=E.parse(str(O/'23b0e785.aproj')).getroot()
for label in ('LH','LH_RG','LH_HipZ','LH_K&A Synergy','LH_Knee Motoneuron'):
 n=next(n for n in old.findall('.//NervousSystem//Node') if n.findtext('Text')==label)
 print('\nSUBSYSTEM',label,'fields',[e.tag for e in n]);print('DIRECT NODES',[(q.findtext('Text'),q.findtext('ClassName')) for q in n.findall('Nodes/Node')]);print('DRAWING',n.findtext('DiagramXml')[:230] if n.find('DiagramXml') is not None else None)
off=next(n for n in old.findall('.//NervousSystem//Node') if 'OffPage' in (n.findtext('ClassName') or ''));print('\nOFFPAGE',E.tostring(off).decode())
link=old.find('.//NervousSystem//Links/Link');print('LINK',E.tostring(link).decode())
n=next(n for n in old.findall('.//NervousSystem//Node') if n.findtext('Text')=='LH_RG');c=E.fromstring(n.findtext('DiagramXml').encode());print('DRAWINGS',E.tostring(c).decode()[:6000])
for fname in ['W2L_VirtualWalker.aproj','../W2L_VirtualWalker_Newer_20261031/W2L_VirtualWalker_Newer.aproj']:
 p=W/fname;r=E.parse(str(p)).getroot();n=next(n for n in r.findall('.//NervousSystem//Node') if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.Subsystem');print('CURRENT',fname,len(n.findall('Nodes/Node')),len(n.findall('Links/Link')));print('NAMES',[(q.findtext('Text'),q.findtext('ClassName')) for q in n.findall('Nodes/Node')]);print('LINKFIELDS',[e.tag for e in n.find('Links/Link')])
