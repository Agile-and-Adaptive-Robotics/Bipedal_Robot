from lxml import etree as E
from pathlib import Path
p=E.parse(r'D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Deng_Biarticular_Rat\Biarticular_net.aproj')
alln=p.findall('.//NervousSystem//Node');pages={n.findtext('Text').strip():n for n in alln if n.find('DiagramXml') is not None}
for key in ['LH_RG','LH_Hip','Biar']:
 pg=pages[key]; nodes={n.findtext('ID'):n for n in pg.findall('Nodes/Node')};print('\nPAGE',key)
 for i,n in nodes.items():print('NODE',i,n.findtext('Text').strip(),n.findtext('ClassName'))
 for e in pg.findall('Links/Link'):
  print('EDGE',nodes[e.findtext('OriginID')].findtext('Text').strip(),'->',nodes[e.findtext('DestinationID')].findtext('Text').strip(),e.findtext('ArrowDestination/Style'),e.findtext('SynapseTypeID'))
