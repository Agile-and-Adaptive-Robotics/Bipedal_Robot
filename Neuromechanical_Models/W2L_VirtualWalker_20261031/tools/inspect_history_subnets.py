from pathlib import Path
import subprocess,json,hashlib,shutil
from lxml import etree as E
R=Path(__file__).resolve().parents[3];W=Path(__file__).resolve().parents[1];O=W/'build/subnetwork_rebuild_20261031';O.mkdir(parents=True,exist_ok=True)
for folder,name in [('W2L_VirtualWalker_20261031','W2L_VirtualWalker'),('W2L_VirtualWalker_Newer_20261031','W2L_VirtualWalker_Newer')]:
 dest=O/folder;dest.mkdir(exist_ok=True)
 for ext in ('.aproj','.asim','_Standalone.asim'):
  p=R/'Neuromechanical_Models'/folder/(name+ext)
  if not (dest/p.name).exists():shutil.copy2(p,dest/p.name)
stats=[]
for rev in ('23b0e785','2e24466e','b29ba681','bcbeab88','6a0fd1eb'):
 listing=subprocess.check_output(['git','ls-tree','-r','--name-only',rev],cwd=R,text=True).splitlines();p=next((s for s in listing if s.endswith('Biped_2xCPG_wSubs.aproj')),None)
 if not p:continue
 raw=subprocess.check_output(['git','show',rev+':'+p],cwd=R);(O/(rev+'.aproj')).write_bytes(raw);root=E.fromstring(raw)
 rows=[]
 for n in root.findall('.//NervousSystem//Node'):
  rows.append(dict(id=n.findtext('ID'),name=n.findtext('Text'),klass=n.findtext('ClassName'),children=len(n.findall('Nodes/Node')),links=len(n.findall('Links/Link'))))
 stats.append(dict(rev=rev,path=p,sha256=hashlib.sha256(raw).hexdigest(),nodes=rows))
(O/'history_inspection.json').write_text(json.dumps(stats,indent=2))
for item in stats:
 print(item['rev'],len(item['nodes']))
 for n in item['nodes']:
  if n['children']:print(' SUB',n)
root=E.parse(str(W/'W2L_VirtualWalker.aproj')).getroot()
print('CURRENT',[(n.findtext('Text'),n.findtext('ClassName')) for n in root.findall('.//NervousSystem//Node')][:8])
