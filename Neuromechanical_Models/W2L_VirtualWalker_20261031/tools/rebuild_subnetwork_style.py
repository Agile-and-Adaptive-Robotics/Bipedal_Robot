"""Reorganize unchanged W2L circuits into the historical uniarticular LH hierarchy."""
from pathlib import Path
from historical_build_guard import require_historical_regeneration
require_historical_regeneration()
from lxml import etree as E
import copy,uuid,json,hashlib,re
from collections import Counter,defaultdict
W=Path(__file__).resolve().parents[1];R=W.parents[1];B=W/'build/subnetwork_rebuild_20261031'
P=E.XMLParser(strip_cdata=False)
old=E.parse(str(B/'23b0e785.aproj'),P).getroot()
template=next(n for n in old.findall('.//NervousSystem//Node') if n.findtext('Text')=='LH')
offtemplate=next(n for n in old.findall('.//NervousSystem//Node') if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.OffPage')
historical_drawings={}
for node in old.findall('.//NervousSystem//Node'):
 if node.find('DiagramXml') is not None:
  for shape in E.fromstring(node.findtext('DiagramXml').encode()).findall('.//AddFlow/Node'):historical_drawings[shape.findtext('Tag')]=shape
def put(n,k,v):
 c=n.find(k)
 if c is None:c=E.SubElement(n,k)
 c.text=str(v);return c
def uid(folder,label):return str(uuid.uuid5(uuid.NAMESPACE_URL,'ben-w2l-unisubsystems-20261031/'+folder+'/'+label))
def canon(e):return E.tostring(e,method='c14n')
def group(name):
 m=re.match(r'^([LR])[ _]',name or '')
 if not m:return None
 s=m.group(1);base=(s+'H',);low=name.lower()
 if ' c1' in low or ' v3' in low:return ('Commissural coordination',)
 if ' rg ' in low:return base+(s+'H_RG',)
 if 'heel' in low or 'toe contact' in low:return base+(s+'H_Contact feedback',)
 if 'hip' in low:return base+(s+'H_HipZ',)
 ka=base+(s+'H_K&A Synergy',)
 if 'knee' in low:return ka+(s+'H_K&A Pattern Formation' if ' pf ' in low else s+'H_Knee Motoneuron',)
 if 'ank' in low or 'foot' in low or 'toe' in low:return ka+(s+'H_Anklez Motoneuron',)
 return None
def build(folder,name):
 d=R/'Neuromechanical_Models'/folder;src=B/folder/(name+'.aproj');tree=E.parse(str(src),P);root=tree.getroot()
 top=next(n for n in root.findall('.//NervousSystem//Node') if n.findtext('ClassName')=='AnimatGUI.DataObjects.Behavior.Nodes.Subsystem')
 nodes=list(top.find('Nodes'));links=list(top.find('Links'));byid={n.findtext('ID'):n for n in nodes};homes={i:group(n.findtext('Text')) for i,n in byid.items()};neighbours=defaultdict(list)
 for link in links:
  a,b=link.findtext('OriginID'),link.findtext('DestinationID');neighbours[a].append(b);neighbours[b].append(a)
 # Numeric adapter labels are resolved from their attached named circuit objects.
 for _ in range(4):
  for i,home in list(homes.items()):
   if home is None:
    votes=Counter(homes.get(j) for j in neighbours[i] if homes.get(j))
    if votes:homes[i]=votes.most_common(1)[0][0]
 unassigned=[byid[i].findtext('Text') for i,h in homes.items() if h is None]
 for i in homes:
  if homes[i] is None:homes[i]=('Auxiliary interfaces',)
 corrections=[]
 for i,n in byid.items():
  label=n.findtext('Text') or ''
  if label=='L Knee MN ext RE':
   near=Counter((byid[j].findtext('Text') or '')[:1] for j in neighbours[i] if j in byid)
   if near['R']>=3 and near['R']>near['L']:
    homes[i]=('RH','RH_K&A Synergy','RH_Knee Motoneuron')
    corrections.append({'id':i,'legacy_label':label,'placement':'RH knee MN; label retained in all files','neighbour_side_counts':dict(near)})
 # All physical bodies and neural module parameters must remain identical.
 bodies_before=[canon(n) for n in root.iter('RigidBody')]
 modules_before=[canon(n) for n in root.findall('.//NervousSystem/NeuralModules')]
 original_neural={i:canon(n) for i,n in byid.items()}
 source_draw=E.fromstring(top.findtext('DiagramXml').encode());source_flow=source_draw.find('.//AddFlow')
 shapes={n.findtext('Tag'):copy.deepcopy(n) for n in source_flow if n.tag in ('Node','Link')}
 fallback_node=next(copy.deepcopy(n) for n in source_flow if n.tag=='Node')
 fallback_link=next(copy.deepcopy(n) for n in source_flow if n.tag=='Link')
 for tag in ('Nodes','Links'):top.find(tag).clear()
 pages={():top}
 def page(path):
  if path in pages:return pages[path]
  parent=page(path[:-1]);n=copy.deepcopy(template)
  for tag in ('Nodes','Links','InLinks','OutLinks'):n.find(tag).clear()
  put(n,'ID',uid(folder,'page/'+('/'.join(path))));put(n,'Text',path[-1]);put(n,'ToolTip','Historical uniarticular LH organization; current functional parameters preserved.')
  parent.find('Nodes').append(n);pages[path]=n;return n
 for i,n in byid.items():page(homes[i]).find('Nodes').append(n)
 aliases={}
 def endpoint(i,path):
  if homes[i]==path:return i
  key=(i,path)
  if key in aliases:return aliases[key]
  n=copy.deepcopy(offtemplate);a=uid(folder,'offpage/'+i+'/'+('/'.join(path)));put(n,'ID',a);put(n,'Text',byid[i].findtext('Text'));put(n,'LinkedNodeID',i)
  put(n,'ToolTip','Same original node on '+ '/'.join(homes[i])+'; this connector adds no neuron or synapse.')
  for tag in ('InLinks','OutLinks'):n.find(tag).clear()
  page(path).find('Nodes').append(n);aliases[key]=a;return a
 preserved=[]
 for link in links:
  a,b=link.findtext('OriginID'),link.findtext('DestinationID');ha,hb=homes[a],homes[b];shared=[]
  for x,y in zip(ha,hb):
   if x!=y:break
   shared.append(x)
  path=tuple(shared);linkid=link.findtext('ID');preserved.append((linkid,a,b,link.findtext('SynapticTypeID')))
  put(link,'OriginID',endpoint(a,path));put(link,'DestinationID',endpoint(b,path));page(path).find('Links').append(link)
 # Rebuild actual drawing references and functional endpoint registrations.
 all_nodes={n.findtext('ID'):n for p in pages.values() for n in p.findall('Nodes/Node')}
 for n in all_nodes.values():
  for tag in ('InLinks','OutLinks'):
   z=n.find(tag)
   if z is None:z=E.SubElement(n,tag)
   z.clear()
 for path,p in pages.items():
  ns=list(p.find('Nodes'));ls=list(p.find('Links'));doc=copy.deepcopy(source_draw)
  if path:put(doc.find('Diagram'),'ID',uid(folder,'diagram/'+('/'.join(path))))
  put(doc.find('Diagram'),'PageName',p.findtext('Text'))
  flow=doc.find('.//AddFlow');flow.clear();flow.attrib.update({'Nodes':str(len(ns)),'Links':str(len(ls))});E.SubElement(flow,'Version').text='1.5.0.1'
  subsystem_index=0;alias_index=0
  for index,n in enumerate(sorted(ns,key=lambda n:(n.findtext('ClassName') or '',n.findtext('Text') or '',n.findtext('ID')))):
   i=n.findtext('ID');klass=n.findtext('ClassName') or '';label=n.findtext('Text') or '';alias='OffPage' in klass
   if 'Subsystem' in klass:x,y=90+(subsystem_index%3)*300,80+(subsystem_index//3)*150;subsystem_index+=1
   else:
    # Paired E/F populations share rows; connectors and physical interfaces stay in side columns.
    is_f=False if 'plantar' in label.lower() else bool(re.search(r'flx|flex|dorsi|\bF\b',label,re.I))
    if label.isdigit():
     words=[byid[j].findtext('Text') or '' for j in neighbours[i] if j in byid];is_f=sum(bool(re.search(r'flx|flex|dorsi',s,re.I)) for s in words)>sum(bool(re.search(r'ext|plantar',s,re.I)) for s in words)
    x=650 if is_f else 350
    if alias:x=40+(alias_index%4)*260;y=1000+(alias_index//4)*85;alias_index+=1
    elif 'NodeToPhysicalAdapter' in klass:x=840 if is_f else 180;y=650
    elif 'PhysicalToNodeAdapter' in klass:x=840 if is_f else 180;y=300+(index%3)*110
    elif 'Muscle' in klass or 'StretchReceptor' in klass:x=1000 if is_f else 20;y=650 if 'Muscle' in klass else 420
    elif ' PF ' in label:y=80+(110 if ' IN' in label else 0)
    elif ' RG ' in label:y=100+(130 if ' IN' in label else 0)
    elif ' MN ' in label:y=650+(110 if ' RE' in label else 0)
    elif 'Ia' in label:y=420+(110 if ' IN' in label else 0)
    elif 'Ib' in label:y=300
    else:x=140+(index%4)*200;y=80+(index//4)*110
   put(n,'Location','').attrib.update({'x':str(x),'y':str(y)});put(n,'Size','').attrib.update({'Width':'150' if alias else '125','Height':'50'})
   default_shape=historical_drawings[offtemplate.findtext('ID')] if alias else historical_drawings[template.findtext('ID')] if 'Subsystem' in klass else fallback_node
   draw=copy.deepcopy(shapes.get(i,default_shape));draw.attrib.update({'Left':str(x),'Top':str(y),'Width':'150' if alias else '125','Height':'50'});put(draw,'Tag',i);put(draw,'Text',label)
   if i not in shapes:
    for field in ('FillColor','DrawColor','TextColor'):
     if n.find(field) is not None:put(draw,field,n.findtext(field))
   flow.append(draw)
  # Keep overview-page references immediately below the subsystem boxes rather
  # than carrying the larger cell-level page's vertical spacing into every view.
  if any('OffPage' in (n.findtext('ClassName') or '') for n in ns):
   nonaliases=[n for n in ns if 'OffPage' not in (n.findtext('ClassName') or '')]
   base=max((float(n.find('Location').get('y'))+float(n.find('Size').get('Height')) for n in nonaliases),default=80)+130
   alias_nodes=sorted((n for n in ns if 'OffPage' in (n.findtext('ClassName') or '')),key=lambda n:n.findtext('Text') or '')
   for j,n in enumerate(alias_nodes):
    y=base+(j//4)*85;x=40+(j%4)*260;n.find('Location').attrib.update({'x':str(x),'y':str(y)})
    shape=next(z for z in flow if z.tag=='Node' and z.findtext('Tag')==n.findtext('ID'));shape.attrib.update({'Left':str(x),'Top':str(y)})
  idx={n.findtext('Tag'):k for k,n in enumerate([z for z in flow if z.tag in ('Node','Link')])}
  for link in ls:
   i=link.findtext('ID');a,b=link.findtext('OriginID'),link.findtext('DestinationID');E.SubElement(all_nodes[a].find('OutLinks'),'ID').text=i;E.SubElement(all_nodes[b].find('InLinks'),'ID').text=i
   draw=copy.deepcopy(shapes.get(i,fallback_link));draw.attrib.update({'Org':str(idx[a]),'Dst':str(idx[b])});put(draw,'Tag',i)
   for pt in draw.findall('Point'):draw.remove(pt)
   flow.append(draw)
  p.find('DiagramXml').text=E.CDATA(E.tostring(doc,encoding='unicode'))
 # Validate every page, effective original edge, and original parameter payload.
 alias_resolve={a:i for (i,path),a in aliases.items()};now={l.findtext('ID'):l for p in pages.values() for l in p.findall('Links/Link')}
 for i,a,b,t in preserved:
  l=now[i];assert alias_resolve.get(l.findtext('OriginID'),l.findtext('OriginID'))==a;assert alias_resolve.get(l.findtext('DestinationID'),l.findtext('DestinationID'))==b;assert l.findtext('SynapticTypeID')==t
 def neuron_params(n):
  n=copy.deepcopy(n)
  for tag in ('Location','Size','InLinks','OutLinks'):
   q=n.find(tag)
   if q is not None:n.remove(q)
  return canon(n)
 assert all(neuron_params(E.fromstring(original_neural[i]))==neuron_params(n) for i,n in byid.items())
 assert bodies_before==[canon(n) for n in root.iter('RigidBody')];assert modules_before==[canon(n) for n in root.findall('.//NervousSystem/NeuralModules')]
 drawing_count=0
 for path,p in pages.items():
  doc=E.fromstring(p.findtext('DiagramXml').encode());flow=doc.find('.//AddFlow');seq=[z for z in flow if z.tag in ('Node','Link')];ns=p.findall('Nodes/Node');ls=p.findall('Links/Link');assert len(ns)==int(flow.get('Nodes'));assert len(ls)==int(flow.get('Links'));tags={q.findtext('Tag'):q for q in seq};assert set(tags)=={q.findtext('ID') for q in ns+ls}
  for n in ns:assert tags[n.findtext('ID')].findtext('Text')==n.findtext('Text')
  for l in ls:
   q=tags[l.findtext('ID')];assert seq[int(q.get('Org'))].findtext('Tag')==l.findtext('OriginID');assert seq[int(q.get('Dst'))].findtext('Tag')==l.findtext('DestinationID')
  drawing_count+=1
 ids=[n.findtext('ID') for n in root.findall('.//NervousSystem//Node')];assert len(ids)==len(set(ids))
 output=d/(name+'.aproj');tree.write(str(output),encoding='utf-8',xml_declaration=True)
 E.parse(str(output),P)
 # Runnable neural and physical objects do not change when only diagram pages change.
 for suffix in ('.asim','_Standalone.asim'):
  original=(B/folder/(name+suffix)).read_bytes();assert (d/(name+suffix)).read_bytes()==original
 report=dict(machine='easteregg2',template_commit='23b0e785',template_scope='historical uniarticular LH subsystem organization, not copied physical model or neuron parameters',project=str(output),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),original_objects=len(nodes),original_links=len(links),added_neurons=0,added_functional_connections=0,offpage_connectors=len(aliases),pages={'/'.join(path) or 'root':{'nodes':len(p.findall('Nodes/Node')),'links':len(p.findall('Links/Link'))} for path,p in pages.items()},validated_drawing_pages=drawing_count,physical_payloads_unchanged=True,neuron_parameters_unchanged=True,effective_edges_unchanged=True,standalone_bytes_unchanged=True,root_y_telemetry='already present: root_y column in Virtual walker verification',legacy_label_findings=corrections,unassigned_before_auxiliary_group=unassigned)
 report['original_object_class_counts']=dict(Counter(n.findtext('ClassName') for n in nodes));report['actual_neuron_count']=sum('.Neurons.' in (n.findtext('ClassName') or '') for n in nodes);report['original_link_class_counts']=dict(Counter(n.findtext('ClassName') for n in links))
 (d/'SUBNETWORK_REBUILD.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
for args in [('W2L_VirtualWalker_20261031','W2L_VirtualWalker'),('W2L_VirtualWalker_Newer_20261031','W2L_VirtualWalker_Newer')]:build(*args)
