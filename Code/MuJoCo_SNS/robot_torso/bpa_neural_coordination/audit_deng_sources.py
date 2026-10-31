"""Read-only inventory of the downloaded ModelDB archive and its XML models."""
from pathlib import Path
import sys
sys.stdout.reconfigure(encoding='utf-8',errors='replace')
import zipfile,hashlib,json,xml.etree.ElementTree as ET
HERE=Path(__file__).resolve().parent
REPO=HERE.parents[3]
archive=REPO/'Code/Matlab/SNS_Simscape/dev/deng2019_modeldb_264844.zip'
OUT=HERE/'runs/deng_source_audit_20261031';OUT.mkdir(parents=True,exist_ok=True)
summary=dict(archive=str(archive),sha256=hashlib.sha256(archive.read_bytes()).hexdigest(),models=[])
with zipfile.ZipFile(archive) as z:
 summary['entries']=z.namelist()
 for name in z.namelist():
  if name.lower().endswith(('.asim','.aproj')):
   b=z.read(name)
   try:root=ET.fromstring(b)
   except ET.ParseError:continue
   interesting=[]
   for e in root.iter():
    t=e.findtext('Name') or e.findtext('Text')
    if t and any(k in t.lower() for k in ('ankle','knee','gas','tib','rect','bicep','femor','biart','syner')):
     interesting.append(dict(tag=e.tag,name=t,type=e.findtext('Type'),id=e.findtext('ID')))
   muscles=[e for e in root.iter('RigidBody') if e.findtext('Type')=='LinearHillMuscle']
   summary['models'].append(dict(entry=name,bytes=len(b),root=root.tag,muscle_count=len(muscles),muscle_names=[e.findtext('Name') for e in muscles],interesting=interesting))
   (OUT/(hashlib.sha256(name.encode()).hexdigest()[:8]+'_'+Path(name).name+'.txt')).write_text(b.decode('utf-8-sig',errors='replace'),encoding='utf-8')
(OUT/'ARCHIVE_INVENTORY.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(dict(archive_sha256=summary['sha256'],entries=len(summary['entries']),models=[{'entry':x['entry'],'bytes':x['bytes'],'interesting':x['interesting'][:30]} for x in summary['models']]),indent=2))
with zipfile.ZipFile(archive) as z:
 name='K_A wirh RG  march 21/K_A synergy.aproj';root=ET.fromstring(z.read(name))
source=Path('D:/GitHub/_zcode_backups/lit_20261026/kaiyu2023.docx')
with zipfile.ZipFile(source) as z:
 root=ET.fromstring(z.read('word/document.xml'));ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 paragraphs=[''.join(t.text or '' for t in p.findall('.//w:t',ns)) for p in root.findall('.//w:p',ns)]
 (OUT/'kaiyu2023_visible_paragraphs.txt').write_text('\n'.join(f'P{i+1}: {p}' for i,p in enumerate(paragraphs)),encoding='utf-8')
 for i,p in enumerate(paragraphs):
  if any(k in p.lower() for k in ('github','modeldb','gastrocnemius','plantar','dorsi','appendix','available')):print('PAPER2023',i+1,p[:700])
with zipfile.ZipFile(archive) as z:
 root=ET.fromstring(z.read('K_A wirh RG  march 21/K_A synergy.aproj'))
nodes={e.findtext('ID'):e for e in root.iter('Node')}
def label(e):return e.findtext('Name') or e.findtext('Text') or ''
def resolve(key):
 visited=set()
 while key in nodes and key not in visited:
  visited.add(key);e=nodes[key]
  if e.findtext('LinkedNodeID'):key=e.findtext('LinkedNodeID')
  else:return label(e)
 return key
edges=[]
for e in root.iter('Link'):
 src=e.findtext('OriginID');dst=e.findtext('DestinationID')
 if src and dst:edges.append(dict(source=resolve(src),destination=resolve(dst),type=e.findtext('ClassName'),conductance=None if e.find('SynapticConductance') is None else e.find('SynapticConductance').attrib))
(OUT/'2019_RESOLVED_LINKS.json').write_text(json.dumps(edges,indent=2))
for e in edges:
 if 'Ankle' in (e['source'] or '') or 'Ankle' in (e['destination'] or ''):print('EDGE2019',json.dumps(e))
for e in root.iter('Link'):
 sid=e.findtext('OriginID');did=e.findtext('DestinationID')
 if sid in nodes and resolve(sid)=='LH_AnkleZ MN E' and did in nodes:
  n=nodes[did]
  print('MN_E_TARGET',json.dumps({c.tag:c.text for c in n if c.tag not in ['DiagramXml','Font']}))
for e in root.iter('Node'):
 if label(e)=='LH_Ankle Ext':print('PHYSICAL_NODE_EXT',json.dumps({c.tag:c.text for c in e if c.tag not in ['DiagramXml','Font']}))
for e in root.iter('RigidBody'):
 if e.findtext('Type')=='LinearHillMuscle':print('MUSCLE2019',label(e))
for e in root.iter('Node'):
 if 'Ankle' in label(e) and any('Adapter' in (c.text or '') for c in e if c.tag=='ClassName'):print('ADAPTER2019',ET.tostring(e,encoding='unicode')[-2500:])
later=ET.parse(OUT/'Biarticular_net_Standalone.asim').getroot()
for tag in ('Neuron','NeuralModule','Adapter'):
 examples=list(later.iter(tag))
 print('LATER_TAG',tag,len(examples),[dict(name=e.findtext('Name'),children=[c.tag for c in e]) for e in examples[:1]])
for e in later.iter('RigidBody'):
 if e.findtext('Type')=='LinearHillMuscle':print('MUSCLE_LATER',label(e))
lnodes={e.findtext('ID'):e.findtext('Name') for e in later.iter('Neuron')}
lsyn={e.findtext('ID'):dict(name=e.findtext('Name'),g=e.findtext('SynAmp'),equil=e.findtext('Equil')) for e in later.iter('SynapseType')}
later_edges=[]
for e in later.iter('Connexion'):
 src=lnodes.get(e.findtext('SourceID'),'');dst=lnodes.get(e.findtext('TargetID'),'')
 if 'PF' in src and any(k in dst for k in ('SO','TA','RF','BFP','GA')):
  row=dict(source=src,target=dst,synapse=lsyn.get(e.findtext('SynapseTypeID')));later_edges.append(row);print('LATER_EDGE',json.dumps(row))
(OUT/'PUBLIC_BIARTICULAR_SELECTED_EDGES.json').write_text(json.dumps(later_edges,indent=2))
with zipfile.ZipFile(source) as z:
 doc=ET.fromstring(z.read('word/document.xml'))
 for i,tbl in enumerate(doc.findall('.//w:tbl',ns)):
  rows=[' | '.join(''.join(t.text or '' for t in c.findall('.//w:t',ns)) for c in row.findall('w:tc',ns)) for row in tbl.findall('w:tr',ns)]
  (OUT/f'2023_table_{i+1}.txt').write_text('\n'.join(rows),encoding='utf-8')
public=[]
for name,blob in [('Biarticular_net_Standalone.asim','afdc829ced0ead90ebf7f0552da6edc6ee504219'),('Biarticular_net.aproj','e0b3123513594407f6f09f1e87773f59d75a5dfb')]:
 p=OUT/name
 if p.exists():
  b=p.read_bytes()
  # apply_patch normalizes line endings. Restore original CRLF only if the
  # resulting bytes exactly match GitHub's immutable blob hash.
  restored=b.replace(b'\r\n',b'\n').replace(b'\n',b'\r\n')
  for candidate in (restored,restored.rstrip(b'\r\n'),restored.rstrip(b'\r\n')+b'\r\n'):
   h=hashlib.sha1(b'blob '+str(len(candidate)).encode()+b'\0'+candidate).hexdigest()
   if h==blob:
    b=candidate;p.write_bytes(b);break
  git_hash=hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()
  public.append(dict(file=str(p),source_url='https://github.com/vipamol/Rat-SNS-design-and-tuning/blob/main/Biarticular_net/'+name,
      source_git_blob_sha=blob,cache_git_blob_sha=git_hash,cache_exact=git_hash==blob,sha256=hashlib.sha256(b).hexdigest(),bytes=len(b)))
(OUT/'PUBLIC_SOURCE_PROVENANCE.json').write_text(json.dumps(public,indent=2))
print('COUNTS',[(m['entry'],m['muscle_count']) for m in summary['models']]);print('PUBLIC_PROVENANCE',json.dumps(public))
