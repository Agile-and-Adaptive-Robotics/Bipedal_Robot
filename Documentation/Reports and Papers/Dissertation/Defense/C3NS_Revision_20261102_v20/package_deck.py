from pathlib import Path
import zipfile,json,hashlib,copy
from lxml import etree as E
W=Path(__file__).resolve().parent
helper=(W.parent/'C3NS_Revision_20261031_v18/transplant_v18.py').read_text(encoding='utf-8').split('src = read_package(SOURCE)')[0].replace("f'v18_s", "f'v20_s")
ctx={'__file__':str(W/'package_deck.py')}
exec(compile(helper,'package_helpers','exec'),ctx)
read=ctx['read_package'];parts=ctx['presentation_slide_parts'];rels=ctx['rels_part'];resolve=ctx['resolve_target'];N=ctx['NS']
B=W.parent/'C3NS_Revision_20261102_v19/output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx'
base=read(B);src=read(W/'build/artifact_draft.pptx');bp=parts(base);sp=parts(src)
changed=json.loads((W/'build/changed_slides.json').read_text())
order=json.loads((W/'build/slide_order.json').read_text())
assert len(order)==128 and set(order)==set(range(1,130))-{14}
pres=E.fromstring(base['ppt/presentation.xml']);presrels=E.fromstring(base['ppt/_rels/presentation.xml.rels'])
ids=pres.find('p:sldIdLst',N);id_nodes=list(ids);nextid=max(int(x.get('id')) for x in ids)+1
ct=E.fromstring(base['[Content_Types].xml']);CT='http://schemas.openxmlformats.org/package/2006/content-types'
for n in range(124,130):
 newpart=f'ppt/slides/slide{900+n}.xml';newnote=f'ppt/notesSlides/notesSlide{900+n}.xml'
 assert newpart not in base and newnote not in base
 template=bp[36];base[newpart]=base[template]
 tr=E.fromstring(base[rels(template)])
 oldnote=next(r for r in tr if ctx['relation_type'](r)=='notesSlide')
 oldnotepart=resolve(template,oldnote.get('Target'))
 base[newnote]=base[oldnotepart]
 nr=E.fromstring(base[rels(oldnotepart)])
 for r in nr:
  if ctx['relation_type'](r)=='slide':r.set('Target',ctx['relative_target'](newnote,newpart))
 base[rels(newnote)]=E.tostring(nr,xml_declaration=True,encoding='UTF-8',standalone=True)
 oldnote.set('Target',ctx['relative_target'](newpart,newnote))
 base[rels(newpart)]=E.tostring(tr,xml_declaration=True,encoding='UTF-8',standalone=True)
 rid=f'rIdV20_{n}'
 E.SubElement(presrels,'{'+N['pr']+'}Relationship',Id=rid,Type=N['r']+'/slide',Target=newpart.removeprefix('ppt/'))
 node=E.Element('{'+N['p']+'}sldId',id=str(nextid));node.set('{'+N['r']+'}id',rid);id_nodes.append(node);nextid+=1;bp.append(newpart)
 for p,t in [(newpart,'slide'),(newnote,'notesSlide')]:
  E.SubElement(ct,'{'+CT+'}Override',PartName='/'+p,ContentType=f'application/vnd.openxmlformats-officedocument.presentationml.{t}+xml')
base['[Content_Types].xml']=E.tostring(ct,xml_declaration=True,encoding='UTF-8',standalone=True)
for n in changed:
 if n!=2:
  ctx['transplant_slide'](n,src,base,sp,bp)
 else:
  # Preserve native animated GIF shape geometry, animation, and media.
  sr=E.fromstring(src[sp[1]]);br=E.fromstring(base[bp[1]])
  def title(root):
   for s in root.findall('./p:cSld/p:spTree/p:sp',N):
    off=s.find('./p:spPr/a:xfrm/a:off',N);ext=s.find('./p:spPr/a:xfrm/a:ext',N)
    if off is not None and ext is not None and int(off.get('y'))<760000 and int(ext.get('cx'))>6000000:return s
  target=title(br);source=title(sr)
  target.find('p:txBody',N).getparent().replace(target.find('p:txBody',N),copy.deepcopy(source.find('p:txBody',N)))
  base[bp[1]]=E.tostring(br,xml_declaration=True,encoding='UTF-8',standalone=True)
  a=next(r for r in E.fromstring(src[rels(sp[1])]) if ctx['relation_type'](r)=='notesSlide')
  b=next(r for r in E.fromstring(base[rels(bp[1])]) if ctx['relation_type'](r)=='notesSlide')
  base[resolve(bp[1],b.get('Target'))]=src[resolve(sp[1],a.get('Target'))]
for child in list(ids):ids.remove(child)
for old in order:ids.append(id_nodes[old-1])
base['ppt/presentation.xml']=E.tostring(pres,xml_declaration=True,encoding='UTF-8',standalone=True)
base['ppt/_rels/presentation.xml.rels']=E.tostring(presrels,xml_declaration=True,encoding='UTF-8',standalone=True)
mapping={old:i for i,old in enumerate(order,1)}
for old in order:
 root=E.fromstring(base[bp[old-1]]);edited=False
 for s in root.findall('./p:cSld/p:spTree/p:sp',N):
  off=s.find('./p:spPr/a:xfrm/a:off',N);texts=s.findall('.//a:t',N)
  if mapping[old]!=old and off is not None and int(off.get('y'))>4700000 and ''.join(x.text or '' for x in texts)==str(old):
   texts[0].text=str(mapping[old]);edited=True
 if edited:base[bp[old-1]]=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
# Duplicate media share unchanged bytes; obsolete slide-6 animal image is removed
# when no active relationship refers to it.
seen={};dupes={}
for p,data in base.items():
 if p.startswith('ppt/media/'):
  key=(Path(p).suffix.lower(),hashlib.sha256(data).hexdigest())
  if key in seen:dupes[p]=seen[key]
  else:seen[key]=p
for rp,data in list(base.items()):
 if not rp.endswith('.rels'):continue
 owner='' if rp=='_rels/.rels' else rp.replace('/_rels/','/').removesuffix('.rels');root=E.fromstring(data);edited=False
 for r in root:
  if r.get('TargetMode')=='External':continue
  target=resolve(owner,r.get('Target'))
  if target in dupes:r.set('Target',ctx['relative_target'](owner,dupes[target]));edited=True
 if edited:base[rp]=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
targets=set()
for rp,data in base.items():
 if not rp.endswith('.rels'):continue
 owner='' if rp=='_rels/.rels' else rp.replace('/_rels/','/').removesuffix('.rels')
 for r in E.fromstring(data):
  if r.get('TargetMode')!='External':targets.add(resolve(owner,r.get('Target')))
orphan=[p for p in base if p.startswith('ppt/media/') and p not in targets]
for p in orphan:del base[p]
# Ensure the original native chart, movie, and master parts are byte-identical.
original=read(B)
for p,data in original.items():
 if p.startswith(('ppt/charts/','ppt/embeddings/','ppt/slideMasters/','ppt/slideLayouts/')) or p.lower().endswith('.mp4'):
  assert base.get(p)==data,p
for n,p in enumerate(bp[:123],1):
 if n not in changed and mapping.get(n)==n:assert base[p]==original[p],n
out=W/'build/candidate.pptx'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for p,data in base.items():z.writestr(p,data)
(W/'build/package_report.json').write_text(json.dumps({'changed_source_slides':changed,'source_to_final':mapping,'main_slides':75,'question_slides':53,'total_slides':128,'orphan_media_removed':len(orphan),'bytes':out.stat().st_size},indent=2))
print('Packaged',out.stat().st_size,'bytes;',len(changed),'changed slides')
