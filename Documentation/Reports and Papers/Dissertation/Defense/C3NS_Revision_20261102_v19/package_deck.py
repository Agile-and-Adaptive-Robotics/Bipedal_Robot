from pathlib import Path
import json, zipfile, re, copy, posixpath, hashlib
from lxml import etree as E
W=Path(__file__).resolve().parent
old=W.parent/'C3NS_Revision_20261031_v18'
helper=(old/'transplant_v18.py').read_text(encoding='utf-8').split('src = read_package(SOURCE)')[0]
ctx={'__file__':str(W/'package_deck.py')}
exec(compile(helper,'transplant_helpers','exec'),ctx)
read=ctx['read_package']; parts=ctx['presentation_slide_parts'];rels=ctx['rels_part'];resolve=ctx['resolve_target'];kind=ctx['relation_type'];NS=ctx['NS']
base=read(old/'output/Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx')
src=read(W/'build/artifact_draft.pptx');bp=parts(base);sp=parts(src)
protected={1,2,16,26,31,44,52,74,97,98,112,123}
changed=set(range(1,71))|{76,77,78,79,82,83,87,89,92,95,103,104,108,109,111,118}
for n in sorted(changed):
 if n not in protected:
  ctx['transplant_slide'](n,src,base,sp,bp)
 else:
  a=E.fromstring(src[sp[n-1]]); b=E.fromstring(base[bp[n-1]])
  def istitle(s):
   off=s.find('./p:spPr/a:xfrm/a:off',NS); ext=s.find('./p:spPr/a:xfrm/a:ext',NS)
   return off is not None and ext is not None and int(off.get('y'))<760000 and int(ext.get('cx'))>6000000
  titles=[s for s in a.findall('./p:cSld/p:spTree/p:sp',NS) if istitle(s)]
  if titles:
   tree=b.find('./p:cSld/p:spTree',NS)
   for t in list(tree):
    if t.tag.endswith('}sp') and istitle(t):tree.remove(t)
   nextid=max([int(t.get('id')) for t in b.findall('.//p:cNvPr',NS)]+[1])+1
   for t in titles:
    node=copy.deepcopy(t)
    for nv in node.findall('.//p:cNvPr',NS):nv.set('id',str(nextid));nextid+=1
    tree.append(node)
  base[bp[n-1]]=E.tostring(b,xml_declaration=True,encoding='UTF-8',standalone=True)
  if n==1:
   tree=b.find('./p:cSld/p:spTree',NS)
   for node in a.findall('./p:cSld/p:spTree/p:sp',NS):
    nv=node.find('./p:nvSpPr/p:cNvPr',NS)
    if nv is not None and nv.get('name')=='title logo mask':
     mask=copy.deepcopy(node);mask.find('./p:nvSpPr/p:cNvPr',NS).set('id','9999');tree.append(mask)
   base[bp[n-1]]=E.tostring(b,xml_declaration=True,encoding='UTF-8',standalone=True)
  sr=E.fromstring(src[rels(sp[n-1])]);br=E.fromstring(base[rels(bp[n-1])])
  sn=next((r for r in sr if kind(r)=='notesSlide'),None);bn=next((r for r in br if kind(r)=='notesSlide'),None)
  if sn is not None and bn is not None:base[resolve(bp[n-1],bn.get('Target'))]=src[resolve(sp[n-1],sn.get('Target'))]
# Remove the oversized leftover title-page logo, preserving small master branding.
t=E.fromstring(base[bp[0]])
for pic in t.findall('.//p:pic',NS):
 off=pic.find('./p:spPr/a:xfrm/a:off',NS);ext=pic.find('./p:spPr/a:xfrm/a:ext',NS)
 if off is not None and ext is not None and int(off.get('y'))>4000000 and int(off.get('x'))>5000000 and int(ext.get('cx'))>3000000:pic.getparent().remove(pic)
base[bp[0]]=E.tostring(t,xml_declaration=True,encoding='UTF-8',standalone=True)
# The newer-model question slide must show only its matching recording, not the old circuit.
t=E.fromstring(base[bp[73]]);rr=E.fromstring(base[rels(bp[73])]);targets={r.get('Id'):r.get('Target') for r in rr}
for node in list(t.find('./p:cSld/p:spTree',NS)):
 if node.tag.endswith('}pic'):
  blip=node.find('.//a:blip',NS);target=targets.get(blip.get('{'+NS['r']+'}embed'),'') if blip is not None else ''
  if target.lower().endswith('.gif'):
   xf=node.find('./p:spPr/a:xfrm',NS)
   xf.find('a:off',NS).set('x','2200000');xf.find('a:off',NS).set('y','1100000')
   xf.find('a:ext',NS).set('cx','4900000');xf.find('a:ext',NS).set('cy','3300000')
  else:node.getparent().remove(node)
 elif node.tag.endswith('}sp'):
  off=node.find('./p:spPr/a:xfrm/a:off',NS)
  if off is not None and 900000<int(off.get('y'))<4500000:node.getparent().remove(node)
base[bp[73]]=E.tostring(t,xml_declaration=True,encoding='UTF-8',standalone=True)
# Keep the original master, native charts/workbooks, animation settings, movies, and hide states.
# New media need content types matching their copied extensions.
CT='http://schemas.openxmlformats.org/package/2006/content-types'
ct=E.fromstring(base['[Content_Types].xml']); seen={x.get('Extension') for x in ct if x.tag.endswith('Default')}
for ext,mime in {'png':'image/png','jpg':'image/jpeg','jpeg':'image/jpeg','gif':'image/gif','svg':'image/svg+xml'}.items():
 if ext not in seen:E.SubElement(ct,'{'+CT+'}Default',Extension=ext,ContentType=mime)
base['[Content_Types].xml']=E.tostring(ct,xml_declaration=True,encoding='UTF-8',standalone=True)
replacements={'transmission model':'joint-torque model','transmission path':'BPA and artificial-tendon path','transmission mechanics':'compliance and bending','Joint transmission':'Joint torque','loaded joint transmission':'measured joint torque'}
hits=[]
for part,data in list(base.items()):
 if part.startswith(('ppt/slides/slide','ppt/notesSlides/notesSlide')) and part.endswith('.xml'):
  root=E.fromstring(data);edited=False
  for t in root.findall('.//a:t',NS):
   if t.text:
    for a,b in replacements.items():
     if a.lower() in t.text.lower():t.text=re.sub(re.escape(a),b,t.text,flags=re.I);edited=True;hits.append(part)
  if edited:base[part]=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
# Multiple uses of the same original figure share one unchanged media part.
seen={};duplicates={}
for part,data in base.items():
 if part.startswith('ppt/media/'):
  key=(Path(part).suffix.lower(),hashlib.sha256(data).hexdigest())
  if key in seen:duplicates[part]=seen[key]
  else:seen[key]=part
for rp,data in list(base.items()):
 if not rp.endswith('.rels'):continue
 owner=rp.replace('/_rels/','/').removesuffix('.rels');root=E.fromstring(data);edited=False
 if rp=='_rels/.rels':owner=''
 for r in root:
  if r.get('TargetMode')=='External':continue
  target=resolve(owner,r.get('Target'))
  if target in duplicates:r.set('Target',ctx['relative_target'](owner,duplicates[target]));edited=True
 if edited:base[rp]=E.tostring(root,xml_declaration=True,encoding='UTF-8',standalone=True)
# Drop orphaned media left by earlier revisions; retained rel targets determine ownership.
targets=set()
for rp,data in base.items():
 if not rp.endswith('.rels'):continue
 owner=rp.replace('/_rels/','/').removesuffix('.rels')
 if rp=='_rels/.rels':owner=''
 for r in E.fromstring(data):
  if r.get('TargetMode')!='External':targets.add(resolve(owner,r.get('Target')))
orphan=[p for p in base if p.startswith('ppt/media/') and p not in targets]
for p in orphan:del base[p]
out=W/'build/candidate.pptx'
with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
 for p,data in base.items():z.writestr(p,data)
(W/'build/package_report.json').write_text(json.dumps({'changed_slides':sorted(changed),'native_preserved':sorted(protected),'rejected_wording_removed_from':hits,'orphan_media_removed':len(orphan),'bytes':out.stat().st_size},indent=2))
print(out, out.stat().st_size)
