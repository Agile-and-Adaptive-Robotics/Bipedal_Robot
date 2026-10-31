"""Rebuild the defense narrative in the original C3NS package, retaining evidence."""
from pathlib import Path
from lxml import etree as E
import copy,json,zipfile,subprocess,math,io,hashlib,re
from PIL import Image
from openpyxl import Workbook
W=Path(__file__).resolve().parent;ROOT=W.parents[4];DISS=W.parents[1]
BASE=W.parent/'C3NS_Revision_20261030_v16/output/Bolen_PhD_Defense_RoseGreen_20261030_v16.pptx'
BUILD=W/'build';ASSETS=W/'assets';OUT=W/'output'
for p in [BUILD,ASSETS,OUT]:p.mkdir(parents=True,exist_ok=True)
P='http://schemas.openxmlformats.org/presentationml/2006/main';A='http://schemas.openxmlformats.org/drawingml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships';PK='http://schemas.openxmlformats.org/package/2006/relationships';CT='http://schemas.openxmlformats.org/package/2006/content-types';C='http://schemas.openxmlformats.org/drawingml/2006/chart';N={'p':P,'a':A,'r':R,'c':C};EMU=914400
with zipfile.ZipFile(BASE) as z:parts={n:z.read(n) for n in z.namelist()}
original=dict(parts);manifest=[];newids=[];uid=100;nextslide=87;ct=E.fromstring(parts['[Content_Types].xml']);pres=E.fromstring(parts['ppt/presentation.xml']);pr=E.fromstring(parts['ppt/_rels/presentation.xml.rels'])
def read(n):return E.fromstring(parts[n])
def put(n,r):parts[n]=E.tostring(r,xml_declaration=True,encoding='UTF-8',standalone=True)
def el(tag,parent=None,**attrs):
 n=E.Element(tag,**{k:str(v) for k,v in attrs.items()})
 if parent is not None:parent.append(n)
 return n
def ident(s):
 n=s.find('.//p:cNvPr',N);return n.get('id') if n is not None else None
def strtext(s):return ' '.join(s.xpath('.//a:t/text()',namespaces=N))
def settext(s,t):
 b=s.find('p:txBody',N);p0=b.find('a:p',N);pp=p0.find('a:pPr',N);rp=p0.find('a:r/a:rPr',N)
 for q in list(b):
  if q.tag==f'{{{A}}}p':b.remove(q)
 q=el(f'{{{A}}}p',b)
 if pp is not None:q.append(copy.deepcopy(pp))
 rr=el(f'{{{A}}}r',q)
 if rp is not None:rr.append(copy.deepcopy(rp))
 el(f'{{{A}}}t',rr).text=t
def bounds(s,x,y,w,h):
 tr=s.find('.//a:xfrm',N)
 if tr is None:tr=el(f'{{{A}}}xfrm',s.find('p:spPr',N))
 for tag,attrs in [('off',dict(x=int(x*EMU),y=int(y*EMU))),('ext',dict(cx=int(w*EMU),cy=int(h*EMU)))]:
  z=tr.find('a:'+tag,N)
  if z is None:z=el(f'{{{A}}}'+tag,tr)
  z.attrib.update({k:str(v) for k,v in attrs.items()})
def font(s,size,color='174D24',bold=False):
 for rp in s.findall('.//a:rPr',N)+s.findall('.//a:defRPr',N):
  rp.set('sz',str(int(size*100)));rp.set('b','1' if bold else '0')
  for x in list(rp):
   if x.tag in [f'{{{A}}}solidFill',f'{{{A}}}latin']:rp.remove(x)
  fill=E.Element(f'{{{A}}}solidFill');el(f'{{{A}}}srgbClr',fill,val=color);rp.insert(0,fill);latin=E.Element(f'{{{A}}}latin',typeface='Calibri');rp.insert(1,latin)
def addtext(tree,t,x,y,w,h,size=18,color='174D24',bold=False,align='l'):
 global uid;uid+=1;s=copy.deepcopy(read('ppt/slides/slide70.xml').find('p:cSld/p:spTree/p:sp',N));s.find('p:nvSpPr/p:cNvPr',N).attrib.update({'id':str(uid),'name':t[:70]});settext(s,t);bounds(s,x,y,w,h);font(s,size,color,bold)
 b=s.find('p:txBody/a:bodyPr',N);b.attrib.update({'lIns':'0','rIns':'0','tIns':'0','bIns':'0','anchor':'t'})
 for p in s.findall('.//a:pPr',N):p.set('algn',align)
 tree.append(s);return s
def line(tree,x1,y1,x2,y2,color='174D24',width=2,arrow=False):
 global uid;uid+=1;s=el(f'{{{P}}}sp');nv=el(f'{{{P}}}nvSpPr',s);el(f'{{{P}}}cNvPr',nv,id=uid,name='Connector');el(f'{{{P}}}cNvSpPr',nv);el(f'{{{P}}}nvPr',nv);sp=el(f'{{{P}}}spPr',s);tr=el(f'{{{A}}}xfrm',sp)
 if x2<x1:tr.set('flipH','1')
 if y2<y1:tr.set('flipV','1')
 el(f'{{{A}}}off',tr,x=int(min(x1,x2)*EMU),y=int(min(y1,y2)*EMU));el(f'{{{A}}}ext',tr,cx=max(1,int(abs(x2-x1)*EMU)),cy=max(1,int(abs(y2-y1)*EMU)));el(f'{{{A}}}avLst',el(f'{{{A}}}prstGeom',sp,prst='line'));ln=el(f'{{{A}}}ln',sp,w=int(width*12700));el(f'{{{A}}}srgbClr',el(f'{{{A}}}solidFill',ln),val=color)
 if arrow:el(f'{{{A}}}tailEnd',ln,type='triangle',w='med',len='med')
 tree.append(s);return s
def box(tree,t,x,y,w,h,size=18,accent=False):
 s=addtext(tree,t,x,y,w,h,size,bold=True,align='ctr');sp=s.find('p:spPr',N)
 for z in list(sp):
  if z.tag in [f'{{{A}}}solidFill',f'{{{A}}}noFill',f'{{{A}}}ln']:sp.remove(z)
 el(f'{{{A}}}srgbClr',el(f'{{{A}}}solidFill',sp),val='E5F3E4' if accent else 'FFFFFF');el(f'{{{A}}}srgbClr',el(f'{{{A}}}ln',sp,w='25400'),val='174D24');s.find('p:txBody/a:bodyPr',N).set('anchor','ctr');return s
def image(s,asset,x,y,w,h,crop=None):
 global uid;uid+=1;data=asset.read_bytes() if isinstance(asset,Path) else parts[asset];ext=asset.suffix.lstrip('.') if isinstance(asset,Path) else asset.rsplit('.',1)[-1];name=asset.rsplit('/',1)[-1] if isinstance(asset,str) and asset.startswith('ppt/media/') else f'v17_{hashlib.sha256(data).hexdigest()[:16]}.{ext}';parts['ppt/media/'+name]=data;rid='rIdImage'+str(uid);el(f'{{{PK}}}Relationship',s['rels'],Id=rid,Type=R+'/image',Target='../media/'+name)
 if not any(x.get('Extension')==ext for x in ct if x.tag.endswith('Default')):el(f'{{{CT}}}Default',ct,Extension=ext,ContentType='image/jpeg' if ext in ['jpg','jpeg'] else 'image/'+ext)
 if crop is None:
  try:iw,ih=Image.open(io.BytesIO(data)).size;ratio=min(w/iw,h/ih);dw,dh=iw*ratio,ih*ratio;x+=(w-dw)/2;y+=(h-dh)/2;w,h=dw,dh
  except Exception:pass
 pic=el(f'{{{P}}}pic');nv=el(f'{{{P}}}nvPicPr',pic);el(f'{{{P}}}cNvPr',nv,id=uid,name=str(asset)[-80:],descr=s['title']);el(f'{{{A}}}picLocks',el(f'{{{P}}}cNvPicPr',nv),noChangeAspect='1');el(f'{{{P}}}nvPr',nv);fill=el(f'{{{P}}}blipFill',pic);el(f'{{{A}}}blip',fill,**{f'{{{R}}}embed':rid})
 if crop:el(f'{{{A}}}srcRect',fill,**{k:str(int(v*100000)) for k,v in zip(['l','t','r','b'],crop)})
 el(f'{{{A}}}fillRect',el(f'{{{A}}}stretch',fill));sp=el(f'{{{P}}}spPr',pic);tr=el(f'{{{A}}}xfrm',sp);el(f'{{{A}}}off',tr,x=int(x*EMU),y=int(y*EMU));el(f'{{{A}}}ext',tr,cx=int(w*EMU),cy=int(h*EMU));el(f'{{{A}}}avLst',el(f'{{{A}}}prstGeom',sp,prst='rect'));s['tree'].append(pic);return pic
def media(slide,shape):
 r=E.fromstring(original[f'ppt/slides/slide{slide}.xml']);pic=next(x for x in r.findall('p:cSld/p:spTree/*',N) if ident(x)==str(shape));b=pic.find('.//a:blip',N);rels=E.fromstring(original[f'ppt/slides/_rels/slide{slide}.xml.rels']);target=next(x.get('Target') for x in rels if x.get('Id')==b.get(f'{{{R}}}embed'));return 'ppt/'+target.replace('../','')
def pdf(f):
 p=DISS/'Figures'/f;dest=ASSETS/p.stem
 if not dest.with_suffix('.png').exists():subprocess.run(['C:/Users/Ben Bolen/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe','-singlefile','-r','220','-png',str(p),str(dest)],capture_output=True,check=True)
 return dest.with_suffix('.png')
def note(no,t):
 r=read('ppt/notesSlides/notesSlide70.xml')
 for sp in r.findall('.//p:sp',N):
  ph=sp.find('p:nvSpPr/p:nvPr/p:ph',N)
  if ph is not None and ph.get('type')=='body':settext(sp,t)
 put(f'ppt/notesSlides/notesSlide{no}.xml',r)
 rels=read('ppt/notesSlides/_rels/notesSlide70.xml.rels')
 for x in rels:
  if x.get('Type').endswith('/slide'):x.set('Target',f'../slides/slide{no}.xml')
 put(f'ppt/notesSlides/_rels/notesSlide{no}.xml.rels',rels)
def start(title,credit='',notes='',old=None):
 global nextslide
 if old is None:
  no=nextslide;nextslide+=1;r=read('ppt/slides/slide70.xml');tree=r.find('p:cSld/p:spTree',N)
  for x in list(tree):
   if x.tag in [f'{{{P}}}nvGrpSpPr',f'{{{P}}}grpSpPr']:continue
   if ident(x) not in ['7','8','18']:tree.remove(x)
  rels=read('ppt/slides/_rels/slide70.xml.rels')
  for x in list(rels):
   if not x.get('Type').endswith(('/slideLayout','/notesSlide')):rels.remove(x)
   elif x.get('Type').endswith('/notesSlide'):x.set('Target',f'../notesSlides/notesSlide{no}.xml')
  for part,typ in [('slides','slide'),('notesSlides','notesSlide')]:el(f'{{{CT}}}Override',ct,PartName=f'/ppt/{part}/{typ}{no}.xml',ContentType=f'application/vnd.openxmlformats-officedocument.presentationml.{typ}+xml')
  el(f'{{{PK}}}Relationship',pr,Id='rIdNew'+str(no),Type=R+'/slide',Target=f'slides/slide{no}.xml');newids.append(no)
 else:no=old;r=read(f'ppt/slides/slide{no}.xml');tree=r.find('p:cSld/p:spTree',N);rels=read(f'ppt/slides/_rels/slide{no}.xml.rels')
 s=dict(no=no,root=r,tree=tree,rels=rels,title=title,notes=notes,credit=credit)
 texts=[x for x in tree.findall('p:sp',N) if x.find('p:txBody',N) is not None]
 for x in texts:
  off=x.find('.//a:xfrm/a:off',N)
  if off is not None and int(off.get('y'))<300000 and len(strtext(x))>10:settext(x,title)
  elif off is not None and 4500000<int(off.get('y'))<4800000 and credit:settext(x,credit)
 return s
def finish(s,section,seconds=30):
 if s['no']==76:
  for sh in s['tree'].findall('p:sp',N):
   if strtext(sh).startswith('Shared effective stiffness pair:'):settext(sh,'One effective stiffness pair is shared across both configurations')
 if s['no']==13:
  cover=box(s['tree'],'',.56,2.48,.25,.22,10)
  for ln in cover.findall('p:spPr/a:ln',N):cover.find('p:spPr',N).remove(ln)
 if s['no']==5:
  clear(s);bpa_pair_group(s);addtext(s['tree'],'BPA length  lₘ  |  Initial length is measured clamp to clamp',.7,4.55,8.8,.4,19)
 if s['title']=='The manufacturer tool overpredicts force in short BPAs':
  clear(s);manufacturer_chart(s);addtext(s['tree'],'10 mm BPA  |  112 mm initial length  |  620 kPa',.8,4.6,8.4,.4,18);s['notes']='The manufacturer predicts 498.6 N for the 112 mm actuator, while the measured force is 350.9 N: a 42.1% overprediction. This is the concrete modeling problem before introducing the new fit. Source: published BPA Discussion, Effects of Resting Length.'
 if s['title']=='A model fitted at one length need not predict another length':
  clear(s)
  for x,t in [(.7,'Parameters fitted at one initial length'),(5.55,'A different initial length in the robot')]:box(s['tree'],t,x,1.65,3.7,1.1,22)
  line(s['tree'],4.5,2.2,5.35,2.2,arrow=True);addtext(s['tree'],'Accuracy across lengths is a separate question',.8,3.65,8.4,.8,26,bold=True,align='ctr')
 if s['title']=='The force model separates initial-length scale from operating state':
  clear(s);s['title']='A geometry-specific model can fail at a different initial length'
  for sh in s['tree'].findall('p:sp',N):
   if ident(sh)=='7':settext(sh,s['title'])
  image(s,pdf('20-muscles/05_Comparison-eps-converted-to.pdf'),.55,1.15,8.9,3.7,crop=(0,.5,0,0))
  s['notes']='The published comparison tests the existing models and the length-aware model against measured forces at 300 and 450 mm resting lengths for the 20 mm BPA. A model that is accurate near its identified geometry can become inaccurate at another length. The length-aware model provides more consistent predictions in these displayed cases. All plotted values and colors are from the published Figure 5; no model is refitted here. Source: BPA model-comparison section.'
 for sh in s['tree'].findall('p:sp',N):
  if not (ident(sh) or '').isdigit() or int(ident(sh))<100:continue
  tr=sh.find('p:spPr/a:xfrm',N)
  if tr is None:continue
  of=tr.find('a:off',N);ex=tr.find('a:ext',N);yy=int(of.get('y'))/EMU;hh=int(ex.get('cy'))/EMU
  if yy>=4.0 and yy<4.95 and yy+hh>4.95:of.set('y',str(int(4.45*EMU)));ex.set('cy',str(int(.45*EMU)))
 put(f"ppt/slides/slide{s['no']}.xml",s['root']);put(f"ppt/slides/_rels/slide{s['no']}.xml.rels",s['rels'])
 if s['notes']:note(s['no'],s['notes'])
 manifest.append(dict(source=s['no'],title=s['title'],section=section,seconds=seconds,credit=s['credit']))
def clear(s):
 for x in list(s['tree']):
  if x.tag in [f'{{{P}}}nvGrpSpPr',f'{{{P}}}grpSpPr']:continue
  off=x.find('.//a:xfrm/a:off',N)
  if off is not None and (int(off.get('y'))<300000 or int(off.get('y'))>4500000):continue
  s['tree'].remove(x)
def old(no,title,section,notes='',credit=''):
 s=start(title,credit,notes,old=no);finish(s,section);return s
def asset_photo(f):return DISS/'Figures/15-background'/f
def bpa_pair_group(s):
 src=E.fromstring(original['ppt/slides/slide5.xml']);items=[]
 for sh in src.findall('p:cSld/p:spTree/*',N):
  tr=sh.find('.//a:xfrm',N)
  if tr is None:continue
  of=tr.find('a:off',N);ex=tr.find('a:ext',N)
  if of is not None and int(of.get('x'))>4400000 and int(of.get('y'))<4100000:items.append(copy.deepcopy(sh))
 boxes=[(int(z.find('.//a:xfrm/a:off',N).get('x')),int(z.find('.//a:xfrm/a:off',N).get('y')),int(z.find('.//a:xfrm/a:ext',N).get('cx')),int(z.find('.//a:xfrm/a:ext',N).get('cy'))) for z in items];x=min(z[0] for z in boxes);y=min(z[1] for z in boxes);w=max(z[0]+z[2] for z in boxes)-x;h=max(z[1]+z[3] for z in boxes)-y;scale=min(8*EMU/w,3.3*EMU/h)
 g=el(f'{{{P}}}grpSp');nv=el(f'{{{P}}}nvGrpSpPr',g);el(f'{{{P}}}cNvPr',nv,id='990',name='BPA length definitions');el(f'{{{P}}}cNvGrpSpPr',nv);el(f'{{{P}}}nvPr',nv);tr=el(f'{{{A}}}xfrm',el(f'{{{P}}}grpSpPr',g));el(f'{{{A}}}off',tr,x=int((10-w*scale/EMU)/2*EMU),y=int(1.15*EMU));el(f'{{{A}}}ext',tr,cx=int(w*scale),cy=int(h*scale));el(f'{{{A}}}chOff',tr,x=x,y=y);el(f'{{{A}}}chExt',tr,cx=w,cy=h)
 for sh in items:
  for nv in sh.findall('.//p:cNvPr',N):nv.set('id',str(500+int(nv.get('id'))))
  if strtext(sh)=='Uncontracted':settext(sh,'Initial length')
  g.append(sh)
 s['tree'].append(g)
def callout(s,t,x,y,w=2.4):addtext(s['tree'],t,x,y,w,.48,15,bold=True)
def loop(s,x=3.7,y=1.65):
 box(s['tree'],'Neural controller',x,y,2.0,.7,18);box(s['tree'],'Actuator + body',x+3,y,2.0,.7,18);line(s['tree'],x+2.05,y+.35,x+2.9,y+.35,arrow=True);addtext(s['tree'],'Motor output',x+2.0,y-.38,1.1,.3,12)
 line(s['tree'],x+4,y+.78,x+4,y+1.65);line(s['tree'],x+4,y+1.65,x+1,y+1.65,arrow=True);line(s['tree'],x+1,y+1.65,x+1,y+.78);addtext(s['tree'],'Length, force, and contact feedback',x+1.35,y+1.78,3.9,.55,17)
def unified(s,focus=None):
 image(s,asset_photo('aarl_bpa_biped.png'),.5,1.15,1.8,3.2)
 image(s,media(5,5),2.8,1.35,2.6,1.4);image(s,media(10,24),2.8,3.05,2.6,1.3)
 box(s['tree'],'Rhythm',6.05,1.75,1.25,.65,17,focus=='neural');box(s['tree'],'Motor\noutput',8.05,1.75,1.25,.65,17,focus=='neural');line(s['tree'],7.4,2.08,7.95,2.08,arrow=True);line(s['tree'],8.65,2.5,8.65,3.45);line(s['tree'],8.65,3.45,6.65,3.45,arrow=True);line(s['tree'],6.65,3.45,6.65,2.5)
 addtext(s['tree'],'Sensory feedback',6.3,3.65,3,.45,17,bold=True)
 for t,x,y,w in [('Humanoid skeleton',.35,4.8,2.25),('Actuator force',2.8,2.72,2.6),('Joint torque',2.8,4.65,2.6),('Neural control',6.05,1.1,3.25)]:addtext(s['tree'],t,x,y,w,.35,16,bold=True)
 line(s['tree'],2.2,2.6,2.7,2.05,arrow=True);line(s['tree'],5.45,2.05,5.95,2.05,arrow=True)
 if focus:box(s['tree'],{'force':'Focus: actuator force','torque':'Focus: joint torque','neural':'Focus: neural control'}[focus],6.05,4.35,3.2,.55,18,True)
def divider(title,focus,section,notes):
 s=start(title,'Dissertation research pathway',notes);unified(s,focus);finish(s,section,12)
def pairfigure(title,f1,f2,section,credit,notes,labels=None):
 s=start(title,credit,notes);image(s,f1,.45,1.35,4.4,3.45);image(s,f2,5.15,1.35,4.4,3.45)
 if labels:
  for t,x in zip(labels,[.5,5.2]):addtext(s['tree'],t,x,1.05,4.2,.4,18,bold=True)
 finish(s,section);return s

# Native scatter charts retain the original human-reference color and editable data.
def human_chart(s,txt,col,label,x,y,w,h,index):
 rows=[list(map(float,l.split())) for l in (ROOT/txt).read_text().splitlines()[7:] if l.strip()];xs=[r[1] for r in rows];ys=[r[col] for r in rows]
 book=Workbook();ws=book.active;ws.title='Reference';ws.append(['Knee angle (deg)',label]);[ws.append([a,b]) for a,b in zip(xs,ys)];buf=io.BytesIO();book.save(buf);bn=f'Microsoft_Excel_Sheet{index}.xlsx';parts['ppt/embeddings/'+bn]=buf.getvalue()
 with zipfile.ZipFile(io.BytesIO(buf.getvalue())) as z:sheet=E.fromstring(z.read('xl/worksheets/sheet1.xml'));literal={c.get('r'):c.find('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v').text for c in sheet.findall('.//{http://schemas.openxmlformats.org/spreadsheetml/2006/main}c') if c.find('{http://schemas.openxmlformats.org/spreadsheetml/2006/main}v') is not None}
 cs=E.Element(f'{{{C}}}chartSpace',nsmap={'c':C,'a':A,'r':R});el(f'{{{C}}}date1904',cs,val='0');chart=el(f'{{{C}}}chart',cs);plot=el(f'{{{C}}}plotArea',chart);el(f'{{{C}}}layout',plot);scatter=el(f'{{{C}}}scatterChart',plot);el(f'{{{C}}}scatterStyle',scatter,val='line');ser=el(f'{{{C}}}ser',scatter);el(f'{{{C}}}idx',ser,val='0');el(f'{{{C}}}order',ser,val='0');el(f'{{{C}}}v',el(f'{{{C}}}tx',ser)).text=label
 sp=el(f'{{{C}}}spPr',ser);ln=el(f'{{{A}}}ln',sp,w='31750');el(f'{{{A}}}srgbClr',el(f'{{{A}}}solidFill',ln),val='000000');el(f'{{{A}}}prstDash',ln,val='sysDot');el(f'{{{C}}}symbol',el(f'{{{C}}}marker',ser),val='none')
 for tag,vals,letter in [('xVal',xs,'A'),('yVal',ys,'B')]:
  ref=el(f'{{{C}}}numRef',el(f'{{{C}}}'+tag,ser));el(f'{{{C}}}f',ref).text=f'Reference!${letter}$2:${letter}${len(vals)+1}';cache=el(f'{{{C}}}numCache',ref);el(f'{{{C}}}formatCode',cache).text='General';el(f'{{{C}}}ptCount',cache,val=len(vals))
  for i,v in enumerate(vals):el(f'{{{C}}}v',el(f'{{{C}}}pt',cache,idx=i)).text=literal[f'{letter}{i+2}']
 for aid in [1000+index,2000+index]:el(f'{{{C}}}axId',scatter,val=aid)
 for aid,other,pos,mi,ma,title in [(1000+index,2000+index,'b',-120,10,'Knee angle (deg)'),(2000+index,1000+index,'l',math.floor(min(ys)/10)*10,math.ceil(max(ys)/10)*10,'Torque (N·m)')]:
  ax=el(f'{{{C}}}valAx',plot);el(f'{{{C}}}axId',ax,val=aid);sc=el(f'{{{C}}}scaling',ax);el(f'{{{C}}}orientation',sc,val='minMax');el(f'{{{C}}}min',sc,val=mi);el(f'{{{C}}}max',sc,val=ma);el(f'{{{C}}}delete',ax,val='0');el(f'{{{C}}}axPos',ax,val=pos);el(f'{{{C}}}numFmt',ax,formatCode='0',sourceLinked='0');el(f'{{{C}}}tickLblPos',ax,val='low');el(f'{{{C}}}crossAx',ax,val=other);el(f'{{{C}}}crosses',ax,val='min');el(f'{{{C}}}crossBetween',ax,val='midCat')
  tx=el(f'{{{C}}}txPr',ax);el(f'{{{A}}}bodyPr',tx);el(f'{{{A}}}lstStyle',tx);pa=el(f'{{{A}}}p',tx);rp=el(f'{{{A}}}defRPr',el(f'{{{A}}}pPr',pa),sz='1400');el(f'{{{A}}}srgbClr',el(f'{{{A}}}solidFill',rp),val='333333');el(f'{{{A}}}latin',rp,typeface='Arial')
 el(f'{{{C}}}plotVisOnly',chart,val='1');el(f'{{{C}}}dispBlanksAs',chart,val='gap');el(f'{{{C}}}externalData',cs,**{f'{{{R}}}id':'rId1'});put(f'ppt/charts/chart{index}.xml',cs);relr=E.Element(f'{{{PK}}}Relationships',nsmap={None:PK});el(f'{{{PK}}}Relationship',relr,Id='rId1',Type=R+'/package',Target='../embeddings/'+bn);put(f'ppt/charts/_rels/chart{index}.xml.rels',relr)
 el(f'{{{CT}}}Override',ct,PartName=f'/ppt/charts/chart{index}.xml',ContentType='application/vnd.openxmlformats-officedocument.drawingml.chart+xml');rid='rIdHuman'+str(index);el(f'{{{PK}}}Relationship',s['rels'],Id=rid,Type=R+'/chart',Target=f'../charts/chart{index}.xml')
 frame=el(f'{{{P}}}graphicFrame');nv=el(f'{{{P}}}nvGraphicFramePr',frame);el(f'{{{P}}}cNvPr',nv,id=str(800+index),name=label);el(f'{{{P}}}cNvGraphicFramePr',nv);el(f'{{{P}}}nvPr',nv);tr=el(f'{{{P}}}xfrm',frame);el(f'{{{A}}}off',tr,x=int(x*EMU),y=int(y*EMU));el(f'{{{A}}}ext',tr,cx=int(w*EMU),cy=int(h*EMU));gd=el(f'{{{A}}}graphicData',el(f'{{{A}}}graphic',frame),uri=C);el(f'{{{C}}}chart',gd,**{f'{{{R}}}id':rid});s['tree'].append(frame)

def manufacturer_chart(s):
 chart=E.fromstring(original['ppt/charts/chart1.xml']);plot=chart.find('c:chart/c:plotArea',N);bars=plot.findall('c:barChart',N)
 for b in bars[1:]:plot.remove(b)
 primary={x.get('val') for x in bars[0].findall('c:axId',N)}
 for ax in list(plot):
  if ax.tag in [f'{{{C}}}valAx',f'{{{C}}}catAx'] and ax.find('c:axId',N).get('val') not in primary:plot.remove(ax)
 ser=bars[0].find('c:ser',N)
 for parent in [ser]+ser.findall('c:dLbls',N):
  for item in list(parent):
   if item.tag not in [f'{{{C}}}dPt',f'{{{C}}}dLbl']:continue
   ix=item.find('c:idx',N)
   if ix is None:continue
   if ix.get('val')=='1':parent.remove(item)
   elif ix.get('val')=='2':ix.set('val','1')
 for tag,data,col in [('cat',['Measured','Manufacturer'],'A'),('val',[350.9,498.6],'B')]:
  e=ser.find('c:'+tag,N)
  for child in list(e):e.remove(child)
  ref=el(f'{{{C}}}'+('strRef' if tag=='cat' else 'numRef'),e);el(f'{{{C}}}f',ref).text=f'Reference!${col}$2:${col}$3';cache=el(f'{{{C}}}'+('strCache' if tag=='cat' else 'numCache'),ref)
  if tag=='val':el(f'{{{C}}}formatCode',cache).text='0.0'
  el(f'{{{C}}}ptCount',cache,val='2')
  for i,v in enumerate(data):el(f'{{{C}}}v',el(f'{{{C}}}pt',cache,idx=i)).text=str(v)
 tx=ser.find('c:tx',N)
 for child in list(tx):tx.remove(child)
 el(f'{{{C}}}v',tx).text='Isometric force (N)'
 leg=chart.find('c:chart/c:legend',N)
 if leg is not None:chart.find('c:chart',N).remove(leg)
 for ex in chart.findall('c:externalData',N):ex.set(f'{{{R}}}id','rId1')
 put('ppt/charts/chart7.xml',chart);wb=Workbook();ws=wb.active;ws.title='Reference';ws.append(['Case','Force (N)']);ws.append(['Measured',350.9]);ws.append(['Manufacturer',498.6]);buf=io.BytesIO();wb.save(buf);parts['ppt/embeddings/Microsoft_Excel_Sheet7.xlsx']=buf.getvalue();rr=E.Element(f'{{{PK}}}Relationships',nsmap={None:PK});el(f'{{{PK}}}Relationship',rr,Id='rId1',Type=R+'/package',Target='../embeddings/Microsoft_Excel_Sheet7.xlsx');put('ppt/charts/_rels/chart7.xml.rels',rr);el(f'{{{CT}}}Override',ct,PartName='/ppt/charts/chart7.xml',ContentType='application/vnd.openxmlformats-officedocument.drawingml.chart+xml')
 frame=copy.deepcopy(E.fromstring(original['ppt/slides/slide8.xml']).find('p:cSld/p:spTree/p:graphicFrame',N));frame.find('p:nvGraphicFramePr/p:cNvPr',N).set('id','887');frame.find('.//c:chart',N).set(f'{{{R}}}id','rIdManufacturer');tr=frame.find('p:xfrm',N);tr.find('a:off',N).attrib.update({'x':str(int(.7*EMU)),'y':str(int(1.15*EMU))});tr.find('a:ext',N).attrib.update({'cx':str(int(8.6*EMU)),'cy':str(int(3.4*EMU))});s['tree'].append(frame);el(f'{{{PK}}}Relationship',s['rels'],Id='rIdManufacturer',Type=R+'/chart',Target='../charts/chart7.xml')

# The main-talk sequence is authored below. Original nonselected slides remain questions.

pub='Bolen et al. (2026), Actuators 15(5):230'
intro='Dissertation Ch. 1, Motivation and Research Objectives'
torque='Dissertation Ch. 4, Isometric Knee Torque'
tools='Dissertation Ch. 6, Table 6.1'

# 1-11: motivation, project, and four gaps.
s=start('Design and Modeling of Biomimetic Humanoid Robots',notes='This research develops actuator force models, joint-torque design factors, sensory-evidence resources, and neuromechanical tools for a biomimetic humanoid robot. The physical prototype on the cover establishes the long-term project. Source: dissertation Chapters 1 and 2; AARL robot photograph.',old=1)
for sh in s['tree'].findall('p:sp',N):
 if ident(sh)=='117':settext(sh,s['title']);bounds(sh,3.45,1.55,5.8,2.15);font(sh,25,'FFFFFF',True)
 elif ident(sh)=='2':bounds(sh,3.45,3.95,5.6,.65)
image(s,asset_photo('aarl_bpa_biped.png'),.85,1.65,2.0,3.56);finish(s,'Opening',15)
old(2,'Prepared floors favor wheels, while obstacles demand adaptation','Motivation','Warehouse robots can exploit predictable floors and efficient wheels. The DARPA clips illustrate the different balance and contact demands posed by obstacles. Our long-term goal is autonomous movement in less structured environments, not replacing wheels where they already work well. Sources: existing Amazon/DARPA clips and dissertation Motivation.')
s=start('Biomimetic robots let us investigate biological movement',intro,'The reciprocal goal is to use a robot as a controllable physical model. Musculoskeletal geometry and muscle-like actuation make some biological hypotheses experimentally accessible. Such models complement biological experiments; their similarity and limits must be established. Sources: dissertation Motivation; Asano et al., Kengoro; original talk slide 3.')
image(s,media(3,9),.6,1.2,4.15,3.6);image(s,asset_photo('aarl_bpa_biped.png'),6.15,1.2,2.0,3.15);addtext(s['tree'],'Biological inspiration',.7,4.8,4,.4,18,bold=True);addtext(s['tree'],'Controllable robot model',5.4,4.8,4.1,.4,18,bold=True);finish(s,'Motivation')
s=start('Compliant actuators give the body a mechanical response to load',intro,'Legs with compliant actuators are a candidate approach to uncertain contacts. The mechanical response is part of the system, rather than being supplied entirely by a joint controller. Energy storage and recovery motivate this approach, but this dissertation does not measure a whole-robot energy advantage. Source: dissertation Motivation and BPA background.')
image(s,asset_photo('festo_bpas_inflated_deflated.jpg'),.5,1.15,5.0,3.75);box(s['tree'],'Applied load',6.05,1.6,2.6,.65,20);line(s['tree'],7.35,2.4,7.35,2.9,arrow=True);box(s['tree'],'Compliant response',5.75,3.1,3.2,.75,20,True);addtext(s['tree'],'Mechanics acts within the control loop',5.65,4.25,3.6,.65,19);finish(s,'Motivation')
s=start('Sensory feedback links neural output to the mechanical environment',intro,'A neural controller receives signals created by the mechanics it controls. Muscle length, force, and contact therefore change the operating conditions of that controller. The question is how to retain this interaction while making its individual contributions testable. Sources: dissertation Chapters 1, 2, and 6.')
image(s,asset_photo('aarl_bpa_biped.png'),.6,1.2,2.0,3.6);loop(s,3.7,1.7);finish(s,'Motivation')
s=start('Robotic models allow controlled changes to mechanics and feedback',intro,'A model allows an actuator route or sensory connection to be changed while other conditions are retained. Comparable interventions in a person can be impractical or unethical. The aim is a repeatable test of a mechanism, with the mechanical model and input protocol stated. Source: dissertation Motivation.')
box(s['tree'],'Same mechanical model',.7,1.5,3.6,.8,21);box(s['tree'],'Same input protocol',5.65,1.5,3.6,.8,21);line(s['tree'],2.5,2.45,2.5,3.0,arrow=True);line(s['tree'],7.45,2.45,7.45,3.0,arrow=True);box(s['tree'],'Feedback intact',.7,3.2,3.6,.8,21);box(s['tree'],'One pathway changed',5.65,3.2,3.6,.8,21,True);addtext(s['tree'],'Compare neural and mechanical responses',1.8,4.5,6.6,.55,24,bold=True,align='ctr');finish(s,'Motivation')
pairfigure('The biology–robotics cycle can also inform assistive-device design',pdf('15-background/human_gait2392_paths.pdf'),pdf('15-background/steeleknee.pdf'),'Motivation',intro+'; Steele et al. (2018)','A tertiary goal is to inform prosthetic and orthotic design. The present contribution is the actuation, mechanics, and experimental framework, rather than a clinically validated assistive device. The human reference and biomimetic knee illustrate the connection. Sources: dissertation Motivation; Steele knee figure.',labels=['Biological mechanics','Biomimetic mechanism'])
s=start('The humanoid project couples actuators, joints, and neural control',intro,'This is the common picture for the talk. The skeleton sets geometry, a BPA supplies tensile force, the joint transmits it as torque, and a neural controller closes the sensory loop. The next sections examine those connections at increasing levels of integration. Sources: dissertation Chapters 1, 3, 4, and 6.');unified(s);finish(s,'Motivation')
s=start('Actuator selection requires force predictions across initial lengths',intro,'The short actuators needed for a robot leg did not produce the force expected from the manufacturer specification. A model fitted at one resting length is therefore insufficient for selecting another length. This is the first gap addressed experimentally. Source: dissertation Motivation and published BPA characterization.')
image(s,media(5,5),.6,1.55,4.2,2.5);box(s['tree'],'Initial length',5.5,1.45,3.3,.65,21);box(s['tree'],'Pressure',5.5,2.45,3.3,.65,21);box(s['tree'],'Contraction',5.5,3.45,3.3,.65,21);addtext(s['tree'],'All three affect available force',.65,4.5,4.3,.6,21,bold=True);finish(s,'Motivation')
s=start('Joint torque also depends on the loaded transmission path',intro,'Force must reach the skeleton through fittings, brackets, tendons, and a path that can bend around the knee. The rigid-body route model does not include all those effects. The second gap is to identify which corrections are needed and test them on more than one configuration. Source: dissertation Motivation and Isometric Knee Torque.')
image(s,media(10,24),.65,1.2,5.0,3.65);box(s['tree'],'BPA force',6.15,1.5,2.7,.65,21);line(s['tree'],7.5,2.35,7.5,2.85,arrow=True);box(s['tree'],'Transmission path',5.8,3.0,3.4,.65,21,True);line(s['tree'],7.5,3.85,7.5,4.2,arrow=True);addtext(s['tree'],'Joint torque',6.1,4.35,3,.5,23,bold=True);finish(s,'Motivation')
s=start('Four knowledge and tool gaps guide this dissertation',intro,'The four contributions answer distinct questions: how actuator length changes force, which mechanical losses change joint torque, where to find sensory evidence with its experimental context, and how to implement and test coupled neural models. SADb addresses evidence access; it is presented in the question section. Source: dissertation Research Objectives.',old=4);clear(s)
for i,(left,right) in enumerate([('Force across actuator lengths','BPA characterization'),('Force transmitted to joint torque','Joint design factors'),('Context for sensory pathways','SADb'),('Circuit behavior across bodies','SNS modeling tools')]):
 y=1.3+i*.9;addtext(s['tree'],left,.55,y,4.2,.55,21,bold=True);line(s['tree'],4.85,y+.25,5.35,y+.25,arrow=True);addtext(s['tree'],right,5.6,y,3.8,.55,21)
finish(s,'Motivation')

# 12-27: actuator operation, model limitations, experiment, and published results.
divider('Braided Pneumatic Actuator Modeling','force','BPA modeling','The first section measures how much force is available before joint geometry and neural control are added. The focus is the BPA itself. Source: dissertation Chapter 3; Bolen et al. (2026).')
s=start('Pressurization shortens the braided actuator and produces tension',pub,'An inflatable bladder expands inside a braided sleeve. The braid converts radial expansion into axial shortening. A BPA pulls in tension; antagonistic paths are needed for bidirectional joint action. Initial length l0 is the unpressurized clamp-to-clamp length; l620 is the unloaded length at 620 kPa, and lm is current length. Source: published BPA Methods and background.',old=5);clear(s);image(s,media(5,5),1.15,1.55,6.8,2.4);addtext(s['tree'],'Pressurized: shorter and wider',1.2,1.1,6.9,.4,20,bold=True);addtext(s['tree'],'Initial, unpressurized length  l₀',1.2,4.2,6.9,.45,20);addtext(s['tree'],'Current BPA length  lₘ',1.2,4.75,6.9,.4,17);line(s['tree'],7.6,2.75,9.0,2.75,arrow=True);addtext(s['tree'],'Tension',8.0,3.05,1.4,.4,17);finish(s,'BPA modeling')
s=start('Antagonistic BPA paths can rotate a biomimetic joint',pub+'; Dissertation Ch. 2','The actuator is useful because its tensile force can act through a moment arm about a joint. Opposed muscle paths supply the two torque directions. The force available along either path depends on its current length and pressure. Source: dissertation BPA and knee background.')
image(s,media(10,24),.5,1.25,5.5,3.7);image(s,asset_photo('aarl_bpa_biped.png'),6.65,1.3,1.9,3.45);addtext(s['tree'],'Muscle-like actuation must be sized for the route',.6,4.8,8.7,.45,20,bold=True);finish(s,'BPA modeling')
s=start('A useful force model must describe pressure, contraction, and length',pub,'The design problem is to choose a resting length and determine force over the operating pressure and contraction range. A single maximum-force rating cannot answer that question. The following comparisons show why the existing models needed further characterization. Source: published BPA model and dissertation background.')
for t,x in [('Initial length',.6),('Pressure',3.65),('Contraction',6.7)]:box(s['tree'],t,x,1.55,2.7,.8,22)
for x in [1.95,5.0,8.05]:line(s['tree'],x,2.55,5.0,3.15,arrow=True)
box(s['tree'],'Available isometric force',2.1,3.4,5.8,.95,25,True);finish(s,'BPA modeling')
s=start('The manufacturer tool overpredicts force in short BPAs',pub,'The dashed manufacturer curve predicts more force as resting length becomes very short. The measured points show the opposite trend. This discrepancy motivated a length-dependent characterization, especially for compact robot actuators. Do not yet attribute the improvement to the new fit; first identify the mismatch. Source: published Figure 2a.')
image(s,pdf('20-muscles/02a_MaxForce10-eps-converted-to.pdf'),.8,1.1,7.9,3.8);addtext(s['tree'],'Largest discrepancy at short initial lengths',5.45,1.45,3.4,.8,19,bold=True);finish(s,'BPA modeling')
s=start('A model fitted at one length need not predict another length',pub+'; Dissertation Ch. 2, BPA force models','Sárosi and Martens–Boblan provide useful models at specified actuator geometries. Selecting a different resting length is a separate generalization question. The published comparison here evaluates existing parameter sets against different tested lengths. Source: dissertation BPA background and published model comparison; see Figure 5 for the quantitative comparisons.')
for x,lab,length in [(.7,'Published parameter set','One actuator length'),(5.55,'Robot design choice','A different length')]:
 image(s,media(5,5),x,1.55,3.7,1.7);addtext(s['tree'],lab,x,1.1,3.7,.5,20,bold=True);addtext(s['tree'],length,x,3.55,3.7,.55,20)
line(s['tree'],4.5,2.35,5.35,2.35,arrow=True);addtext(s['tree'],'Prediction across lengths must be tested',1.5,4.5,7.0,.55,24,bold=True,align='ctr');finish(s,'BPA modeling')
s=start('A test jig measured actuator force while length was held fixed',pub,'The jig fixes muscle length while pressure and force are recorded. The load cell measures axial force, and the pressure sensor records the input. Repeating at different fixed lengths and contractions supplies the force surface needed for design. Source: published isometric-force Methods and Figure 1.')
image(s,media(5,4),.5,1.1,4.2,3.9);addtext(s['tree'],'Hold length',5.35,1.55,3.5,.5,24,bold=True);addtext(s['tree'],'Set pressure',5.35,2.65,3.5,.5,24,bold=True);addtext(s['tree'],'Measure force',5.35,3.75,3.5,.5,24,bold=True);finish(s,'BPA modeling')
s=start('Tests span two diameters and multiple operating conditions',pub,'The normalized-force dataset contains 537 points: 321 for 10 mm BPAs and 216 for 20 mm BPAs. The published fit used 80% for training and 20% for validation. Initial length, contraction, and pressure were varied. These counts refer to that normalized-force dataset, not to independent actuators. Source: published Figure 4 caption.')
for x,dia,count in [(.7,'10 mm','321'),(5.6,'20 mm','216')]:
 addtext(s['tree'],dia,x,1.35,3.7,.55,27,bold=True,align='ctr');addtext(s['tree'],count,x,2.1,3.7,1.0,48,bold=True,align='ctr');addtext(s['tree'],'force–pressure–contraction points',x,3.3,3.7,.75,19,align='ctr')
addtext(s['tree'],'80% fit data  |  20% validation data',1.5,4.45,7,.55,24,bold=True,align='ctr');finish(s,'BPA modeling')
for dia,file in [('10 mm','02a_MaxForce10-eps-converted-to.pdf'),('20 mm','02b_MaxForce20-eps-converted-to.pdf')]:
 s=start(f'Shorter {dia} BPAs deliver less maximum isometric force',pub,'The measured points establish the resting-length dependence at 620 kPa. The fitted curve follows the declining force at short lengths, while the manufacturer curve does not capture that trend. This determines the force scale available to a robot design. Source: published Figure 2, corresponding diameter panel.');image(s,pdf('20-muscles/'+file),.65,1.05,8.7,3.95);finish(s,'BPA modeling')
s=start('Pressure changes force across the tested initial lengths',pub,'The 10 mm data extend the maximum-force relationship across pressure. Each isocline corresponds to a specified pressure. Both pressure and initial length are therefore inputs to the force prediction, rather than applying a constant force rating. Source: published Figure 3a.');image(s,pdf('20-muscles/03a_MaxForce10_all-eps-converted-to.pdf'),.65,1.05,8.7,3.95);finish(s,'BPA modeling')
s=start('Maximum free contraction still needs an individual measurement',pub,'The maximum-contraction data vary substantially. The experiment did not establish a reliable initial-length relationship for this quantity, so the actuator-specific free contraction must be measured for the force normalization. The displayed uncertainty comes from the published length-measurement analysis, not a new measurement campaign. Source: published Figure 3b.');image(s,pdf('20-muscles/03b_kmax10-eps-converted-to.pdf'),.65,1.05,8.7,3.95);finish(s,'BPA modeling')
s=start('Normalization puts different operating conditions on common scales',pub,'Force is expressed relative to maximum force at 620 kPa, pressure relative to 620 kPa, and contraction relative to the measured free contraction. These definitions explain the axes of the following plots. The full equations and coefficients remain available for questions. Source: published actuator force model.')
for i,(a,b) in enumerate([('Force','Fraction of the force scale at 620 kPa'),('Pressure','Fraction of 620 kPa'),('Contraction','Fraction of measured free contraction')]):
 y=1.45+i*1.1;addtext(s['tree'],a,.65,y,2.0,.55,24,bold=True);line(s['tree'],2.8,y+.25,3.35,y+.25,arrow=True);addtext(s['tree'],b,3.7,y,5.7,.75,23)
finish(s,'BPA modeling')
s=start('Normalized force curves describe contraction across pressures',pub,'With the normalization defined, each color denotes a pressure, and the measured points can be compared with the solid fitted curves and dashed manufacturer predictions. Separate coefficients describe 10 and 20 mm actuators. Initial length supplies the dimensional force scale. Source: published Figure 4. The long equation is omitted from the main-talk visual.',old=7);clear(s);image(s,pdf('20-muscles/04a_FStar10-eps-converted-to.pdf'),.5,1.2,4.4,3.75);image(s,pdf('20-muscles/04b_FStar20-eps-converted-to.pdf'),5.1,1.2,4.4,3.75);finish(s,'BPA modeling')
old(8,'At 112 mm, the new model reduces force prediction error to 4.5%','BPA modeling','This published example compares a measured 350.9 N, the model prediction of 335.3 N, and the manufacturer prediction of 498.6 N. Absolute errors relative to the measurement are 4.5% and 42.1%. The overlaid black bars use the right percentage axis; they are prediction-error magnitudes, not statistical intervals. Source: published BPA Discussion.',pub)
s=start('The force model separates initial-length scale from operating state',pub,'The dimensional force prediction combines a length-dependent force scale with the normalized pressure–contraction surface. That makes initial length an explicit design variable over the tested range. The next question is whether this actuator force reaches the joint as the rigid route model predicts. Source: published BPA force model; dissertation Chapters 3 and 4.')
box(s['tree'],'Initial length',.6,1.45,3.6,.65,23);box(s['tree'],'Pressure + contraction',5.6,1.45,3.8,.65,23);line(s['tree'],2.4,2.3,2.4,2.8,arrow=True);line(s['tree'],7.5,2.3,7.5,2.8,arrow=True);box(s['tree'],'Force scale',.6,3.0,3.6,.65,23);box(s['tree'],'Fraction of that force',5.6,3.0,3.8,.65,23);addtext(s['tree'],'×',4.65,3.0,.65,.7,32,bold=True,align='ctr');addtext(s['tree'],'Predicted actuator force',2.05,4.35,5.9,.6,26,bold=True,align='ctr');finish(s,'BPA modeling')

# 28-48: human targets, measured mismatch, corrections, validation, and redesign.
divider('Isometric Knee Torque','torque','Joint mechanics','The actuator force model is now available. The second section asks how a route transmits that force to the knee, and how the prediction can guide a design against human-muscle torque targets. Source: dissertation Chapter 4.')
pairfigure('Human muscle paths provide a geometric reference for robot design',pdf('15-background/human_gait2392_paths.pdf'),pdf('15-background/robot_pam_routes.pdf'),'Joint mechanics','Bolen and Hunt (2019); dissertation Chs. 2 and 4','The earlier work matched robot geometry to the Gait2392 reference and compared muscle paths. The red paths show how muscle force is routed over the skeleton. The same skeletal geometry does not imply the same force–length behavior when muscles are replaced with BPAs. Sources: original talk and dissertation background.',labels=['Human muscle paths','Robot BPA routes'])
old(72,'Joint torque depends on actuator force and its moment arm','Joint mechanics','The shortest distance from the joint axis to the force line is the moment arm. Torque about that axis is force multiplied by that perpendicular distance, with direction given by the cross product. Both force and moment arm vary as the knee moves. This geometry makes the human-muscle torque envelope a meaningful design target. Sources: original talk slide 8; dissertation moment-arm methods.','Bolen and Hunt (2019); dissertation joint-torque methods')
s=start('OpenSim supplies muscle-specific torque targets across knee angle','OpenSim Gait2392 reference; Testing_Data/2022_02_Festo','These reference queries provide the knee-torque profiles of biceps femoris short head and vastus medialis. They are model-derived isometric targets, not measured maximum voluntary human torques. The stored knee-angle column uses degrees in the established evaluators. Sources: OpenSim_Bifem_Results.txt, column bifemsh_r, and OpenSim_Vasti_Results.txt, column vas_med_r; dissertation Chapter 4.')
human_chart(s,'Testing_Data/2022_02_Festo/OpenSim_Bifem_Results.txt',3,'Biceps femoris short head',.55,1.65,4.25,2.9,5);human_chart(s,'Testing_Data/2022_02_Festo/OpenSim_Vasti_Results.txt',4,'Vastus medialis',5.15,1.65,4.25,2.9,6)
for x,t in [(.65,'Biceps femoris short head'),(5.25,'Vastus medialis')]:addtext(s['tree'],t,x,1.1,4.1,.5,20,bold=True);addtext(s['tree'],'Knee angle (deg)',x+.6,4.65,3.2,.35,16)
addtext(s['tree'],'Torque (N·m)',.65,1.5,3,.3,14);addtext(s['tree'],'Torque (N·m)',5.25,1.5,3,.3,14);finish(s,'Joint mechanics')
s=start('Human attachment locations do not reproduce human torque with BPAs','Bolen and Hunt (2019), Fig. 8; dissertation Ch. 2','The earlier biceps-femoris-long-head comparison illustrates the mismatch even when a robot route is based on human geometry. The knee-torque surfaces differ because the actuator force–length relation and route geometry differ. The figure is a prior model comparison; the new measurements that follow test the joint prediction directly. Source: masters torque comparison, human left and robot right.')
image(s,pdf('15-background/masters_torque_comparison.pdf'),.65,1.65,8.7,2.95,crop=(0,.332,0,.337));addtext(s['tree'],'Human reference',.8,1.1,3.8,.45,21,bold=True);addtext(s['tree'],'BPA-actuated robot',5.2,1.1,4.2,.45,21,bold=True);finish(s,'Joint mechanics')
s=start('A modular stand tests joint geometry and actuator routes separately',torque,'The stand accommodates pinned and biomimetic knees, flexor and extensor routes, and different actuator dimensions. Load-cell measurements can be transformed to joint torque at the instantaneous knee center. This makes the mechanics experimentally testable before whole-robot controller tuning. Source: dissertation Figure 4.1 and Methods.')
image(s,pdf('30-torque/testJigs_reviewed.pdf'),.5,1.08,9.0,3.95);finish(s,'Joint mechanics')
old(9,'Improved BPA characterization does not explain the joint-torque mismatch','Joint mechanics','The four-panel pinned-flexor comparison uses the characterized force model before the mechanical correction terms. Predicted and measured torque still differ. This is the trigger for modeling length offsets, effective attachment-system compliance, and bending loss. The complete dissertation four-panel layout and the actual test photograph are retained. Source: dissertation pinned-flexor pre-optimization figure.',torque)
s=start('A constant offset changes the usable actuator length',torque,'The first correction is a constant length offset, chi0. It corrects the difference between the nominal path-based length and the actuator length available to generate force. It is distinct from force-dependent compliance. The length convention is defined before fitting the remaining mechanics. Source: Chapter 4 length-correction equation and Appendix C.')
image(s,media(5,5),.7,1.9,5.9,1.25,crop=(0,.5,0,0));addtext(s['tree'],'Initial length  l₀',.8,3.9,4.6,.5,22);line(s['tree'],1.6,3.45,5.7,3.45);line(s['tree'],1.6,3.28,1.6,3.62);line(s['tree'],5.7,3.28,5.7,3.62);line(s['tree'],5.7,3.45,6.2,3.45);line(s['tree'],6.2,3.28,6.2,3.62);addtext(s['tree'],'χ₀',5.7,3.65,.7,.65,26,bold=True);addtext(s['tree'],'Constant length offset',6.45,2.2,2.8,1.0,24,bold=True);finish(s,'Joint mechanics')
s=start('Loaded attachments change effective muscle length',torque,'The fitted chi1 and chi2 values describe effective system stiffness. They include the combined effects of brackets, fixtures, tendons, and loading hardware in these tests; they are not literal beam-stiffness measurements of a single bracket. The load-dependent displacement changes muscle length and therefore force. Source: Chapter 4 compliance model and Appendix C.')
image(s,media(10,24),.45,1.1,5.0,3.9);box(s['tree'],'Actuator force',5.8,1.4,3.4,.65,22);line(s['tree'],7.5,2.25,7.5,2.7,arrow=True);box(s['tree'],'Attachment-system deflection',5.5,2.9,3.95,.9,21,True);addtext(s['tree'],'χ₁, χ₂: effective system stiffness',5.5,4.3,4.0,.65,19,bold=True);finish(s,'Joint mechanics')
old(10,'Length offsets, compliance, and bending alter delivered torque','Joint mechanics','The complete design-factor picture distinguishes a constant length offset, effective system stiffness, and length loss due to bending. The bending term is shown along the arc R alpha, with delta-l defined in the schematic. The stiffness pair is shared across flexor and extensor tests. Source: Chapter 4 and Appendix C. These corrections are fitted mechanical descriptions, not direct measurements of individual bracket stiffness.',torque)
old(11,'The design factors improve the pinned-flexor torque prediction','Joint mechanics','The corrected prediction is compared with the same four pinned-flexor configurations shown before correction. The geometry and force characterization have not been replaced by an arbitrary torque fit. The design factors act through effective length and force transmission. The dissertation four-panel layout and marginal test photograph are retained. Source: pinned-flexor optimized figure.',torque)
old(76,'Separate tests evaluate the fitted mechanical corrections','Joint mechanics','For the adopted pinned-flexor fold, two configurations are used for fitting and two for validation; the 47 cm test is not used. For the pinned-extensor fold, three are used for fitting, three for validation, and three are not used. The same effective stiffness pair is retained across flexor and extensor configurations. Source: dissertation Appendix C adopted-fold definitions. Detailed parameter-pick identifiers remain in the question material.',torque)
old(12,'The shared stiffness pair also supports pinned-extensor predictions','Joint mechanics','The extensor comparisons test whether the corrections remain useful for a different force path, rather than choosing a separate stiffness pair for every panel. The full six-panel dissertation layout is retained with the rotated test photograph. Source: Chapter 4 pinned-extensor validation.',torque)
pairfigure('The biomimetic knee changes its moment arm as the center moves',pdf('15-background/steeleknee.pdf'),pdf('30-torque/KneeICR-eps-converted-to.pdf'),'Joint mechanics','Steele et al. (2018); dissertation biomimetic-knee Methods','The four-bar knee has a moving instantaneous center of rotation. Both the torque transformation and the actuator moment arm must respect that changing geometry. This is why validation on a pinned hinge alone would not establish the prediction for the biomimetic joint. Sources: Steele knee and Chapter 4 ICR figure.')
old(13,'The design factors improve predictions on the biomimetic knee','Joint mechanics','These flexor and extensor comparisons test the design factors on the biomimetic joint. The earlier route optimization had predicted that the required human torque would be met, but the physical route did not achieve that prediction. The route optimization therefore had to be modified to include the corrected mechanics. Source: dissertation Figures 4.10 and 4.11; the distinction between model prediction and physical validation is essential.',torque)
s=start('Corrected mechanics can be used inside the route-design search',torque,'A candidate attachment and wrapping route determines muscle length and moment arm. The characterized force model and mechanical corrections predict the torque envelope. The search then compares that envelope with the selected human-muscle target subject to geometric and actuator constraints. The identification and redesign algorithms serve different tasks; optimizer settings are question material. Source: Chapter 4 placement search and Appendix C.')
for i,t in enumerate(['Route geometry','BPA force','Mechanical corrections','Joint torque','Human target']):
 x=.5+i*1.9;box(s['tree'],t,x,1.7,1.65,.95,18,i==2)
 if i<4:line(s['tree'],x+1.7,2.18,x+1.85,2.18,arrow=True)
addtext(s['tree'],'Constrain contraction, pressure, clearance, and range of motion',.65,3.6,8.7,.8,23,bold=True,align='ctr');finish(s,'Joint mechanics')
old(14,'Route redesign increases whole-leg torque coverage','Joint mechanics','The whole-leg comparison summarizes what changes when actuator routes are redesigned rather than simply copied from human attachments. The native chart and its source workbook are retained. Coverage and range-of-motion implications are interpreted using the definitions in the dissertation; this is a design comparison, not a demonstration of autonomous walking. Source: whole-leg comparison and Chapter 4.',torque)
old(15,'All fifteen flexor measurements meet the selected human target','Joint mechanics','The two-BPA flexor is tested at fifteen knee angles. The measured points meet or exceed the selected biceps-femoris-short-head reference at those angles. Distinguish the human reference, optimized 620 kPa prediction, measured torque, and pressure-specific model prediction. Exact installed lengths and reaction-point details are retained in the question material and Appendix C. Source: final flexor physical tests.',torque)
s=start('A cam changes the extensor route through the knee range',torque,'The extensor redesign uses a cam to alter moment arm and actuator travel across knee angle. The two real configurations show the implemented route, not a generic ideal pulley. The design prediction must still be checked on the hardware. Source: Chapter 4 cam-route figure and geometry.')
image(s,pdf('30-torque/cam_pose_A.pdf'),.5,1.2,4.4,3.55);image(s,pdf('30-torque/cam_pose_H.pdf'),5.1,1.2,4.4,3.55);finish(s,'Joint mechanics')
s=start('All nine cam-route measurements exceed the selected human target',torque,'Nine measured extensor configurations meet or exceed the selected vastus-medialis reference. This is evidence of isometric torque capacity at the tested angles and pressures. It does not establish dynamic control or walking. The exact actuator sizes, pressure range, dimensions, and the pressure-specific prediction discrepancy at -86.5 degrees remain in notes and the detailed questions. Source: final extensor physical tests.',old=16)
for sh in s['tree'].findall('p:sp',N):
 if strtext(sh).startswith('Two 20 mm') or strtext(sh).startswith('Pressure '):settext(sh,'Physical knee test')
finish(s,'Joint mechanics')
s=start('Measured force and transmission mechanics make torque targets testable',torque+'; '+pub,'The mechanical contribution is a chain of experimentally checked predictions: actuator force, transmission through the joint, and physical torque against a selected human-muscle reference. The remaining problem is to use that mechanical model in a neural-control experiment. Sources: dissertation Chapters 3 and 4.')
for x,t,img in [(.5,'Actuator force',media(6,4)),(3.65,'Joint prediction',media(11,4)),(6.8,'Physical validation',media(15,4))]:image(s,img,x,1.4,2.7,2.9);addtext(s['tree'],t,x,4.5,2.7,.55,19,bold=True,align='ctr')
finish(s,'Joint mechanics')

# 49-64: a compact tools narrative in Table 6.1 order.
divider('Neuromechanical Models and Tools','neural','Neural models and tools','The mechanical model supplies a testable body. Adding a neural controller creates a coupled experiment, whose behavior depends on the controller, sensory signals, mechanics, support, and implementation. Source: Chapter 6.')
s=start('A neuromechanical model closes the loop through the body',tools,'Motor activity produces actuator force, which moves and loads the mechanical model. Length, force, and contact signals return to the controller. This is the distinction between a biomechanical model and a neuromechanical model, and it explains why transplanting a circuit can change its behavior. Source: dissertation Chapter 6.')
image(s,media(29,5),.5,1.1,2.2,3.8);loop(s,3.7,1.7);finish(s,'Neural models and tools')
old(18,'We combined Deng’s neural controller with Li’s separate mechanical model','Neural models and tools','The first walker attempt combines the neural controller from one neuromechanical model with the mechanics of another. It is not a direct replay of either original complete model. This distinction sets up why successful rhythmic activity did not establish a successful controller transplant. Sources: Li et al. and Deng et al. (2019); dissertation first W2L attempt.',tools+' · Li–Deng/W2L')
old(19,'Supported leg cycling does not establish self-supported walking','Neural models and tools','The paired native recordings contrast a fixed pelvis on a translating platform with release under the model’s weight. This is a new support-condition comparison, distinct from the historical suspended-in-air recording. It shows why rhythmic leg activity alone is insufficient. Source: Ben’s frozen-root treadmill.mp4 and free_drop.mp4, supplied October 31, 2026; original W2L attempt.',tools+' · Li–Deng/W2L')
old(22,'Controller transfer required a faster and maintainable simulation workflow','Neural models and tools','The failed transfer motivates the tools. The next task was to add and inspect controllers on a more complicated model such as Gait2392 while exchanging neural commands and mechanical sensory signals. AnimatLab’s long development gap was a practical limitation. The tools that follow support specific comparisons, rather than establishing unsupported walking.',tools)
old(21,'Sensory feedback changes bilateral foot-loading coordination','Neural models and tools','With the saved controller, drive, and support fixed, removing all feedback shifts contact phase from 0.498 to 0.223 cycles. Cycling and forward progression persist, but the feedback-disabled case fails the stated supported-stepping phase criterion. This establishes a contribution to coordination, not proof that feedback alone creates locomotion. Source: Chapter 6 W2L feedback comparison.',tools+' · 1. W2L controller')
old(29,'Gait2392 supplies paired human and BPA-driven robot models','Neural models and tools','A human reference and BPA-driven counterpart allow mechanics and actuation choices to be tested on related coordinates. Conversion and routing repairs require their own checks. These are model implementations, not a claim of complete equivalence between physics engines or of harness-free gait. Source: Chapter 6 Gait2392-based models.',tools+' · 2. Gait2392-based models')
old(32,'A persistent rhythm does not guarantee the correct knee motion','Neural models and tools','The joint-angle comparison shows why rhythm alone is not the acceptance criterion. The human-model simulation retains kinematic differences relative to the OpenSim reference even while cycles persist. Detailed hip, contact, and force comparisons are available in the questions. Source: Chapter 6 Gait2392 verification.',tools+' · 2. Gait2392 verification')
s=start('The graphical editor exposes controller and actuator definitions',tools+' · 3. SNS-Simscape library and editor','The library supplies sensory, neuron, synapse, activation, and actuator components, and the editor assembles them into a circuit. The displayed interface is the actual editor. The labeled diagram describes the signal exchange without treating decorative icons as evidence. The BPA force blocks use the characterized force model. Source: Chapter 6 SNS-Simscape library and editor.')
image(s,media(40,4),.5,1.2,5.6,3.5);box(s['tree'],'Sensory signal',6.5,1.45,2.7,.65,18);box(s['tree'],'Neural circuit',6.5,2.65,2.7,.65,18);box(s['tree'],'Actuator force',6.5,3.85,2.7,.65,18);line(s['tree'],7.85,2.2,7.85,2.5,arrow=True);line(s['tree'],7.85,3.4,7.85,3.7,arrow=True);finish(s,'Neural models and tools')
old(42,'Sensory feedback reduces sag in the loaded-arm verification','Neural models and tools','The loaded-arm case checks the full sensory-to-motor pathway. The reported sag is 2.00 degrees with feedback and 7.85 degrees without it in this specified test. This is a corrective reflex demonstration, not a general balance or locomotion result. Source: Chapter 6 loaded-arm verification.',tools+' · 3. Coupled reflex verification')
old(44,'The CAD workflow preserves assembly geometry, joints, and inertias','Neural models and tools','The CAD-derived route makes a separate mechanical implementation available. Assembly geometry, joints, frames, and inertias must be carried into simulation before neural and BPA interfaces are added. Component checks precede the coupled experiments. Source: Chapter 6 SOLIDWORKS-derived assembly.',tools+' · 4. CAD-derived robot')
old(47,'An abductor intervention prevents leg crossing in the supported CAD test','Neural models and tools','This is a controlled intervention in the supported AH001 test: the added abductor command changes the separation diagnostic and prevents crossing in the stated window. The support and hybrid control conditions are essential to the interpretation. It is not an unsupported walking result. Source: CAD-model integration verification.',tools+' · 4. CAD-model verification')
old(51,'The MuJoCo–Simulink bridge preserves the tested sensory signals','Neural models and tools','The bridge exchanges actuator commands and measured sensory signals between a graphical controller and MuJoCo. The overlay reports agreement for the specified five channels and tested command protocol. It establishes this interface check, not complete agreement of every circuit or mechanics implementation. Source: Chapter 6 MuJoCo–Simulink verification.',tools+' · 5. MuJoCo–Simulink integration')
old(53,'The experiment tools make Hunt’s four training steps observable','Neural models and tools','Hunt’s sequence separates intrinsic rhythm, feedback with the CPG, pattern formation with motor output, and the final sensory-motor refinements. The curriculum tools provide measured inputs and outputs for those stages rather than tuning the full network as one opaque object. Source: A. J. Hunt, Real-time control of dynamic walking in a dog-like robot, INRC2 Meeting, Intel Labs, February 15, 2019.','A. J. Hunt, INRC2 Meeting, Intel Labs (Feb 15, 2019) · Dissertation Ch. 6')
old(62,'Targeted stimuli distinguish rhythm changes from motor-output changes','Neural models and tools','An intervention at the rhythm-generator level can remove bursts and shift later phase. Downstream interventions can change motor output without the same phase effect, as detailed in the question comparisons. The tools record both neural and mechanical responses under fixed conditions. Source: Chapter 6 experiment drivers and Appendix D perturbation protocols.',tools+' · 6. Controlled interventions')
old(64,'Neural encoding is another test variable in a coupled model','Neural models and tools','The graded-versus-spiking comparison reports different timing in the specified implementation and input protocol. Literature circuit motifs guide such mechanistic tests; translated circuit rules are not presented as fully verified running reproductions. The result shows why a circuit name alone does not establish equivalent behavior. Source: Chapter 6 literature circuits and spiking variants.',tools+' · 7. Literature circuits and spiking variants')

# 65-70: overview, enabled future experiments, and conclusion.
s=start('The contributions make linked design choices experimentally testable','Dissertation Overview','The actuator model establishes available force, the joint corrections establish delivered torque, the sensory resource retains evidence with its context, and the simulation tools support controlled neural-mechanical tests. These are foundations for increasingly biomimetic robots and experiments, not a claim that the complete autonomous humanoid is finished. Sources: dissertation Overview and Conclusion.');unified(s);finish(s,'Overview and future work')
old(67,'Dynamic loading is the next test of the characterized joint routes','Overview and future work','The torque measurements establish isometric capacity at tested conditions. Dynamic loads, pressure-flow dynamics, and motion control are the next experiments. This separates what the present mechanical validation establishes from what remains to be demonstrated. Source: dissertation future work.','Dissertation future work · Dynamic actuation')
old(68,'Measured support reduction can test what sustains loaded movement','Overview and future work','Future work should reduce assistance in measured stages while monitoring posture, contact, and joint motion. Rhythmic stepping under a harness cannot substitute for this test. The original and newer controller comparisons can help identify which sensory and coordination mechanisms matter under load. Source: dissertation future work.','Dissertation future work · Support and balance')
old(69,'Robot interventions can test biological hypotheses','Overview and future work','With mechanics and neural pathways separately accessible, experiments can test how a change in routing, compliance, or feedback alters function. The results can generate biological hypotheses and inform future assistive mechanisms. The biological or clinical transfer remains a further validation step. Source: dissertation future work.','Dissertation future work · Biological and assistive experiments')
s=start('Other researchers can vary controllers and mechanics in the same workflow','Dissertation Overview and Future Work','The enabled use is a controlled model comparison: choose a mechanical model, use an inspectable neural controller, specify sensory and actuator interfaces, then record a response to a stated intervention. The database can support pathway selection with its experimental context; its detailed demonstration is at the end of the questions. Source: dissertation Overview and model-construction workflow.')
for x,t in [(.5,'Mechanical model'),(3.65,'Neural controller'),(6.8,'Specified intervention')]:box(s['tree'],t,x,1.6,2.7,.85,21)
for x in [1.85,5.0,8.15]:line(s['tree'],x,2.65,5.0,3.25,arrow=True)
box(s['tree'],'Neural + mechanical measurements',1.5,3.55,7.0,.9,25,True);finish(s,'Overview and future work')
old(70,'Measured actuation and controlled models make biomimetic design choices testable','Conclusion','The dissertation contributes an experimentally grounded BPA force model, joint design factors and physically checked torque targets, a sensory-evidence resource, and tools for controlled neuromechanical experiments. These advance the design and testing of biomimetic robots. Dynamic and unsupported behavior remain future validation tasks. Thank you. Sources: dissertation Conclusion and the evidence presented in the talk.','Dissertation Conclusion')

assert len(manifest)==70,(len(manifest),[(i+1,x['title']) for i,x in enumerate(manifest)])

# Remove the criticized decorative block artwork from retained technical questions.
for no in [17,30,35,37,38,39,45]:
 s=start(next((strtext(x) for x in read(f'ppt/slides/slide{no}.xml').findall('p:cSld/p:spTree/p:sp',N) if len(strtext(x))>35),'Neural and actuator interfaces'),tools,old=no)
 if no in [17,30,37,38,39,45]:
  if no==45:
   for sh in list(s['tree'].findall('p:sp',N)):
    if strtext(sh).startswith('Length and vel'):s['tree'].remove(sh)
  for pic in list(s['tree'].findall('p:pic',N)):
   bl=pic.find('.//a:blip',N)
   if bl is None:continue
   rid=bl.get(f'{{{R}}}embed');target=next((q.get('Target') for q in s['rels'] if q.get('Id')==rid),'')
   labels={'image50.png':'BPA force models\n10 mm  |  20 mm','image36.png':'Neural\ncontroller','image37.png':'Length / velocity\nfeedback','image54.png':'Synapse'}
   name=target.rsplit('/',1)[-1]
   if name in labels:
    tr=pic.find('p:spPr/a:xfrm',N);of=tr.find('a:off',N);ex=tr.find('a:ext',N);s['tree'].remove(pic);box(s['tree'],labels[name],int(of.get('x'))/EMU,int(of.get('y'))/EMU,int(ex.get('cx'))/EMU,max(.55,int(ex.get('cy'))/EMU),14)
 else:
  clear(s)
  for i,t in enumerate(['Sensory state','Neuron / synapse','Motor activation','Actuator force']):
   x=.5+i*2.35;box(s['tree'],t,x,2.0,2.1,.85,18)
   if i<3:line(s['tree'],x+2.12,2.42,x+2.3,2.42,arrow=True)
  addtext(s['tree'],'Signals and equations remain inspectable in the model',.7,3.8,8.5,.75,23,bold=True,align='ctr')
 put(f'ppt/slides/slide{no}.xml',s['root']);put(f'ppt/slides/_rels/slide{no}.xml.rels',s['rels'])

# Ben's later recording of the newer saved controller, with Root fixed.
s=start('The newer bilateral-RG controller also moves with the root fixed','Native AnimatLab recording supplied by Ben; newer saved W2L controller',old=20)
clear(s)
image(s,media(20,4),.45,1.3,5.5,3.1)
image(s,Path('C:/Users/Ben Bolen/Videos/Screen Recordings/w2l_bilateralRGfreeze.gif'),6.25,1.25,2.55,3.3)
addtext(s['tree'],'Original shared-RG circuit',.6,4.65,5.0,.35,17,bold=True)
addtext(s['tree'],'Newer controller, Root fixed',5.9,4.65,3.7,.35,17,bold=True)
put('ppt/slides/slide20.xml',s['root']);put('ppt/slides/_rels/slide20.xml.rels',s['rels'])
note(20,'The diagram at left is the original shared-RG circuit, included for architectural context. The animation at right is Ben’s w2l_bilateralRGfreeze.gif: the newer saved model with type II and contact feedback and separate left/right RG layers, with its Root fixed. This does not demonstrate self-supported walking. Recording predates the subsequently requested body/joint alignment review; preserve it as recorded, without substituting a corrected-model claim.')

# Preserve every existing slide and put all dedicated SADb slides last.
chosen=[x['source'] for x in manifest];sadb=[71,23,24,25,26,27,28];old_order=[x['source_slide'] for x in json.loads((BASE.parents[1]/'build/slide_order.json').read_text())]
questions=[x for x in old_order if x not in chosen and x not in sadb]+sadb;order=chosen+questions
assert len(order)==len(set(order)) and set(range(1,87)).issubset(order)
oldids=list(pres.find('p:sldIdLst',N));oldrelmap={x.get('Id'):x.get('Target') for x in pr};idmap={int(oldrelmap[x.get(f'{{{R}}}id')].split('slide')[-1].split('.')[0]):x for x in oldids}
lst=pres.find('p:sldIdLst',N)
for x in list(lst):lst.remove(x)
for no in newids:idmap[no]=E.Element(f'{{{P}}}sldId',id=str(5000+no),**{f'{{{R}}}id':'rIdNew'+str(no)})
mapping=[]
for pos,no in enumerate(order,1):
 lst.append(idmap[no]);r=read(f'ppt/slides/slide{no}.xml');r.set('show','1' if pos<=70 else '0')
 for sh in r.findall('p:cSld/p:spTree/p:sp',N):
  of=sh.find('.//a:xfrm/a:off',N);tt=sh.findall('.//a:t',N)
  if of is not None and int(of.get('y'))>4900000 and len(tt)==1 and (tt[0].text or '').isdigit():tt[0].text=str(pos)
 put(f'ppt/slides/slide{no}.xml',r);mapping.append(dict(position=pos,source_slide=no,hidden=pos>70,sadb=no in sadb))
put('ppt/presentation.xml',pres);put('ppt/_rels/presentation.xml.rels',pr);put('[Content_Types].xml',ct)
app=read('docProps/app.xml')
for x in app:
 if x.tag.endswith('}Slides'):x.text=str(len(order))
 elif x.tag.endswith('}HiddenSlides'):x.text=str(len(questions))
put('docProps/app.xml',app)
candidate=BUILD/'candidate.pptx'
with zipfile.ZipFile(candidate,'w',zipfile.ZIP_DEFLATED) as z:
 for n,b in parts.items():z.writestr(n,b)
(BUILD/'main_story.json').write_text(json.dumps(manifest,indent=2));(BUILD/'slide_order.json').write_text(json.dumps(mapping,indent=2))
owners=[i+1 for i,no in enumerate(order) if read(f'ppt/slides/slide{no}.xml').find('.//c:chart',N) is not None]
(BUILD/'native_chart_owners.json').write_text(json.dumps(owners));(BUILD/'build_summary.json').write_text(json.dumps(dict(slides=len(order),main=70,questions=len(questions),seconds=sum(s['seconds'] for s in manifest),new_slides=len(newids),native_chart_owners=owners,baseline_sha256=hashlib.sha256(BASE.read_bytes()).hexdigest()),indent=2))
import shutil
shutil.copyfile(BASE.parent/'sadb_app.html',OUT/'sadb_app.html')
print(json.dumps(dict(slides=len(order),main=70,questions=len(questions),new=len(newids),minutes=sum(x['seconds'] for x in manifest)/60,native_chart_owners=owners)))
