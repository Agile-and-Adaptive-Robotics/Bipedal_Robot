"""Run after apply_revision.ps1. Preserve existing slides and native chart parts."""
from pathlib import Path
from lxml import etree as E
from PIL import Image
import zipfile,copy,json
W=Path(__file__).resolve().parent;ROOT=Path.cwd();C=ROOT/'Neuromechanical_Models/W2L_Comparison_20261030'
P='http://schemas.openxmlformats.org/presentationml/2006/main';A='http://schemas.openxmlformats.org/drawingml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships';PK='http://schemas.openxmlformats.org/package/2006/relationships';CT='http://schemas.openxmlformats.org/package/2006/content-types';N={'p':P,'a':A,'r':R}
path=W/'build/candidate.pptx'
with zipfile.ZipFile(path) as z:d={n:z.read(n) for n in z.namelist()}
def read(n):return E.fromstring(d[n])
def write(n,r):d[n]=E.tostring(r,xml_declaration=True,encoding='UTF-8',standalone=True)
def text(s,t):
 b=s.find('p:txBody',N);para=b.find('a:p',N);pp=copy.deepcopy(para.find('a:pPr',N));rp=copy.deepcopy(para.find('a:r/a:rPr',N))
 for p in list(b):
  if p.tag==f'{{{A}}}p':b.remove(p)
 p=E.SubElement(b,f'{{{A}}}p')
 if pp is not None:p.append(pp)
 r=E.SubElement(p,f'{{{A}}}r')
 if rp is not None:r.append(rp)
 E.SubElement(r,f'{{{A}}}t').text=t
def shape_id(s):
 n=s.find('.//p:cNvPr',N);return n.get('id') if n is not None else None
def picture(rid,name,x,y,w,h,crop=None):
 p=E.Element(f'{{{P}}}pic');nv=E.SubElement(p,f'{{{P}}}nvPicPr');E.SubElement(nv,f'{{{P}}}cNvPr',id='50',name=name,descr=name);cp=E.SubElement(nv,f'{{{P}}}cNvPicPr');E.SubElement(cp,f'{{{A}}}picLocks',noChangeAspect='1');E.SubElement(nv,f'{{{P}}}nvPr')
 f=E.SubElement(p,f'{{{P}}}blipFill');E.SubElement(f,f'{{{A}}}blip',{f'{{{R}}}embed':rid})
 if crop:E.SubElement(f,f'{{{A}}}srcRect',**{k:str(v) for k,v in crop.items()})
 st=E.SubElement(f,f'{{{A}}}stretch');E.SubElement(st,f'{{{A}}}fillRect');sp=E.SubElement(p,f'{{{P}}}spPr');tr=E.SubElement(sp,f'{{{A}}}xfrm');E.SubElement(tr,f'{{{A}}}off',x=str(x),y=str(y));E.SubElement(tr,f'{{{A}}}ext',cx=str(w),cy=str(h));g=E.SubElement(sp,f'{{{A}}}prstGeom',prst='rect');E.SubElement(g,f'{{{A}}}avLst');return p
def rel(relroot,rid,kind,target):E.SubElement(relroot,f'{{{PK}}}Relationship',Id=rid,Type=R+'/'+kind,Target=target)
def add_notes(no,prose):
 n=read('ppt/notesSlides/notesSlide70.xml')
 for s in n.findall('.//p:sp',N):
  ph=s.find('p:nvSpPr/p:nvPr/p:ph',N)
  if ph is not None and ph.get('type')=='body':text(s,prose)
 write(f'ppt/notesSlides/notesSlide{no}.xml',n)
 nr=read('ppt/notesSlides/_rels/notesSlide70.xml.rels')
 for x in nr:
  if x.get('Type').endswith('/slide'):x.set('Target',f'../slides/slide{no}.xml')
 write(f'ppt/notesSlides/_rels/notesSlide{no}.xml.rels',nr)
content=read('[Content_Types].xml')
slide=read('ppt/slides/slide19.xml');tree=slide.find('p:cSld/p:spTree',N)
for s in list(tree):
 if shape_id(s)=='14':tree.remove(s)
 elif shape_id(s)=='15':text(s,'Native AnimatLab capture at 2.14 s')
 elif shape_id(s)=='16':text(s,'Released onto the platform')
tree.append(picture('rIdNativeDrop','Actual W2L free-drop in AnimatLab at 2.141 s',4750594,1830000,4071938,2072326,dict(l=40625,r=30208,t=43269,b=29327)))
d['ppt/media/w2l_native_drop.png']=(C/'screenshots/legacy_drop_clean.png').read_bytes();rels=read('ppt/slides/_rels/slide19.xml.rels');rel(rels,'rIdNativeDrop','image','../media/w2l_native_drop.png');write('ppt/slides/_rels/slide19.xml.rels',rels);write('ppt/slides/slide19.xml',slide)
notes=read('ppt/notesSlides/notesSlide19.xml')
for s in notes.findall('.//p:sp',N):
 ph=s.find('p:nvSpPr/p:nvPr/p:ph',N)
 if ph is not None and ph.get('type')=='body':text(s,'Left: Ben\'s original 2021 native AnimatLab recording. Right: new native AnimatLab screenshot at t=2.141 s using the saved legacy controller and original initial pose, with Root free and a static platform. No controller retuning or extra support. The new matched campaign records MN voltages, joint angles, contact counts, and Root position. This capture is a separate test, not an archived released segment of the 2021 video. The recording project is Neuromechanical_Models/W2L_Comparison_20261030/recording_projects/legacy_drop/legacy_drop.aproj. A moving released video can replace this still after recording.')
write('ppt/notesSlides/notesSlide19.xml',notes)
items=[
 ('Half-weight support delays the fall but does not stabilize either model','defense_support','Native AnimatLab | matched mechanics | unchanged saved controllers','Eight baselines, two saved controllers by four support/contact conditions. Partial support is a constant 205.4031 N positive-Y force at Root, 50% of declared 41.919 kg weight at gravity 9.8 m/s^2. Root is free in both fall conditions. There is no balance or rotational constraint. Both controllers fall. Root first passes below 0.5 m at 1.103 versus 1.790 s (legacy) and 0.999 versus 2.077 s (newer), without versus with support. The fixed-root air and belt trials retain Root Y=0.99298 m. The newer partly supported model later moves below belt height. No unsupported walking is demonstrated.'),
 ('PF and RG pulses produce different motor responses in the saved controllers','defense_deng2019','Deng et al. (2019), Biomimetics 4(1):21 | pulse protocol transferred to W2L','Protocol transfer from Deng 2019 figures 8-10: -10 nA for 100 ms at t=2.0-2.1 s, applied separately to L Hip PF ext and L RG ext. Plot shows L Hip MN ext in air for both saved controllers, with the unperturbed trace in gray. Both Root and floors are configured so all contact counts are zero. All four mechanical conditions were also tested and are available in the campaign folder. This is a transfer of the pulse protocol to W2L, not an exact replication of the rat model. Source https://www.mdpi.com/2313-7673/4/1/21'),
 ('The saved controllers show different phase responses with mechanical contact','defense_deng2022','Deng et al. (2022), Biomimetics 7(4):226 | pulse protocol transferred to W2L','Protocol transfer: +/-0.5 nA to L RG flx IN, duration 5% of the measured baseline period, ten phases 0.05 to 0.95, both saved controllers in air and on the moving belt. Phase reference: L RG ext adaptive-midrange upward crossing. Positive shift is an advance of the first post-stimulus cycle relative to baseline, divided by the baseline period. Baseline periods: legacy air/belt 0.443 s; newer air 0.493 s; newer belt 0.4645 s. The newer controller has four Type II adapters, four contact adapters, and separate L/R RGs. Gains are unchanged. The comparison does not isolate individual feedback contributions. We did not replicate the original rat model, 1/2-Hz retuning, or synaptic-strength sweep. Source https://www.mdpi.com/2313-7673/7/4/226')]
pres=read('ppt/presentation.xml');lst=pres.find('p:sldIdLst',N);pr=read('ppt/_rels/presentation.xml.rels');mapping=json.loads((W/'build/slide_order.json').read_text());newmaps=[]
for k,(title,fig,credit,prose) in enumerate(items):
 no=84+k;r=read('ppt/slides/slide70.xml');r.set('show','0');t=r.find('p:cSld/p:spTree',N)
 for s in list(t):
  sid=shape_id(s)
  if s.tag in [f'{{{P}}}nvGrpSpPr',f'{{{P}}}grpSpPr']:continue
  if sid not in ['7','8','18']:t.remove(s)
  elif sid=='7':text(s,title)
  elif sid=='8':text(s,credit)
  else:text(s,str(77+k))
 image=C/'figures'/(fig+'.png');w,h=Image.open(image).size;scale=min(8429625/w,3400000/h);dw,dh=int(w*scale),int(h*scale);t.append(picture('rIdEvidence',title,321470+(8429625-dw)//2,1050000+(3400000-dh)//2,dw,dh));d['ppt/media/'+fig+'.png']=image.read_bytes()
 write(f'ppt/slides/slide{no}.xml',r);rr=E.Element(f'{{{PK}}}Relationships',nsmap={None:PK});oldrels=read('ppt/slides/_rels/slide70.xml.rels')
 for x in oldrels:
  if x.get('Type').endswith('/slideLayout'):rr.append(copy.deepcopy(x))
 rel(rr,'rIdEvidence','image','../media/'+fig+'.png');rel(rr,'rIdNotes','notesSlide',f'../notesSlides/notesSlide{no}.xml');write(f'ppt/slides/_rels/slide{no}.xml.rels',rr);add_notes(no,prose)
 rid='rIdW2LEvidence'+str(no);rel(pr,rid,'slide',f'slides/slide{no}.xml');node=E.Element(f'{{{P}}}sldId',id=str(4000+no),**{f'{{{R}}}id':rid});lst.insert(76+k,node)
 for part,typ in [('slides','slide'),('notesSlides','notesSlide')]:E.SubElement(content,f'{{{CT}}}Override',PartName=f'/ppt/{part}/{typ}{no}.xml',ContentType=f'application/vnd.openxmlformats-officedocument.presentationml.{typ}+xml')
 newmaps.append(dict(position=77+k,source_slide=no,hidden=True,sadb=False,new_evidence=True))
mapping=mapping[:76]+newmaps+mapping[76:]
for pos,row in enumerate(mapping,1):
 row['position']=pos
 r=read(f"ppt/slides/slide{row['source_slide']}.xml")
 for s in r.findall('p:cSld/p:spTree/p:sp',N):
  off=s.find('.//a:xfrm/a:off',N);tt=s.findall('.//a:t',N)
  if off is not None and int(off.get('y'))>4900000 and len(tt)==1 and tt[0].text.isdigit():tt[0].text=str(pos)
 write(f"ppt/slides/slide{row['source_slide']}.xml",r)
write('ppt/presentation.xml',pres);write('ppt/_rels/presentation.xml.rels',pr);write('[Content_Types].xml',content)
app=read('docProps/app.xml');nsapp={'e':'http://schemas.openxmlformats.org/officeDocument/2006/extended-properties'}
for tag,v in [('Slides','86'),('HiddenSlides','22')]:
 el=app.find('e:'+tag,nsapp)
 if el is not None:el.text=v
write('docProps/app.xml',app)
tmp=path.with_suffix('.tmp')
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
 for n,b in d.items():z.writestr(n,b)
tmp.replace(path);(W/'build/slide_order.json').write_text(json.dumps(mapping,indent=2));print('86 slides; native drop capture on 19; experiment questions 77-79; SADb 80-86')
