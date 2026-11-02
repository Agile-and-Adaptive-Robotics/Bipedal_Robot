from pathlib import Path
import zipfile,json,hashlib,posixpath,re
from lxml import etree as E
from pypdf import PdfReader
W=Path(__file__).resolve().parent
B=W.parent/'C3NS_Revision_20261031_v18/output/Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx'
N={'p':'http://schemas.openxmlformats.org/presentationml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
with zipfile.ZipFile(B) as z:
 rels={r.get('Id'):r.get('Target') for r in E.fromstring(z.read('ppt/_rels/presentation.xml.rels'))}
 slides=[]
 for i,sid in enumerate(E.fromstring(z.read('ppt/presentation.xml')).findall('p:sldIdLst/p:sldId',N),1):
  part=posixpath.normpath('ppt/'+rels[sid.get('{'+N['r']+'}id')]);root=E.fromstring(z.read(part))
  sr={r.get('Id'):r.get('Target') for r in E.fromstring(z.read(posixpath.dirname(part)+'/_rels/'+posixpath.basename(part)+'.rels'))}
  texts=[' '.join(x.itertext()) for x in root.findall('.//a:t',N)]
  ims=[]
  for j,pic in enumerate(root.findall('.//p:pic',N)):
   blip=pic.find('.//a:blip',N)
   if blip is None:continue
   rid=blip.get('{'+N['r']+'}embed')
   if rid not in sr:continue
   ip=posixpath.normpath(posixpath.dirname(part)+'/'+sr[rid]);name=f's{i:03d}_img{j}'+Path(ip).suffix
   (W/'assets'/name).write_bytes(z.read(ip));ims.append(name)
  slides.append(dict(slide=i,text=texts,images=ims,part=part))
 (W/'build/slides.json').write_text(json.dumps(slides,indent=2),encoding='utf8')
 for x in slides:print(x['slide'], ' / '.join(x['text'])[:200], x['images'])
pdf=Path('C:/Users/Ben Bolen/Documents/ME603 Dissertation/Defense/Bolen_Dissertation.pdf')
doc=PdfReader(pdf);matches=[]
for i,page in enumerate(doc.pages):
 t=page.extract_text();flat=re.sub(r'\s+',' ',t)
 for m in re.finditer('transmission',flat,re.I):matches.append(dict(pdf_page=i+1,context=flat[max(0,m.start()-160):m.end()+200]))
report=dict(pdf=str(pdf),sha256=hashlib.sha256(pdf.read_bytes()).hexdigest(),pages=len(doc.pages),matches=matches)
(W/'build/dissertation_wording_audit.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print('DISSERTATION AUDIT',json.dumps(report))
