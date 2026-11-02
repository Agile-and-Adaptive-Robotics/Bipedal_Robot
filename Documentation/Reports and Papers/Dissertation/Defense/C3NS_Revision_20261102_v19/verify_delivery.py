from pathlib import Path
import zipfile, json, hashlib, re, io
from lxml import etree as E
from PIL import Image
from pypdf import PdfReader
W=Path(__file__).resolve().parent
target=W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pptx'
baseline=W.parent/'C3NS_Revision_20261031_v18/output/Bolen_PhD_Defense_RoseGreen_20261031_v18.pptx'
N={'p':'http://schemas.openxmlformats.org/presentationml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
def sha(data):return hashlib.sha256(data).hexdigest()
with zipfile.ZipFile(target) as z,zipfile.ZipFile(baseline) as b:
 p=E.fromstring(z.read('ppt/presentation.xml'));rr=E.fromstring(z.read('ppt/_rels/presentation.xml.rels'));rm={x.get('Id'):x.get('Target') for x in rr}
 import posixpath
 parts=[posixpath.normpath('ppt/'+rm[x.get('{'+N['r']+'}id')]) for x in p.findall('p:sldIdLst/p:sldId',N)]
 parts=[x.replace('ppt/ppt/','ppt/') for x in parts]
 hidden=[];rejected=[];titles=[]
 for i,part in enumerate(parts,1):
  root=E.fromstring(z.read(part));texts=root.findall('.//a:t',N);words=' '.join(t.text or '' for t in texts)
  if root.get('show')=='0':hidden.append(i)
  if re.search('transmission',words,re.I):rejected.append({'slide':i,'text':words})
  titles.append({'slide':i,'text':words})
 for part in z.namelist():
  if part.startswith('ppt/notesSlides/') and part.endswith('.xml'):
   words=' '.join(E.fromstring(z.read(part)).xpath('//a:t/text()',namespaces=N))
   if re.search('transmission (model|path)',words,re.I):rejected.append({'part':part,'text':words})
 charts=[n for n in z.namelist() if re.fullmatch(r'ppt/charts/chart\d+.xml',n)]
 chart_equal={n:z.read(n)==b.read(n) for n in charts}
 books=[n for n in z.namelist() if n.startswith('ppt/embeddings/') and n.endswith('.xlsx')]
 book_equal={n:z.read(n)==b.read(n) for n in books}
 videos={sha(z.read(n)) for n in z.namelist() if n.lower().endswith('.mp4')}
 oldvideos={sha(b.read(n)) for n in b.namelist() if n.lower().endswith('.mp4')}
 gifs=[]
 for n in z.namelist():
  if n.lower().endswith('.gif'):
   im=Image.open(io.BytesIO(z.read(n)));gifs.append({'part':n,'frames':im.n_frames})
submitted=Path(r'C:\Users\Ben Bolen\Documents\ME603 Dissertation\Defense\Bolen_Dissertation.pdf')
report={'deck':str(target),'sha256':sha(target.read_bytes()),'bytes':target.stat().st_size,'slides':len(parts),'hidden_slides':hidden,'rejected_wording':rejected,'native_charts_unchanged':chart_equal,'embedded_workbooks_unchanged':book_equal,'native_video_hashes_preserved':videos==oldvideos,'videos':len(videos),'animated_gifs':gifs,'review_pdf_pages':len(PdfReader(W/'output/Bolen_PhD_Defense_RoseGreen_20261102_v19.pdf').pages),'submitted_pdf_unchanged':sha(submitted.read_bytes())=='447c2cc5a9dd219fa0857e298727fb68d4797af8f406251865eb22f4b126361f'}
(W/'build/delivery_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
(W/'build/final_slide_text.json').write_text(json.dumps(titles,indent=2,ensure_ascii=False),encoding='utf-8')
assert len(parts)==123 and hidden==list(range(71,124))
assert not rejected,rejected
assert all(chart_equal.values()) and all(book_equal.values()) and videos==oldvideos
assert report['submitted_pdf_unchanged'] and report['review_pdf_pages']==123
print(json.dumps({k:report[k] for k in ['sha256','bytes','slides','videos','native_video_hashes_preserved','submitted_pdf_unchanged','review_pdf_pages']}))
