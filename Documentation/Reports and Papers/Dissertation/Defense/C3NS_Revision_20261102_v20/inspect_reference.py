from pathlib import Path
from lxml import etree as E
import zipfile,posixpath,json,re
W=Path(__file__).resolve().parent
P=W.parent/'Movement and Control of Biomimetic Humanoid Robots.pptx'
N={'p':'http://schemas.openxmlformats.org/presentationml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships'}
records=[]
with zipfile.ZipFile(P) as z:
 rels={r.get('Id'):r.get('Target') for r in E.fromstring(z.read('ppt/_rels/presentation.xml.rels'))}
 for i,s in enumerate(E.fromstring(z.read('ppt/presentation.xml')).findall('p:sldIdLst/p:sldId',N),1):
  part=posixpath.normpath('ppt/'+rels[s.get('{'+N['r']+'}id')]);root=E.fromstring(z.read(part))
  texts=root.xpath('//a:t/text()',namespaces=N);text=' '.join(texts)
  if not re.search('torque|moment|length|correction|offset|χ',text,re.I):continue
  sr={r.get('Id'):r.get('Target') for r in E.fromstring(z.read(posixpath.dirname(part)+'/_rels/'+posixpath.basename(part)+'.rels'))}
  pics=[]
  for j,pic in enumerate(root.findall('.//p:pic',N)):
   b=pic.find('.//a:blip',N)
   if b is None:continue
   rid=b.get('{'+N['r']+'}embed');ip=posixpath.normpath(posixpath.dirname(part)+'/'+sr[rid]);out=W/'assets'/f'old_talk_s{i:02d}_{j}{Path(ip).suffix}';out.write_bytes(z.read(ip));pics.append(out.name)
  records.append({'slide':i,'text':text,'images':pics});print(i,text[:800],pics)
(W/'build/reference_slides.json').write_text(json.dumps(records,indent=2,ensure_ascii=False),encoding='utf-8')
