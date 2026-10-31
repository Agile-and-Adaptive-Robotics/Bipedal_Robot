from pathlib import Path
from lxml import etree as E
import zipfile,copy,json,hashlib
W=Path(__file__).resolve().parent;p=W/'build/candidate.pptx';BASE=W.parent/'C3NS_Revision_20261030_v16/output/Bolen_PhD_Defense_RoseGreen_20261030_v16.pptx'
P='http://schemas.openxmlformats.org/presentationml/2006/main';A='http://schemas.openxmlformats.org/drawingml/2006/main';R='http://schemas.openxmlformats.org/officeDocument/2006/relationships';PK='http://schemas.openxmlformats.org/package/2006/relationships';N={'p':P,'a':A};EMU=914400
with zipfile.ZipFile(p) as z:d={n:z.read(n) for n in z.namelist()}
with zipfile.ZipFile(BASE) as z:old=E.fromstring(z.read('ppt/slides/slide19.xml'));oldrels=E.fromstring(z.read('ppt/slides/_rels/slide19.xml.rels'))
r=E.fromstring(d['ppt/slides/slide19.xml']);tree=r.find('p:cSld/p:spTree',N);rels=E.fromstring(d['ppt/slides/_rels/slide19.xml.rels']);template=next(x for x in old.findall('p:cSld/p:spTree/p:pic',N) if x.find('.//a:videoFile',N) is not None)
for x in list(tree):
 if x.tag==f'{{{P}}}pic':tree.remove(x)
 elif x.tag==f'{{{P}}}sp':
  sid=x.find('p:nvSpPr/p:cNvPr',N).get('id')
  if sid in ['6','16']:tree.remove(x);continue
  txt=x.findall('.//a:t',N);whole=' '.join(t.text or '' for t in txt)
  replacements={'Original recording · 2021':'Root fixed, moving platform','Native AnimatLab capture at 2.14 s':'Root released'}
  if whole in replacements:
   txt[0].text=replacements[whole]
   for t in txt[1:]:t.text=''
tim=copy.deepcopy(old.find('p:timing',N));rootchildren=tim.find('p:tnLst/p:par/p:cTn/p:childTnLst',N);mainseq=rootchildren.find('p:seq/p:cTn/p:childTnLst',N);auto=copy.deepcopy(mainseq[0]);vid=copy.deepcopy(rootchildren.find('p:video',N));inter=copy.deepcopy(rootchildren.findall('p:seq',N)[1])
for branch in [auto,vid,inter]:
 for t in branch.findall('.//p:cTn',N):t.set('id',str(int(t.get('id'))+100))
 for t in branch.findall('.//p:spTgt',N):t.set('spid','14')
# Both play commands share the first automatic parallel group. A second main
# sequence entry waits for the first indefinitely looping movie to complete.
first_group=mainseq[0].find('.//p:cTn[@id="4"]/p:childTnLst',N)
second_effect=auto.find('.//p:cTn[@id="104"]/p:childTnLst/p:par',N)
first_group.append(copy.deepcopy(second_effect));rootchildren.append(vid);rootchildren.append(inter)
for oldtim in r.findall('p:timing',N):r.remove(oldtim)
r.append(tim)
records=[]
for sid,key,x,w in [(4,'treadmill',.65,4.1),(14,'free_drop',5.25,4.1)]:
 pic=copy.deepcopy(template);pic.find('p:nvPicPr/p:cNvPr',N).attrib.update({'id':str(sid),'name':key+' native recording','descr':'Native AnimatLab recording: '+('fixed Root and moving platform' if key=='treadmill' else 'Root released and model falling')});mp4='v17_'+key+'.mp4';poster='v17_'+key+'_poster.png';d['ppt/media/'+mp4]=(W/'assets'/(key+'_recording.mp4')).read_bytes();d['ppt/media/'+poster]=(W/'assets'/(key+'_poster.png')).read_bytes()
 mapping={}
 for name in ['videoFile','media','blip']:
  matches=pic.xpath('.//*[local-name()=$nm]',nm=name)
  for node in matches:
   for attr,val in list(node.attrib.items()):
    if attr.startswith('{'+R+'}'):
     if val not in mapping:
      oldrel=next(q for q in oldrels if q.get('Id')==val);rid='rId17'+key+str(len(mapping));mapping[val]=rid;target='../media/'+(poster if name=='blip' else mp4);E.SubElement(rels,f'{{{PK}}}Relationship',Id=rid,Type=oldrel.get('Type'),Target=target)
     node.set(attr,mapping[val])
 tr=pic.find('p:spPr/a:xfrm',N);of=tr.find('a:off',N);ex=tr.find('a:ext',N);height=3.1;width=(500/480 if key=='treadmill' else 418/480)*height
 of.attrib.update({'x':str(int((x+(w-width)/2)*EMU)),'y':str(int(1.25*EMU))});ex.attrib.update({'cx':str(int(width*EMU)),'cy':str(int(height*EMU))})
 crop=pic.find('p:blipFill/a:srcRect',N)
 if crop is not None:pic.find('p:blipFill',N).remove(crop)
 tree.append(pic);records.append(dict(shape_id=sid,source=str(Path('C:/Users/Ben Bolen/Videos/Screen Recordings')/('frozen-root treadmill.mp4' if key=='treadmill' else 'free_drop.mp4')),embedded=mp4,sha256=hashlib.sha256(d['ppt/media/'+mp4]).hexdigest(),muted=True,looped=True,autoplay=True))
for name,xml in [('ppt/slides/slide19.xml',r),('ppt/slides/_rels/slide19.xml.rels',rels)]:d[name]=E.tostring(xml,xml_declaration=True,encoding='UTF-8',standalone=True)
tmp=p.with_suffix('.tmp')
with zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as z:
 for n,b in d.items():z.writestr(n,b)
tmp.replace(p);(W/'build/native_recordings.json').write_text(json.dumps(records,indent=2));print('Embedded two muted, looping native recordings on main slide 52')
