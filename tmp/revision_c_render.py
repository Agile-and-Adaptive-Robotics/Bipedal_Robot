from pathlib import Path
import fitz,re,json
R=Path(r'D:\Github\Bipedal_Robot');C=R/'Documentation/Reports and Papers/Dissertation/Notes/advisor_revision_20261022/Revision_C';B=Path((C/'review/latest_build.txt').read_text())
d=fitz.open(B/'main.pdf');toc=(B/'main.toc').read_text(errors='replace')
print('PAGES',len(d));print(toc[-1800:].encode('ascii','backslashreplace').decode())
out=C/'review/rendered';out.mkdir(exist_ok=True)
for old in out.glob('*.png'):old.unlink()
groups={}
for n,p in enumerate(d):
    t=p.get_text()
    if 'The Sensory Afferent Database' in t and re.search(r'Chapter\s+5',t):groups['SADb']=n
    if 'Synthetic Nervous System Simulation' in t and re.search(r'Chapter\s+6',t):groups['SNS']=n
    if 'Appendix D: Supplementary Neuromechanical' in t:groups['Appendix_D']=n
    if 'Appendix E: Supporting SADb Records' in t:groups['Appendix_E']=n
    if 'Supporting Route Searches and Transmission Sweeps' in t and 'C.' in t and len(t)<3500:groups.setdefault('Appendix_C_routes',n)
groups={}
for level,title,page in d.get_toc():
    if level==1 and title.startswith('The Sensory'):groups['SADb']=page-1
    if level==1 and title.startswith('Synthetic Nervous'):groups['SNS']=page-1
    if level==1 and title.startswith('Appendix D:'):groups['Appendix_D']=page-1
    if level==1 and title.startswith('Appendix E:'):groups['Appendix_E']=page-1
    if 'Supporting Route Searches' in title:groups['Appendix_C_routes']=page-1
    if 'Cam-Route Pose Sequence' in title:groups['Cam_sequence']=page-1
print('GROUPS',groups)
# Render all new main-chapter and appendix pages, plus targeted mechanical/background figures.
selected=set()
for name,start in groups.items():
    end=min(len(d),start+({'SADb':5,'SNS':11,'Appendix_D':15,'Appendix_E':5,'Appendix_C_routes':10,'Cam_sequence':4}[name]))
    for n in range(start,end):selected.add(n)
for n,p in enumerate(d):
    t=p.get_text()
    if any(k in t for k in ["robot PAM routes.",'Cam-extensor route at','Cam-extensor poses','Initial cam-route holds and','Predicted insertion-torque magnitude versus','Torque margins at 620','Route geometry of the direct extensor']):selected.add(n)
for n in sorted(selected):d[n].get_pixmap(matrix=fitz.Matrix(1,1)).save(out/f'page_{n+1:03}.png')
# Contact sheets: page overview for visual layout review, original page images retained.
from PIL import Image,ImageOps,ImageDraw
pages=sorted(selected)
for batch in range(0,len(pages),9):
    sheet=Image.new('RGB',(936,1260),'#bbbbbb');draw=ImageDraw.Draw(sheet)
    for i,n in enumerate(pages[batch:batch+9]):
        im=Image.open(out/f'page_{n+1:03}.png');im.thumbnail((300,392));x=(i%3)*312+6;y=(i//3)*420+20;sheet.paste(im,(x,y));draw.text((x,y-15),f'PDF page {n+1}',fill='black')
    sheet.save(out/f'overview_{batch//9+1:02}.png')
(C/'review/render_manifest.json').write_text(json.dumps({'pdf_pages':len(d),'chapter_start_pdf_pages':{k:v+1 for k,v in groups.items()},'rendered_pdf_pages':[n+1 for n in pages]},indent=2))
print('RENDERED',len(pages),'pages')
