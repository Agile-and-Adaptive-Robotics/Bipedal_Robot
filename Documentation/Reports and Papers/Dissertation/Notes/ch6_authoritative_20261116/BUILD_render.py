"""Render target chapter pages and extract actual compiled text sizes.

Requires a verified successful BUILD_RESULT_<revision>.json from the canonical
build. Images and contact sheets are inspection aids, not a replacement PDF.
"""
from pathlib import Path
import argparse, hashlib, json, math
import fitz
from PIL import Image, ImageDraw

NOTES=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--revision',default='r1')
args=parser.parse_args()
result=json.loads((NOTES/f'BUILD_RESULT_{args.revision}.json').read_text())
if result['status']!='compiled' or not result['sources_conserved_during_build']:
    raise SystemExit('A source-conserved successful canonical build is required')
pdf=Path(result['pdf'])
if hashlib.sha256(pdf.read_bytes()).hexdigest()!=result['pdf_sha256']:
    raise SystemExit('Compiled PDF hash does not match build record')
doc=fitz.open(pdf)
toc=doc.get_toc()
targets={}
for token,name in [('Tools for Constructing','chapter6'),('Appendix D:','appendixD')]:
    hits=[(n,item) for n,item in enumerate(toc) if token in item[1]]
    if len(hits)!=1: raise RuntimeError(f'Expected one bookmark for {token}: {hits}')
    n,item=hits[0]
    end=next((other[2]-1 for other in toc[n+1:] if other[0]<=item[0]),len(doc))
    targets[name]={'title':item[1],'pdf_pages_1based':[item[2],end]}
out=pdf.parent/'render_targets'
out.mkdir(exist_ok=True)
page_info=[]
for name,info in targets.items():
    thumbpages=[]
    for pno in range(info['pdf_pages_1based'][0]-1,info['pdf_pages_1based'][1]):
        page=doc[pno]
        path=out/f'{name}_pdfpage_{pno+1:03d}.png'
        page.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False).save(path)
        spans=[span for block in page.get_text('dict')['blocks'] if 'lines' in block
            for line in block['lines'] for span in line['spans'] if span['text'].strip()]
        small=[{'text':s['text'],'font':s['font'],'size_pt':s['size'],'bbox':list(s['bbox'])}
               for s in spans if s['size']<9.98]
        page_info.append({'group':name,'pdf_page_1based':pno+1,'printed_label':page.get_label(),
            'png':str(path),'minimum_text_span_pt':min(s['size'] for s in spans) if spans else None,
            'below_10pt_spans':small,'extracted_text':page.get_text()})
        thumbpages.append(path)
    for start in range(0,len(thumbpages),6):
        group=thumbpages[start:start+6]
        canvas=Image.new('RGB',(3*420,2*570),'#dddddd')
        draw=ImageDraw.Draw(canvas)
        for position,path in enumerate(group):
            image=Image.open(path)
            image.thumbnail((405,540))
            x=(position%3)*420+(420-image.width)//2
            y=(position//3)*570+25
            canvas.paste(image,(x,y))
            draw.text(((position%3)*420+12,(position//3)*570+6),path.stem,fill='black')
        canvas.save(out/f'{name}_contact_{start//6+1:02d}.png')
report={'status':'rendered_for_visual_review','pdf':str(pdf),'pdf_sha256':result['pdf_sha256'],
        'render_scale':2,'targets':targets,'toc':toc,'render_directory':str(out),
        'pages':page_info,'font_check_limit':'Actual compiled spans include valid mathematical superscripts/subscripts; review small spans in context. Raster figure text requires visual inspection.'}
(NOTES/f'BUILD_RENDER_{args.revision}.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k not in ('pages','toc')},indent=2))
