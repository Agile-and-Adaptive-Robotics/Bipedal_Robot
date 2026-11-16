"""Read-only comparison of already rendered r6/r7 review pages."""
from pathlib import Path
import json,hashlib
from PIL import Image,ImageChops
import pymupdf

HERE=Path(__file__).resolve().parent
ROOT=next(p for p in HERE.parents if p.name=='Bipedal_Robot')
TEMP=ROOT/'tmp/ch6_authoritative_20261116'
old=TEMP/'build_final_r6';new=TEMP/'build_final_r7'
a=pymupdf.open(old/'main.pdf');b=pymupdf.open(new/'main.pdf')
rows=[]
for oldpage in range(293,313):
    newpage=oldpage+2
    x=Image.open(old/'render_targets'/f'appendixD_pdfpage_{oldpage}.png').convert('RGB')
    y=Image.open(new/'render_targets'/f'appendixD_pdfpage_{newpage}.png').convert('RGB')
    box=(0,232,x.width,x.height)
    diff=ImageChops.difference(x.crop(box),y.crop(box))
    at=a[oldpage-1].get_text().split('\n',1)[1]
    bt=b[newpage-1].get_text().split('\n',1)[1]
    rows.append({'old_pdf_page':oldpage,'final_pdf_page':newpage,'body_image_excluding_page_number_exact':diff.getbbox()is None,
                 'difference_bbox_in_body_crop':diff.getbbox(),'text_excluding_first_page_number_exact':at==bt})
result={'final_pdf':str(new/'main.pdf'),'sha256':hashlib.sha256((new/'main.pdf').read_bytes()).hexdigest(),
        'image_comparison_scope':'Full144dpi page width, pixels y232..1584; excludes top printed page-number header only. No rescaling.',
        'pages':rows}
(HERE/'CAD_FINAL_PAGE_COMPARISON.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
