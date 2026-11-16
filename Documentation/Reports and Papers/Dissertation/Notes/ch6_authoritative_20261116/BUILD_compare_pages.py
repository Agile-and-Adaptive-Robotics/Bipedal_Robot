"""Find unchanged publication-size page bodies between reviewed and final PDFs.

PNG pixels are read without resampling. Only the standard page-number/header
and blank bottom margin are excluded from equality comparisons. Changed page
bodies still require visual review; equality is not claimed for those pages.
"""
from pathlib import Path
import hashlib, json
from PIL import Image

NOTES=Path(__file__).resolve().parent
old=json.loads((NOTES/'BUILD_RENDER_r6.json').read_text())
new=json.loads((NOTES/'BUILD_RENDER_r7.json').read_text())
crop=(0,239,1224,1450)
def body_sha(path):
    with Image.open(path) as image:
        if image.size!=(1224,1584): raise ValueError('Unexpected render dimensions')
        return hashlib.sha256(image.crop(crop).tobytes()).hexdigest()
old_hashes={}
for page in old['pages']:
    if page['group']=='chapter6': old_hashes.setdefault(body_sha(page['png']),[]).append(page['pdf_page_1based'])
comparisons=[]
for page in new['pages']:
    if page['group']!='chapter6': continue
    matches=old_hashes.get(body_sha(page['png']),[])
    comparisons.append({'final_pdf_page':page['pdf_page_1based'],'pixel_identical_reviewed_r6_pages':matches,
                        'requires_new_visual_review':not bool(matches)})
report={'old_pdf_sha256':old['pdf_sha256'],'final_pdf_sha256':new['pdf_sha256'],
        'exact_pixel_comparison_crop':crop,'resampling':False,'pages':comparisons,
        'new_visual_review_pages':[p['final_pdf_page'] for p in comparisons if p['requires_new_visual_review']]}
(NOTES/'BUILD_PAGE_COMPARISON_r6_r7.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
