import sys

import fitz

doc = fitz.open(sys.argv[1])
pno = int(sys.argv[2]) - 1
page = doc[pno]
print("page rect:", page.rect)
for b in page.get_text("dict")["blocks"]:
    for l in b.get("lines", []):
        for s in l.get("spans", []):
            x1 = s["bbox"][2]
            if x1 > 545:
                print(f"span x1={x1:.1f} y={s['bbox'][1]:.1f} size={s['size']:.1f} text={s['text'][:80]!r}")
# also images/drawings bboxes
for img in page.get_images(full=True):
    try:
        r = page.get_image_bbox(img)
        print("image bbox:", r)
    except Exception as e:
        print("img err", e)
for d in page.get_drawings():
    r = d["rect"]
    if r.x1 > 545:
        print(f"drawing x1={r.x1:.1f} y0={r.y0:.1f} y1={r.y1:.1f} type={d['type']}")
