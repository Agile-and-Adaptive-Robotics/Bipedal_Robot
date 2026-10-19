"""Assertion-evidence audit of the defense deck (see skill: assertion-evidence)."""
import math, os, sys
from pptx import Presentation
from pptx.util import Emu

PPTX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Bolen_Dissertation_Defense.pptx")
SW, SH = 13.333 * 914400, 7.5 * 914400
AREA = SW * SH

prs = Presentation(PPTX)
print(f"slides: {len(prs.slides.__iter__.__self__._sldIdLst)}" if False else f"slides: {len(prs.slides._sldIdLst)}")
rows = []
for i, slide in enumerate(prs.slides, 1):
    pics = 0.0
    texts = []  # (area_of_box, words, size, text)
    for sh in slide.shapes:
        if sh.shape_type == 13:  # PICTURE
            try:
                pics += (sh.width * sh.height) / AREA
            except Exception:
                pass
        elif sh.has_chart if hasattr(sh, "has_chart") else False:
            try:
                pics += (sh.width * sh.height) / AREA
            except Exception:
                pass
        elif sh.has_text_frame:
            t = sh.text_frame.text.strip()
            if not t:
                continue
            size = 0
            for p in sh.text_frame.paragraphs:
                for r in p.runs:
                    if r.font.size:
                        size = max(size, r.font.size.pt)
            words = len(t.split())
            try:
                barea = (sh.width * sh.height) / AREA
            except Exception:
                barea = 0
            texts.append((barea, words, size, t.replace("\n", " | ")[:100]))
    if not texts:
        rows.append((i, 0.0, 0, 0, "", 0))
        continue
    # assertion = largest-font text box (ties -> widest)
    a = max(texts, key=lambda x: (x[2], x[0]))
    achars = len(a[3])
    # est lines: assume ~0.5em avg char width
    cpl = max(10.0, (sh_width := 12.25) * 72 / (0.5 * max(a[2], 1)))
    alines = math.ceil(achars / cpl) if a[2] else 0
    other_words = sum(w for (_, w, _, _) in texts) - a[1]
    rows.append((i, pics, other_words, a[2], a[3][:88], alines))

print(f"{'#':>3} {'pic%':>5} {'wrds':>4} {'asz':>4} {'ln':>3}  assertion")
for r in sorted(rows, key=lambda x: (-(x[1] < 0.45), x[2])):
    flag = ""
    if r[1] < 0.45: flag += " LOWPIC"
    if r[3] < 24: flag += " SMALLASSERT"
    print(f"{r[0]:>3} {r[1]*100:>4.0f}% {r[2]:>4} {r[3]:>4} {r[5]:>3}  {r[4]}{flag}")
