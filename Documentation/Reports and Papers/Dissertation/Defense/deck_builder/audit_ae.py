"""Assertion-evidence audit of the defense deck (see skill: assertion-evidence).

Thresholds follow the official Alley 16:9 template, scaled x4/3 from its 10.00x5.62 in
canvas to this 13.333x7.5 canvas: headline 26-28pt -> floor 32pt; call-outs 18pt ->
floor 24pt for non-assertion emphasis text; citation 14pt -> floor 12pt absolute
(unprojectable below); evidence region ~80% of slide area -> picture+chart area floor
45% (cards/diagrams legitimately replace pictures on structure slides)."""
import math, os, sys
from pptx import Presentation
from pptx.util import Emu

PPTX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "Bolen_Dissertation_Defense.pptx")
SW, SH = 13.333 * 914400, 7.5 * 914400
AREA = SW * SH

# slide numbers exempt from the picture-area floor (title, pure-text structural)
NO_PIC_OK = {1, 42, 43}

prs = Presentation(PPTX)
print(f"slides: {len(prs.slides._sldIdLst)}")
rows = []
for i, slide in enumerate(prs.slides, 1):
    pics = 0.0
    texts = []  # (area_of_box, words, size, text)
    smallest = 99.0
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
                        smallest = min(smallest, r.font.size.pt)
            words = len(t.split())
            try:
                barea = (sh.width * sh.height) / AREA
            except Exception:
                barea = 0
            texts.append((barea, words, size, t.replace("\n", " | ")[:100]))
    if not texts:
        rows.append((i, 0.0, 0, 0, "", 0, 99.0))
        continue
    # assertion = largest-font text box (ties -> widest)
    a = max(texts, key=lambda x: (x[2], x[0]))
    achars = len(a[3])
    # est lines: assume ~0.5em avg char width over the 12.65in assertion span
    cpl = max(10.0, 12.65 * 72 / (0.5 * max(a[2], 1)))
    alines = math.ceil(achars / cpl) if a[2] else 0
    other_words = sum(w for (_, w, _, _) in texts) - a[1]
    rows.append((i, pics, other_words, a[2], a[3][:88], alines, smallest))

print(f"{'#':>3} {'pic%':>5} {'wrds':>4} {'asz':>4} {'ln':>3} {'min':>4}  assertion")
for r in sorted(rows, key=lambda x: (-(x[1] < 0.45), x[2])):
    flag = ""
    if r[1] < 0.45 and r[0] not in NO_PIC_OK: flag += " LOWPIC"
    if r[3] < 32 and r[0] not in NO_PIC_OK: flag += " SMALLASSERT"
    if r[5] > 2: flag += " LONGASSERT"
    if r[6] < 11 and r[0] not in NO_PIC_OK: flag += " TINYTEXT"
    print(f"{r[0]:>3} {r[1]*100:>4.0f}% {r[2]:>4} {r[3]:>4} {r[5]:>3} {r[6]:>4.0f}  {r[4]}{flag}")
