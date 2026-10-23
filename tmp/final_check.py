import collections
import re

import fitz

doc = fitz.open(r"D:\temp\diss_build_live_20261022\main.pdf")
print("pages:", len(doc))

# locate spike_animatlab table page and check span bounds + rules
for i, page in enumerate(doc):
    t = page.get_text()
    if "synchronous" in t and "11.8" in t and "Two-layer walker" in t:
        print("spike table PDF page:", i + 1)
        worst = 0
        for b in page.get_text("dict")["blocks"]:
            for l in b.get("lines", []):
                if l["dir"] != (1, 0):
                    continue
                for s in l.get("spans", []):
                    worst = max(worst, s["bbox"][2])
        for d in page.get_drawings():
            if d["rect"].y1 - d["rect"].y0 < 1:  # horizontal rule
                worst = max(worst, d["rect"].x1)
        print("worst upright x1:", round(worst, 1), "(text block right edge ~540.7)")
        break

# final-pass warnings
log = open(r"D:\Github\Bipedal_Robot\tmp\build_live4.log", encoding="utf-8", errors="replace").read()
tail = log[log.rfind("=== PASS 4 ==="):]
print("final-pass undefined:", len(re.findall(r"[Uu]ndefined", tail)))
print("final-pass missing files:", tail.count("not found"))
new = collections.Counter(re.findall(r"Overfull \\hbox \(([\d.]+)pt", tail))
print("final-pass overfull total:", sum(new.values()))
base = collections.Counter(
    re.findall(
        r"Overfull \\hbox \(([\d.]+)pt",
        open(r"D:\Github\Bipedal_Robot\tmp\build_base.log", encoding="utf-8", errors="replace").read().rsplit("=== PASS 4 ===", 1)[1],
    )
)
print("added vs baseline:", sum((new - base).values()), "removed:", sum((base - new).values()))
