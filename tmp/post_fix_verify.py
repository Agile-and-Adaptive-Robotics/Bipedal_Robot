import collections
import re

import fitz

doc = fitz.open(r"D:\temp\diss_build_live_20261022\main.pdf")
print("pages:", len(doc))

log = open(r"D:\Github\Bipedal_Robot\tmp\build_live5.log", encoding="utf-8", errors="replace").read()
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

# spot-verify the edited sentences appear
for needle in (
    "Appendix",
    "just below the low end",
    "spiking-neuron mirrors",
    "improved first-pass robot design is now achievable",
    "32 running trials from eight subjects",
):
    hits = [i + 1 for i in range(len(doc)) if needle in doc[i].get_text()]
    print(f"{needle!r}: pages {hits[:4]}")
