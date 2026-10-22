import collections
import re

import fitz

doc = fitz.open(r"D:\temp\diss_build_live_20261022\main.pdf")
print("pages:", len(doc))

log = open(r"D:\Github\Bipedal_Robot\tmp\build_live8.log", encoding="utf-8", errors="replace").read()
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

# contradiction fixes present?
t = "\n".join(doc[i].get_text() for i in range(len(doc)))
print("sentinel-above claim:", "sits numerically above every genuine gait score" in t)
print("old wrong claim gone:", "sits below every genuine gait" not in t)
print("DOI qualified:", "unique identifier wherever a record carries one" in t)
print("app primary:", "primary form is the packaged explorer" in t)
print("archived papers in app:", "papers the prune archived from the working base" in t)
print("no 'the author':", "the author" not in t)
