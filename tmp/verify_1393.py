import re

import fitz

doc = fitz.open(r"D:\temp\diss_build_live_20261022\main.pdf")
print("pages:", len(doc))
t = " ".join(doc[i].get_text() for i in range(140, 150))
t = re.sub(r"\s+", " ", t)
print("1,393 in prose:", "1,393 records before the prune" in t.replace("1;393", "1,393") or "1,393" in t)
log = open(r"D:\Github\Bipedal_Robot\tmp\build_live9.log", encoding="utf-8", errors="replace").read()
tail = log[log.rfind("=== PASS 4 ==="):]
print("final-pass undefined:", len(re.findall(r"[Uu]ndefined", tail)))
print("final-pass missing:", tail.count("not found"))
