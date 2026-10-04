"""Round 2: hyperelastic/FEA vocabulary + docx/pdf/doc filenames in likely author folders."""
import os
import re

ROOTS = [
    r"D:\Github\Academic_Writing",
    r"D:\Github\Bipedal_Robot\Documentation",
    r"D:\Github\Bipedal_Robot\Notes",
    r"C:\Users\Ben Bolen\Desktop",
    r"C:\Users\Ben Bolen\Documents",
    r"C:\Users\Ben Bolen\Zotero\storage",
]
PATS = re.compile(r"(hyperelastic|mooney|ogden|yeoh|finite[- ]element|\bFEA\b|ansys|abaqus|comsol|solidworks simulation)", re.I)
CTX = re.compile(r"(bpa|braided|bladder|pneumatic|mckibben|muscle|actuator|festo)", re.I)
DOC_EXT = {".docx", ".doc", ".pdf", ".pptx", ".odt", ".rtf"}
TEXT_EXT = {".md", ".txt", ".tex", ".bib"}

name_hits, ctx_hits = [], []
for root in ROOTS:
    if not os.path.isdir(root):
        continue
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if not d.startswith(".") and d not in {"node_modules", "pdf_staging"}]
        for f in fns:
            p = os.path.join(dp, f)
            ext = os.path.splitext(f)[1].lower()
            if PATS.search(f) or (ext in DOC_EXT and re.search(r"fea|finite|bpa|character|future|propos", f, re.I)):
                name_hits.append(p)
            if ext in TEXT_EXT and os.path.getsize(p) < 3_000_000:
                try:
                    t = open(p, encoding="utf8", errors="ignore").read()
                except OSError:
                    continue
                for m in PATS.finditer(t):
                    seg = t[max(0, m.start() - 300):m.start() + 300]
                    if CTX.search(seg):
                        ctx_hits.append((p, t.count("\n", 0, m.start()) + 1, seg.replace("\n", " ")[:200]))
                        break

print("=== DOC/FILENAME HITS ===")
for h in name_hits[:50]:
    print(h)
print("\n=== TEXT HITS ===")
seen = set()
for p, ln, seg in ctx_hits:
    if p in seen:
        continue
    seen.add(p)
    print(f"{p}:{ln}\n    ...{seg}...")
