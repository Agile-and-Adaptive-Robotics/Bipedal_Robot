"""Hunt for Ben's 'advanced FEA for BPA characterization' document.
Filename pass + text-content pass over likely locations."""
import os
import re

ROOTS = [
    r"D:\Github\Bipedal_Robot\Documentation",
    r"D:\Github\Bipedal_Robot\Code\Matlab\Mesh_Optimization",
    r"C:\Users\Ben Bolen\Desktop",
    r"C:\Users\Ben Bolen\Documents",
    r"C:\Users\Ben Bolen\Downloads",
    r"D:\temp",
    r"D:\Github",
]
NAME_PAT = re.compile(r"(fea|finite[-_ ]?element|ansys|abaqus|comsol)", re.I)
TEXT_PAT = re.compile(r"(finite[- ]element|\bFEA\b|ansys|abaqus)", re.I)
CTX_PAT = re.compile(r"(bpa|braided|pneumatic|muscle|bladder|actuator|festo|mcKibben|mckibben)", re.I)
SKIP_DIRS = {"node_modules", ".git", "pdf_staging", "AppData", "__pycache__", "gait_lib_staging"}
TEXT_EXT = {".md", ".txt", ".tex", ".bib", ".py", ".m", ".csv", ".json"}

name_hits, ctx_hits = [], []
for root in ROOTS:
    if not os.path.isdir(root):
        continue
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if d not in SKIP_DIRS and not d.startswith(".")]
        if len(name_hits) + len(ctx_hits) > 400:
            break
        for f in fns:
            p = os.path.join(dp, f)
            if NAME_PAT.search(f):
                name_hits.append(p)
            ext = os.path.splitext(f)[1].lower()
            if ext in TEXT_EXT and os.path.getsize(p) < 3_000_000:
                try:
                    t = open(p, encoding="utf8", errors="ignore").read()
                except OSError:
                    continue
                for m in TEXT_PAT.finditer(t):
                    seg = t[max(0, m.start() - 400):m.start() + 400]
                    if CTX_PAT.search(seg):
                        line_no = t.count("\n", 0, m.start()) + 1
                        ctx_hits.append((p, line_no, seg.replace("\n", " ")[:220]))
                        break

print("=== FILENAME HITS ===")
for h in name_hits[:60]:
    print(h)
print("\n=== CONTENT HITS (FEA near BPA/actuator terms) ===")
seen = set()
for p, ln, seg in ctx_hits:
    if p in seen:
        continue
    seen.add(p)
    print(f"{p}:{ln}\n    ...{seg}...")
    if len(seen) > 40:
        break
