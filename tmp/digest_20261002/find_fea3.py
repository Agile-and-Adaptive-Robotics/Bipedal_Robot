"""Round 3: prospectus tex grep + docx internals scan for FEA-for-BPA content."""
import os
import re
import zipfile

PROS = r"D:\Github\Academic_Writing\Ben_Bolen\Bolen_prospectus"
PAT = re.compile(r"(finite[- ]element|\bFEA\b|hyperelastic|mooney|ogden|ansys|abaqus|comsol| bladder simulation|advanced simulation)", re.I)

print("=== prospectus tex/md hits ===")
for dp, _, fns in os.walk(PROS):
    for f in fns:
        if os.path.splitext(f)[1].lower() in {".tex", ".md", ".txt", ".bib"}:
            p = os.path.join(dp, f)
            try:
                for i, ln in enumerate(open(p, encoding="utf8", errors="ignore"), 1):
                    if PAT.search(ln):
                        print(f"{p}:{i}: {ln.strip()[:160]}")
            except OSError:
                pass

print("\n=== docx scan (Desktop/Documents/Downloads/Academic_Writing/Bipedal Documentation) ===")
ROOTS = [r"C:\Users\Ben Bolen\Desktop", r"C:\Users\Ben Bolen\Documents",
         r"D:\Github\Academic_Writing", r"D:\Github\Bipedal_Robot\Documentation"]
for root in ROOTS:
    if not os.path.isdir(root):
        continue
    for dp, dns, fns in os.walk(root):
        dns[:] = [d for d in dns if not d.startswith(".") and d not in {"node_modules", "pdf_staging", "Zotero"}]
        for f in fns:
            if not f.lower().endswith((".docx", ".docm")):
                continue
            p = os.path.join(dp, f)
            try:
                with zipfile.ZipFile(p) as z:
                    xml = z.read("word/document.xml").decode("utf8", errors="ignore")
            except Exception:
                continue
            m = PAT.search(xml)
            if m and re.search(r"(bpa|braided|bladder|pneumatic|actuator|muscle)", xml, re.I):
                seg = re.sub(r"<[^>]+>", " ", xml[max(0, m.start() - 300):m.start() + 300])
                print(f"{p}\n    ...{' '.join(seg.split())[:240]}...")
