"""Scan pptx slide XML for FEA/finite-element content."""
import re
import zipfile
import os

PAT = re.compile(r"(finite[- ]element|\bFEA\b|hyperelastic|mooney|ogden|ansys|abaqus|comsol)", re.I)
FILES = [
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Sigma_XI presentation.pptx",
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Bolen_Dissertation_Defense_20260910.pptx",
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Movement and Control of Biomimetic Humanoid Robots.pptx",
    r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\PSU Template.pptx",
]
for p in FILES:
    if not os.path.exists(p):
        print("missing:", p)
        continue
    hits = []
    with zipfile.ZipFile(p) as z:
        for n in z.namelist():
            if n.startswith("ppt/slides/slide") and n.endswith(".xml"):
                xml = z.read(n).decode("utf8", errors="ignore")
                if PAT.search(xml):
                    slide = n.split("/")[-1].replace(".xml", "")
                    text = " ".join(re.sub(r"<[^>]+>", " ", xml).split())
                    hits.append((slide, text[:400]))
    print("\n===", os.path.basename(p), "->", len(hits), "matching slides")
    for s, t in hits:
        print(f"  {s}: {t[:350]}")
