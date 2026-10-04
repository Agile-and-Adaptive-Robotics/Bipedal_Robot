"""Compare Overleaf zip chapters vs local ProofFinal chapters, ignoring graphics-path lines.
Reports: per-file line counts, and whether the ONLY differences are graphicspath/includegraphics lines.
"""
import difflib
import re
from pathlib import Path

ZIP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002")
PF = Path(r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\ProofFinal")

FILES = ["main.tex"] + [f"chapters/{n}" for n in [
    "01-copyright.tex", "02-abstract.tex", "04-dedication.tex", "06-ack.tex",
    "10-introduction.tex", "15-background.tex", "20-methods.tex", "30-results.tex",
    "40-discussion.tex", "50-futurework.tex", "60-conclusion.tex", "90-bib.tex",
    "92-AppendixA.tex", "93-AppendixB.tex", "94-AppendixC.tex"]]

GRAPHICS_RE = re.compile(r"(\\graphicspath|\\includegraphics)")

for rel in FILES:
    z = (ZIP / rel).read_text(encoding="utf-8", errors="replace").splitlines()
    p = (PF / rel).read_text(encoding="utf-8", errors="replace").splitlines()
    # normalize: drop graphicspath lines and bare path prefixes in includegraphics
    def norm(lines):
        out = []
        for ln in lines:
            if "\\graphicspath" in ln:
                continue
            ln = re.sub(r"\\includegraphics(\[[^\]]*\])?\{[^}]*/", r"\\includegraphics\1{", ln)
            out.append(ln.rstrip())
        return out
    zn, pn = norm(z), norm(p)
    if zn == pn:
        print(f"CONTENT-SAME  {rel}  ({len(zn)} vs {len(pn)} lines)")
        continue
    diff = list(difflib.unified_diff(pn, zn, lineterm="", n=0))
    changed = [d for d in diff if d.startswith(("+", "-")) and not d.startswith(("+++", "---"))]
    # how many of the changed lines are graphics-related?
    non_graphics = [d for d in changed if not GRAPHICS_RE.search(d)]
    print(f"CONTENT-DIFF  {rel}  zip={len(zn)} pf={len(pn)} changed={len(changed)} non-graphics-changed={len(non_graphics)}")
