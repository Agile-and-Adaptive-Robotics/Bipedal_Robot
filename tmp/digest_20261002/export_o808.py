"""Export all ProofFinal .tex from git 80873647 (origin tip), then normalized-diff against the
Overleaf zip chapters (graphicspath lines stripped) to enumerate repo-only vs zip-only content."""
import re
import subprocess
import difflib
from pathlib import Path

GIT = r"C:\Users\Ben Bolen\AppData\Local\GitHubDesktop\app-3.6.5\resources\app\git\cmd\git.exe"
REPO = r"D:\Github\Bipedal_Robot"
ZIP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002")
OUT = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\o808_export")
OUT.mkdir(exist_ok=True)

FILES = ["main.tex"] + [f"chapters/{n}" for n in [
    "01-copyright.tex", "02-abstract.tex", "04-dedication.tex", "06-ack.tex",
    "10-introduction.tex", "15-background.tex", "20-methods.tex", "30-results.tex",
    "40-discussion.tex", "50-futurework.tex", "60-conclusion.tex", "90-bib.tex",
    "92-AppendixA.tex", "93-AppendixB.tex", "94-AppendixC.tex"]]

for rel in FILES:
    blob = subprocess.run([GIT, "show", f"80873647:Documentation/Reports and Papers/Dissertation/ProofFinal/{rel}"],
                          cwd=REPO, capture_output=True, check=True).stdout
    (OUT / rel).parent.mkdir(parents=True, exist_ok=True)
    (OUT / rel).write_bytes(blob)

def norm(p):
    return [ln.rstrip() for ln in p.read_text(encoding="utf8", errors="replace").splitlines()
            if "\\graphicspath" not in ln]

for rel in FILES:
    z, o = norm(ZIP / rel), norm(OUT / rel)
    if z == o:
        print(f"SAME      {rel}")
        continue
    diff = [d for d in difflib.unified_diff(o, z, lineterm="", n=0)
            if d.startswith(("+", "-")) and not d.startswith(("+++", "---"))]
    minus = sum(1 for d in diff if d.startswith("-"))
    plus = sum(1 for d in diff if d.startswith("+"))
    print(f"DIFF      {rel}  repo-only(-)={minus}  zip-only(+)={plus}  (repo {len(o)} / zip {len(z)} lines)")
