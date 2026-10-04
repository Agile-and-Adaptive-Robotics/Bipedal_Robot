"""Show unified diffs (repo=o808 as -old, zip as +new) for the conflicted chapters."""
import difflib
from pathlib import Path

ZIP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002")
REP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\o808_export")

def norm(p):
    return [ln.rstrip() for ln in p.read_text(encoding="utf8", errors="replace").splitlines()
            if "\\graphicspath" not in ln]

for rel in ["chapters/15-background.tex", "chapters/20-methods.tex", "chapters/40-discussion.tex",
            "chapters/50-futurework.tex", "chapters/60-conclusion.tex", "chapters/92-AppendixA.tex",
            "chapters/94-AppendixC.tex", "chapters/02-abstract.tex", "chapters/06-ack.tex",
            "chapters/10-introduction.tex"]:
    print(f"\n{'#'*30} {rel} {'#'*30}")
    diff = list(difflib.unified_diff(norm(REP / rel), norm(ZIP / rel),
                                     fromfile="repo", tofile="zip", lineterm="", n=1))
    for d in diff:
        if d.startswith(("---", "+++")):
            continue
        # print removed(repo-only) and added(zip-only) with markers
        print(d if len(d) < 300 else d[:300] + " ...[TRUNC]")
