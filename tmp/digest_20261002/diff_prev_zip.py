"""Diff PREVIOUS Overleaf zip (git 935ac739) vs CURRENT zip (working tree, committed 89b9cbef).
Shows per-chapter line deltas and prints the actual changed blocks (unified, n=1) so erasures
and insertions are visible. Content-level only (graphicspath lines filtered)."""
import difflib
import re
from pathlib import Path

OLD = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_prev")
NEW = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002")

FILES = ["main.tex"] + [f"chapters/{n}" for n in [
    "01-copyright.tex", "02-abstract.tex", "04-dedication.tex", "06-ack.tex",
    "10-introduction.tex", "15-background.tex", "20-methods.tex", "30-results.tex",
    "40-discussion.tex", "50-futurework.tex", "60-conclusion.tex", "90-bib.tex",
    "92-AppendixA.tex", "93-AppendixB.tex", "94-AppendixC.tex"]]

def norm(p):
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    return [ln.rstrip() for ln in lines if "\\graphicspath" not in ln]

report = []
for rel in FILES:
    o, n = norm(OLD / rel), norm(NEW / rel)
    if o == n:
        report.append(f"UNCHANGED {rel}")
        continue
    diff = list(difflib.unified_diff(o, n, lineterm="", n=1))
    minus = sum(1 for d in diff if d.startswith("-") and not d.startswith("---"))
    plus = sum(1 for d in diff if d.startswith("+") and not d.startswith("+++"))
    report.append(f"CHANGED   {rel}  -{minus} +{plus}  (old {len(o)} -> new {len(n)} lines)")
print("\n".join(report))
print("\n================ DIFF DETAILS (files with removals) ================")
for rel in FILES:
    o, n = norm(OLD / rel), norm(NEW / rel)
    if o == n:
        continue
    diff = list(difflib.unified_diff(o, n, lineterm="", n=1))
    minus = [d for d in diff if d.startswith("-") and not d.startswith("---")]
    if not minus:
        continue
    print(f"\n########## {rel} ##########")
    for d in diff:
        if d.startswith(("---", "+++")):
            continue
        print(d)
