"""Scan the ORIGINAL Overleaf zip chapters for trailing comment blocks (advisor/colleague
notes) so we know what exists and where."""
from pathlib import Path

ZIP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\overleaf_20261002\chapters")
for p in sorted(ZIP.glob("*.tex")):
    ls = p.read_text(encoding="utf8", errors="replace").splitlines()
    # walk from the end over blank/comment-only tail
    i = len(ls)
    while i > 0 and (not ls[i - 1].strip() or ls[i - 1].lstrip().startswith("%")):
        i -= 1
    tail = ls[i:]
    comments = [l for l in tail if l.lstrip().startswith("%")]
    if comments:
        print(f"\n### {p.name}: {len(comments)} trailing comment lines (of {len(tail)} tail lines)")
        for c in comments[:12]:
            print("   ", c.strip()[:150])
    else:
        print(f"{p.name}: no trailing comments")
