"""Scan the REPO (o808 export) chapters for trailing comment blocks — the advisor notes
that the Overleaf-side pastes may have dropped."""
from pathlib import Path

REP = Path(r"C:\Users\Ben Bolen\AppData\Local\Temp\o808_export\chapters")
for p in sorted(REP.glob("*.tex")):
    ls = p.read_text(encoding="utf8", errors="replace").splitlines()
    i = len(ls)
    while i > 0 and (not ls[i - 1].strip() or ls[i - 1].lstrip().startswith("%")):
        i -= 1
    comments = [l for l in ls[i:] if l.lstrip().startswith("%")]
    if comments:
        print(f"\n### {p.name}: {len(comments)} trailing comment lines (tail starts at line {i+1} of {len(ls)})")
        for c in comments[:15]:
            print("   ", c.strip()[:160])
        if len(comments) > 15:
            print(f"    ... +{len(comments)-15} more")
    # also scan for inline '% [advisor' / ajh markers anywhere
    hits = [k + 1 for k, l in enumerate(ls) if ("ajh" in l.lower() and l.lstrip().startswith("%"))]
    if hits:
        print(f"{p.name}: % comments mentioning 'ajh' at lines {hits[:10]}")
