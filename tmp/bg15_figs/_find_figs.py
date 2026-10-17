"""Locate figure-caption pages in extracted PDF text (txtwrite, \f page breaks)."""
import re
import sys

for name in sys.argv[1:]:
    txt = open(name, encoding="utf-8", errors="replace").read()
    pages = txt.split("\f")
    print(f"=== {name}: {len(pages)} pages ===")
    for i, pg in enumerate(pages, start=1):
        hits = re.findall(
            r"(Figure\s+\d+[.:][^\n]{0,110}|Fig\.?\s*\d+[.:][^\n]{0,110})",
            pg)
        if hits:
            for h in hits[:6]:
                print(f"  p{i}: {h.strip()}")
