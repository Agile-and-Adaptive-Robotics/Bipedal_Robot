"""Print windows around keyword hits in an extracted paper text.

Usage: _kw_read.py <txt> <kw1,kw2,...> <window>
"""
import io
import re
import sys

T = io.open(sys.argv[1], encoding="utf-8", errors="replace").read()
KWS = sys.argv[2].split(",")
WIN = int(sys.argv[3]) if len(sys.argv) > 3 else 700

shown = set()
for kw in KWS:
    for m in list(re.finditer(re.escape(kw), T, re.I))[:4]:
        a = max(0, m.start() - WIN // 2)
        key = a // 150
        if key in shown:
            continue
        shown.add(key)
        print(f"--- [{kw}] @{m.start()} ---")
        print(T[a:a + WIN].replace("\n", " "))
        print()
