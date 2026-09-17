"""Grep the Rybak 2015 full text for Renshaw connectivity statements."""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
T = io.open("_rybak2015_plain.txt", encoding="utf-8", errors="replace").read()
seen = set()
for m in re.finditer(r"Renshaw", T):
    a = max(0, m.start() - 250)
    key = a // 200
    if key in seen:
        continue
    seen.add(key)
    print("..." + T[a:m.start() + 350].replace("\n", " ") + "...\n")
