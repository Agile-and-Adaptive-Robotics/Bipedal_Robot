"""Exhaustive contact-model search in Di Russo 2023 (local copy of
Ben's attached PDF): every occurrence of contact/heel/toe/MTP/point
with surrounding context."""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from pypdf import PdfReader

reader = PdfReader("_dirusso_2023.pdf")
pages = [(i + 1, p.extract_text() or "") for i, p in
         enumerate(reader.pages)]

for word in ("contact", "heel", "toe", "MTP", "metatar", "cutaneous",
             "sole"):
    print(f"===== '{word}' =====")
    n = 0
    for pno, txt in pages:
        flat = re.sub(r"-\s*\n", "", txt)
        flat = re.sub(r"\s+", " ", flat)
        for m in re.finditer(word, flat, re.I):
            a, b = max(0, m.start() - 150), min(len(flat), m.end() + 150)
            n += 1
            print(f"  p{pno}: ...{flat[a:b]}...")
    print(f"  ({n} occurrences)\n")
