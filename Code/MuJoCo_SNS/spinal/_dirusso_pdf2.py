"""Targeted context search in the Di Russo 2023 PDF: contact geometry
(MTP/metatarsal/forefoot/contact points/spheres)."""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from pypdf import PdfReader

reader = PdfReader("_dirusso_2023.pdf")
text = "\n".join((p.extract_text() or "") for p in reader.pages)
flat = re.sub(r"-\s*\n", "", text)          # de-hyphenate line breaks
flat = re.sub(r"\s+", " ", flat)

for pat in ("MTP", "metatar", "forefoot", "contact point",
            "contact sphere", "three contact", "two contact",
            "heel contact", "foot model", "cutaneous"):
    for m in re.finditer(pat, flat, re.I):
        a, b = max(0, m.start() - 220), min(len(flat), m.end() + 220)
        print(f"--- [{pat}] ...{flat[a:b]}...")
        print()
