"""Extract Di Russo 2023 sections: the five motor primitives (PF
feedforward) and the proprioceptive-feedback / interneuron / Renshaw
rules (for LIT_CIRCUIT_AUDIT.md, Ben 2026-09-21)."""
import io
import re
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

from pypdf import PdfReader

reader = PdfReader("_dirusso_2023.pdf")
pages = [(i + 1, p.extract_text() or "") for i, p in
         enumerate(reader.pages)]

for word in ("primitive", "synerg", "Renshaw", "interneuron",
             "Ia ", "Ib ", "group II", "positive force", "length "
             "feedback", "feedback rule"):
    print(f"===== '{word}' =====")
    n = 0
    for pno, txt in pages:
        flat = re.sub(r"-\s*\n", "", txt)
        flat = re.sub(r"\s+", " ", flat)
        for m in re.finditer(re.escape(word), flat, re.I):
            a, b = max(0, m.start() - 200), min(len(flat), m.end() + 250)
            n += 1
            if n <= 14:
                print(f"  p{pno}: ...{flat[a:b]}...")
    print(f"  ({n} occurrences)\n")
