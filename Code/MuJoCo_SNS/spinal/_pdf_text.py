"""Extract text from a PDF with pypdf; print hits for V-class keywords.

Usage: _pdf_text.py <pdf_path> <out_txt> [grep...]
"""
import io
import re
import sys

from pypdf import PdfReader


def main():
    pdf, out = sys.argv[1], sys.argv[2]
    r = PdfReader(pdf)
    chunks = []
    for p in r.pages:
        try:
            chunks.append(p.extract_text() or "")
        except Exception:
            chunks.append("")
    t = "\n".join(chunks)
    io.open(out, "w", encoding="utf-8").write(t)
    print(f"{len(r.pages)} pages -> {out} ({len(t)} chars)")
    for kw in sys.argv[3:]:
        idxs = [m.start() for m in re.finditer(re.escape(kw), t, re.I)][:10]
        print(f"{kw!r}: {idxs}")


if __name__ == "__main__":
    main()
