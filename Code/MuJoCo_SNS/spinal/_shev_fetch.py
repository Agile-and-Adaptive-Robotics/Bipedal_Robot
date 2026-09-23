"""Fetch Shevtsova 2026 eLife 107480 PDF (open access) + extract."""
import io
import sys
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
for u in ("https://cdn.elifesciences.org/articles/107480.pdf",
          "https://elifesciences.org/articles/107480.pdf"):
    try:
        d = urllib.request.urlopen(u, timeout=90).read()
        open("shevtsova_2026_elife107480.pdf", "wb").write(d)
        print("OK", u, len(d), "bytes")
        break
    except Exception as e:
        print("fail", u, e)
else:
    sys.exit(1)

from pypdf import PdfReader
r = PdfReader("shevtsova_2026_elife107480.pdf")
t = "\n".join((p.extract_text() or "") for p in r.pages)
open("shevtsova_2026_elife107480_fulltext.txt", "w",
     encoding="utf-8").write(t)
print(f"extracted {len(r.pages)} pages, {len(t)} chars")
