"""Download Shevtsova 2026 (user lib 6LM33WX6) + Shinohara 2025 (group
4PXKC5UR) PDFs and extract full text for the replication passes."""
import io
import sys
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
KEY = "34tQExoeiKS1yRe3Z1fbTF5z"
PDFS = [
    ("shevtsova_2026_elife107480.pdf",
     "https://api.zotero.org/users/631450/items/6LM33WX6/file"),
    ("shinohara_2025_biophiv687930.pdf",
     "https://api.zotero.org/groups/735051/items/4PXKC5UR/file"),
]
for name, url in PDFS:
    req = urllib.request.Request(f"{url}?key={KEY}&v=3")
    try:
        with urllib.request.urlopen(req, timeout=90) as r:
            data = r.read()
        open(name, "wb").write(data)
        print(f"downloaded {name}: {len(data)} bytes")
    except Exception as e:
        print(f"{name} FAILED: {e}")

from pypdf import PdfReader
for name, _ in PDFS:
    try:
        r = PdfReader(name)
        txt = "\n".join((p.extract_text() or "")
                        for p in r.pages)
        out = name.replace(".pdf", "_fulltext.txt")
        open(out, "w", encoding="utf-8").write(txt)
        print(f"extracted {out}: {len(r.pages)} pages, "
              f"{len(txt)} chars")
    except Exception as e:
        print(f"extract {name} FAILED: {e}")
