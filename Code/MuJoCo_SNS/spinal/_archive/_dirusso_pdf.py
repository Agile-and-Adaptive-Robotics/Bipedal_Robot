"""Download the Di Russo 2023 JNE PDF from the AARL Zotero group and
extract every sentence about the foot-contact model (Ben's question:
heel contact + two MTP contacts?)."""
import io
import re
import sys
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

KEY = "34tQExoeiKS1yRe3Z1fbTF5z"
URL = ("https://api.zotero.org/groups/735051/items/7E7F6Q7I/file"
       f"?key={KEY}&v=3")

req = urllib.request.Request(URL)
with urllib.request.urlopen(req, timeout=60) as r:
    pdf = r.read()
print(f"pdf bytes: {len(pdf)}")
with open("_dirusso_2023.pdf", "wb") as f:
    f.write(pdf)

from pypdf import PdfReader

reader = PdfReader("_dirusso_2023.pdf")
text = "\n".join((p.extract_text() or "") for p in reader.pages)
print(f"pages: {len(reader.pages)}, chars: {len(text)}")

# sentences mentioning the contact model
sents = re.split(r"(?<=[.!?])\s+", text.replace("\n", " "))
hits = 0
pat = re.compile(r"contact|heel|metatarsal|MTP|calcan|sole|ground "
                 r"reaction", re.I)
for s in sents:
    if pat.search(s) and re.search(r"foot|heel|MTP|metatarsal|contact "
                                   r"point|sole", s, re.I):
        hits += 1
        print(f"* {s.strip()[:400]}")
        if hits > 40:
            break
