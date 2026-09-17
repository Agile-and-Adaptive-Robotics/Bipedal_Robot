"""Extract McCrea 1980 text; grep RC mutual inhibition + IaIN statements."""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
from pypdf import PdfReader

path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "_mccrea1980_renshaw.pdf")
r = PdfReader(path)
text = "\n".join((p.extract_text() or "") for p in r.pages)
low = text.lower()
print(f"{len(r.pages)} pages")
for t in ["inhibit other Renshaw", "inhibition of other Renshaw",
          "Renshaw cells inhibit", "mutual", "each other",
          "Ia inhibitory interneur", "recurrent inhibitory postsynaptic",
          "during fictive locomotion Renshaw"]:
    idx = low.find(t.lower())
    if idx >= 0:
        snip = text[max(0, idx - 180): idx + 300].replace("\n", " ")
        print(f"[{t}] ...{snip}...\n")
    else:
        print(f"[{t}] not found")
