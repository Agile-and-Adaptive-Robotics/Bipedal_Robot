"""Extract text from audit PDFs and grep the connectivity claims that
decide audit rows. Writes lit_fulltext_findings.txt."""
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
from pypdf import PdfReader

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lit_pdfs")

CLAIMS = {
    "mccrea1980_renshaw.pdf": [
        "mutual inhibition of Renshaw cells", "Renshaw cells inhibit each",
        "recurrent inhibition of Ia inhibitory", "Ia inhibitory interneur",
        "inhibits other Renshaw"],
    "pratt1987_iain_renshaw.pdf": [
        "Ia inhibitory interneur", "Renshaw", "recurrent"],
    "perreault1999_grI_extensor.pdf": [
        "excitation replaced", "oligosynaptic excitation", "weight support",
        "group I"],
    "perreault2011_grii_reset.pdf": [
        "reset", "group II", "flexor phase", "extensor phase"],
    "jankowska2010_ib_ii.pdf": [
        "group II", "group Ib", "overlapping", "distinct populations"],
    "dominguez2020_reset_ins.pdf": [
        "candidate", "interneurons", "V2a", "Ib", "non-reciprocal",
        "polysynaptic", "prolong"],
}

for pdf, terms in CLAIMS.items():
    path = os.path.join(OUT, pdf)
    if not os.path.exists(path):
        continue
    try:
        r = PdfReader(path)
        text = "\n".join((p.extract_text() or "") for p in r.pages[:14])
        low = text.lower()
        print(f"===== {pdf} ({len(r.pages)} pages) =====")
        for t in terms:
            idx = low.find(t.lower())
            if idx >= 0:
                snippet = text[max(0, idx - 160): idx + 260].replace("\n", " ")
                print(f"  [{t}] ...{snippet}...")
            else:
                print(f"  [{t}] not found in first 14 pages")
        print()
    except Exception as e:
        print(f"===== {pdf}: extraction ERROR {e}")
