"""Cross-reference pdf_inventory.csv against the Airtable corpus DOIs.
Outputs author_fix/pdfs_from_zotero.csv (corpus papers that HAVE a Zotero PDF,
with staging path) and counts the corpus papers still lacking any PDF.
"""
import csv, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
SAD = os.path.dirname(HERE)

def norm_doi(d):
    if not d:
        return ""
    d = str(d).strip().lower()
    for p in ("https://doi.org/", "http://doi.org/", "doi:"):
        if d.startswith(p):
            d = d[len(p):]
    return d

corpus = {}  # doi -> (source, display author)
with open(os.path.join(SAD, "airtable_papers_slim.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        d = norm_doi(r.get("DOI"))
        if d:
            corpus[d] = (r.get("author", ""), "orig99")
with open(os.path.join(SAD, "airtable_rest_import_clean.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        d = norm_doi(r.get("DOI"))
        if d:
            corpus[d] = (r.get("author", ""), "rest383")
with open(os.path.join(SAD, "batch3", "rr_digest_shortlist.json"), encoding="utf-8") as f:
    for s in json.load(f):
        corpus[norm_doi(s["doi"])] = ("", "digest20")
corpus[norm_doi("10.7554/elife.107480")] = ("", "digest20")

zot = {}  # doi -> list of (library, path, surname, year)
with open(os.path.join(SAD, "pdf_inventory.csv"), encoding="utf-8-sig") as f:
    for r in csv.DictReader(f):
        if r["file_exists"].lower() != "true" or not r["doi"]:
            continue
        zot.setdefault(r["doi"], []).append((r["library"], r["local_path"], r["surname"], r["year"]))

have, missing = [], []
for d, (disp, src) in sorted(corpus.items()):
    if d in zot:
        lib, path, sn, yr = zot[d][0]
        have.append({"doi": d, "source": src, "old_display_author": disp,
                     "copies": len(zot[d]), "libraries": ",".join(sorted({z[0] for z in zot[d]})),
                     "staging_or_local": path})
    else:
        missing.append({"doi": d, "source": src, "old_display_author": disp})

outp = os.path.join(HERE, "pdfs_from_zotero.csv")
with open(outp, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["doi", "source", "old_display_author", "copies", "libraries", "staging_or_local"])
    w.writeheader()
    w.writerows(have)
outm = os.path.join(HERE, "pdfs_missing_no_zotero.csv")
with open(outm, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["doi", "source", "old_display_author"])
    w.writeheader()
    w.writerows(missing)

print("corpus DOI records:", len(corpus))
print("have Zotero PDF:", len(have), "->", outp)
print("missing (need browser/OA):", len(missing), "->", outm)
