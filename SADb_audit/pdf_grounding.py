"""PDF-first grounding: extract full text for flagged papers with attached
PDFs (Ben's rule, 2026-10-01: the actual PDF beats the abstract when attached).

Sources, in order: the local Zotero storage (matched via pdf_inventory.csv by
DOI, then title), then Airtable attachment URLs (NOT used here — needs a
schema read; add when API budget allows). Output per paper:
  curation_queue/pdf_grounding/<recordId>.txt  (first 20k chars of text)
  curation_queue/pdf_grounding/_index.json     (what matched, from where)
Stdlib + pypdf. Usage: myo python pdf_grounding.py [n_papers]
"""
import csv, json, os, re, sys
from pypdf import PdfReader

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "curation_queue", "pdf_grounding")
os.makedirs(OUT, exist_ok=True)

records = {r["id"]: r for r in json.load(
    open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))}
flags = json.load(open(os.path.join(HERE, "curation_out", "_flags_flat.json"), encoding="utf-8"))

inv = []
for row in csv.DictReader(open(os.path.join(HERE, "pdf_inventory.csv"), encoding="utf-8-sig")):
    if row["file_exists"] == "True" and os.path.exists(row["local_path"]):
        inv.append(row)
by_doi = {r["doi"].strip().lower(): r for r in inv if r["doi"]}
print(f"local PDFs available: {len(inv)}")


def norm(t):
    return re.sub(r"[^a-z0-9]+", " ", (t or "").lower()).strip()


by_title = {norm(r["title"]): r for r in inv}

n_target = int(sys.argv[1]) if len(sys.argv) > 1 else 5
done = 0
index = {}
if os.path.exists(os.path.join(OUT, "_index.json")):
    index = json.load(open(os.path.join(OUT, "_index.json"), encoding="utf-8"))

for f in flags:
    rid = f["id"]
    r = records.get(rid, {})
    if not r.get("has_pdf") or rid in index:
        continue
    if done >= n_target:
        break
    hit = by_doi.get((r.get("doi") or "").strip().lower()) or by_title.get(norm(r["title"]))
    if not hit:
        # fuzzy: title prefix match
        nt = norm(r["title"])[:40]
        for cand in inv:
            if norm(cand["title"]).startswith(nt):
                hit = cand
                break
    if not hit:
        index[rid] = {"status": "no-local-pdf", "title": r["title"]}
        continue
    try:
        reader = PdfReader(hit["local_path"])
        text = " ".join((p.extract_text() or "") for p in reader.pages[:8])
        text = re.sub(r"\s+", " ", text)[:20000]
        open(os.path.join(OUT, rid + ".txt"), "w", encoding="utf-8").write(text)
        index[rid] = {"status": "extracted", "chars": len(text),
                      "pdf": hit["local_path"], "via": hit["doi"] and "doi" or "title"}
        done += 1
    except Exception as e:
        index[rid] = {"status": f"extract-failed: {str(e)[:60]}", "pdf": hit["local_path"]}

json.dump(index, open(os.path.join(OUT, "_index.json"), "w"), indent=1)
n_ok = sum(1 for v in index.values() if v["status"] == "extracted")
n_no = sum(1 for v in index.values() if v["status"] == "no-local-pdf")
print(f"extracted {done} this run | total extracted {n_ok} | no local PDF {n_no}")
