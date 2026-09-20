"""Reconcile rest-import queue rows against the exported corpus state."""
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))
exp = json.load(open(os.path.join(HERE, "export", "sadb_export.json"), encoding="utf-8"))
by_doi = {}
by_title = {}
for r in exp:
    d = (r["doi"] or "").strip().lower().replace("https://doi.org/", "").replace("http://doi.org/", "")
    if d:
        by_doi[d] = r
    t = "".join(ch for ch in r["title"].lower() if ch.isalnum())
    by_title[t] = r

rows = list(csv.DictReader(open(os.path.join(HERE, "airtable_rest_import_clean.csv"), encoding="utf-8-sig")))
lo, hi = (41, 70)
need = []
for i, row in enumerate(rows[lo - 1:hi], start=lo):
    doi = (row["DOI"] or "").strip().lower()
    rec = by_doi.get(doi)
    how = "doi"
    if rec is None:
        t = "".join(ch for ch in row["title"].lower() if ch.isalnum())
        rec = by_title.get(t)
        how = "title"
    if rec is None:
        print(f"{i:3d} {row['zotero_key']} *** NO MATCH *** | {row['title'][:60]}")
        continue
    flag = "NEEDS" if not rec["has_notes"] else "done "
    if not rec["has_notes"]:
        need.append((i, row, rec))
    print(f"{i:3d} {row['zotero_key']} {flag} note={int(rec['has_notes'])} an={len(rec['animals'])} "
          f"fb={len(rec['feedback'])} pdf={int(rec['has_pdf'])} id={rec['id']} ({how}) | {row['title'][:55]}")

print(f"\n{len(need)} of rows {lo}-{hi} still need curation")
json.dump([{"queue_row": i, "zotero_key": row["zotero_key"], "doi": row["DOI"],
            "title": row["title"], "author": row["author"], "year": row["year"],
            "airtable_id": rec["id"]}
           for i, row, rec in need],
          open(os.path.join(HERE, "batch6", "batch6_targets.json"), "w", encoding="utf-8"),
          ensure_ascii=False, indent=1)
