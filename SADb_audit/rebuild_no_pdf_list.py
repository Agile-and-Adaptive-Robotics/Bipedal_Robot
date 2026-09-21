"""Rebuild author_fix/remaining_no_pdf.csv as THE definitive no-PDF hunt list.

Schema: doi, record_id, title, hunt_status (verdict from pass 2, 2026-09-20).
Supersedes (and allows deleting): paperhunt_status.csv, pdfs_missing_no_zotero.csv,
pdfs_paywalled.csv, pdfs_to_attach.json, upload_status.csv, fileio_status.csv,
tunnel_status.csv, pdf_hunt_pass2_20260920.csv.
"""
import csv, json, os

SAD = os.path.dirname(os.path.abspath(__file__))
EXP = json.load(open(os.path.join(SAD, "export", "sadb_export.json"), encoding="utf-8"))

# verdict merge order (later overrides): previous hunt list -> pass-2 sweep -> retry
old_statuses, old_dois = {}, set()
with open(os.path.join(SAD, "author_fix", "remaining_no_pdf.csv"), encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        d = row["doi"].strip().lower()
        old_dois.add(d)
        old_statuses[d] = row.get("hunt_status") or row.get("state") or ""
verdict = dict(old_statuses)
for fn in ("pdf_hunt_pass2_20260920.csv", "pdf_hunt_retry_20260920.csv"):
    p = os.path.join(SAD, "author_fix", fn)
    if os.path.exists(p):
        with open(p, encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                verdict[row["doi"].strip().lower()] = row["status"]

rows_out, attached_since = [], 0
for r in EXP:
    d = (r["doi"] or "").strip().lower()
    if not d or r["has_pdf"]:
        continue
    v = verdict.get(d, "not-yet-hunted (new record)")
    if v.startswith("ATTACHED"):
        v = "attach-dropped (Airtable accepted PATCH but never stored the file)"  # silent-drop trap
    rows_out.append({"doi": d, "record_id": r["id"],
                     "title": (r["title"] or "")[:100], "hunt_status": v})

# records a previous hunt list contained but that have a PDF now:
attached_since = sum(1 for r in EXP
                     if (r["doi"] or "").strip().lower() in old_dois and r["has_pdf"])

out = os.path.join(SAD, "author_fix", "remaining_no_pdf.csv")
with open(out, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["doi", "record_id", "title", "hunt_status"])
    w.writeheader()
    w.writerows(rows_out)

from collections import Counter
print(f"definitive list: {len(rows_out)} records -> {out}")
print(Counter(r["hunt_status"].split(':')[0] for r in rows_out))
print(f"on the old Sept-15 list but have a PDF now: {attached_since}")
