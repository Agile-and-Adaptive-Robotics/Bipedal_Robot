import csv, os

path = r"D:\Github\Bipedal_Robot\SADb_audit\curation_log.csv"
rows = [
    ["2026-09-16", "T5", "DONE", "PROCESS", "TASK 5 COMPLETE: 360 Airtable records created via REST API (10/batch, 36 batches) from auto-classified gap papers; 12 created earlier = 372 total (2 dupes skipped). 230 with auto-generated notes (grounded in own text), 130 without (insufficient extraction quality - flagged for manual curation batches). 35 Review Papers records + 115 Models records auto-created and linked. Airtable Papers table now ~930 records. 266/374 had extractable text; 108 no PDF (flagged INSUFFICIENT). Auto-classification: keyword-based (review/model/animal), imperfect but labeled per paper for Ben to correct. All records have Name/PrimaryAuthor/Year/DOI/SecondaryAuthors", "", "", "", "Notes quality is mixed - old scanned PDFs extract garbled text; flagged for batch-curation like batches 1-5. Ben must review and correct classifications."],
    ["2026-09-16", "T5", "ZOTERO", "PENDING", "Ben said he will handle Zotero uploads himself for the 374; personal Zotero is FULL; AARL group PDFs already exist for most; staged copies retained at D:\\sadb_pdf_staging as source (2213 files remain after 173 consumed by tunnel attach)", "", "", "", ""],
]
with open(path, "a", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    for row in rows:
        w.writerow(row)
print("log written")
