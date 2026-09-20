"""Append WORKFLOW row: definitive no-PDF list consolidation + full-corpus hunt."""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(os.path.dirname(HERE), "curation_log.csv")

rows = [
    ["2026-09-20", "PDFS", "WORKFLOW", "consolidate",
     "PDF hunt extended to the FULL corpus (the Sept-15 list only covered the 456 records that existed then; 417 task5-era records had never been hunted). 44 more PDFs attached via verified OA URLs (274->318); the 48 'attach failures' were a fields[]-GET 422 bug in my own script, all succeeded on retry with a plain PATCH. THE definitive no-PDF list is now author_fix/remaining_no_pdf.csv alone: columns doi,record_id,title,hunt_status; 581 records (488 no OA copy anywhere, 88 candidate URLs failing %PDF verification = paywalled/bot-blocked, 4 attach stragglers pending Airtable async population, 1 supplementary-PDF false positive correctly rejected). Regenerate with rebuild_no_pdf_list.py. DELETED as superseded (Ben's no-bloat ruling): paperhunt_status.csv, pdfs_missing_no_zotero.csv, pdfs_paywalled.csv, pdfs_to_attach.json, upload_status.csv, fileio_status.csv, tunnel_status.csv, newref_attach_status.csv, pdf_match_pairs.csv, and the dated pass-2/retry files. KEPT: remaining_no_pdf.csv (the hunt list) and pdfs_from_zotero.csv (inventory of where Zotero PDF copies live, for drag-drop). Remaining PDFs need Ben's PSU library access", "", "", "", ""],
]

with open(LOG, "a", encoding="utf-8", newline="") as fh:
    csv.writer(fh).writerows(rows)
print("appended 1 WORKFLOW row")
