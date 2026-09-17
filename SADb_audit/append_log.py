import csv, os

row = [
    "2026-09-16", "MORNING", "PROCESS", "",
    "Ben overnight items processed: (1) paperhunt CSV: Ben marked 26 rows Downloaded (17 x + 8 no + 1 chapter note); all 17 x-rows matched to his 19 staged PDFs via Crossref surname+year (author_fix/pdf_match_pairs.csv) and attached via tunnel bridge 17/17 ok; 2 unmatched files remain: Dissertation_EKBunz.pdf (no hunt row) and Singh-2025-Shox2 JPhysiol PDF (no hunt row) - both flagged for Ben; 8 'no' rows stay on hunt list",
    "yes", "", "",
    "(2) 33 new AARL Zotero refs since 2026-09-15 20:00Z: 3 webpages skip, 11 already in Airtable (digest dupes), 19 NEW imported with full curation - grounded notes from own abstracts, Animals allowed-list only, single-surname Primary Author + Secondary Authors bubbles; Review Papers created: Al-Tameemi 2025 rechuDy6oaqDMCEhv, Koseki 2026 recAK6lmElMGNGzeF; Models created: Ton 2026 recQ590IiNqWml6pA, Schumacher 2025 recXDYv9QBijWlpJI, Mo 2025 recFVP7Gt908vdjm5, Bunz 2025 recANrpGbl2cqCHl9, Steffen 2026 rec8RA461yHxDCdCy, Bekhiti 2025 recOCSKgbkCymcLkt, Molkov 2026 recg6wwqL41yB8DQ2, Molkov 2025 rec7RoxSF7CtSBbcB, Ramalingasetty 2023 recObQF8VVMb7vb0N, Ramalingasetty 2023b recpp3mzBNMWb6ftF (conference; suffix b), Lockhart 2023 recWQpGm5mXdjeJkk, Shevtsova 2025 recvpsfPbh73YSlUx; 18/19 PDFs attached from Zotero storage (Mo 2025 NO-LOCAL-PDF: ScienceDirect link only); Papers table ~570 records",
]
row += [""] * (9 - len(row))
path = r"D:\Github\Bipedal_Robot\SADb_audit\curation_log.csv"
with open(path, "a", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(row)
print("log appended")
