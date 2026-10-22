import csv
import json
import os

d = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
outdir = r"D:\Github\Bipedal_Robot\Documentation\Reports and Papers\Dissertation\Notes\advisor_revision_20261022"
rows = []
for r in d:
    rows.append(
        {
            "doi": (r.get("doi") or "").strip().lower().replace("https://doi.org/", ""),
            "title": r.get("title") or "",
            "primary_author": r.get("primary") or "",
            "year": r.get("year") or "",
            "archived": "yes" if r.get("archived") is True else "no",
        }
    )
rows.sort(key=lambda x: x["doi"])

csv_path = os.path.join(outdir, "corpus_dois_943.csv")
with open(csv_path, "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["doi", "title", "primary_author", "year", "archived"])
    w.writeheader()
    w.writerows(rows)

with_doi = sum(1 for r in rows if r["doi"])
txt_path = os.path.join(outdir, "corpus_dois_943.txt")
with open(txt_path, "w", encoding="utf-8", newline="\n") as f:
    f.write("\n".join(r["doi"] for r in rows if r["doi"]))

print("records:", len(rows), "with DOI:", with_doi, "(no DOI:", len(rows) - with_doi, ")")
print("archived:", sum(1 for r in rows if r["archived"] == "yes"))
print("wrote:", csv_path)
print("wrote:", txt_path)
