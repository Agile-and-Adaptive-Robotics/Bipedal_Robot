import csv
import json

ROOT = r"D:\Github\Bipedal_Robot\SADb_audit"

d = json.load(open(ROOT + r"\deletion_manifest.json", encoding="utf-8"))
print("papers_total:", d.get("papers_total"))
pgc = d.get("paper_group_counts")
print("paper_group_counts:", pgc if not isinstance(pgc, dict) else {k: v for k, v in list(pgc.items())[:8]})
sat = d.get("satellite")
print("satellite:", sat if not isinstance(sat, dict) else {k: (len(v) if isinstance(v, list) else v) for k, v in sat.items()})
lte = d.get("links_that_will_empty")
print("links_that_will_empty:", lte if not isinstance(lte, dict) else {k: (len(v) if isinstance(v, list) else v) for k, v in lte.items()})

rows = list(csv.DictReader(open(ROOT + r"\deletion_log.csv", encoding="utf-8-sig")))
for r in rows:
    print("log:", r["table"], "planned", r["planned"], "deleted", r["deleted"])
