import csv
import json

d = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
print("export json type:", type(d).__name__, "len:", len(d))
if isinstance(d, list) and d:
    print("keys of first:", sorted(d[0].keys())[:20])

rows = list(csv.DictReader(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.csv", encoding="utf-8-sig")))
print("export csv rows:", len(rows))
if rows:
    print("csv columns:", sorted(rows[0].keys())[:25])

# the built app payload count
html = open(r"D:\Github\Bipedal_Robot\SADb_audit\app\sadb_app.html", encoding="utf-8", errors="replace").read()
print("app html size MB:", round(len(html) / 1e6, 1))
for marker in ('"papers"', "SADB_DATA", "EMBED", "__DATA__", "const DATA"):
    i = html.find(marker)
    if i >= 0:
        print("marker", marker, "at", i)
