# Probe SADb export distributions to design the Ch5 pivot + bubble figures.
import csv, json, sys
from collections import Counter, defaultdict

EXPORT = r"D:\Github\Bipedal_Robot\SADb_audit\export"
rows = list(csv.DictReader(open(EXPORT + r"\sadb_export.csv", encoding="utf-8")))
print("total records:", len(rows))

retained = [r for r in rows if r["prune"] != "True" and r["archived"] != "True"]
archived = [r for r in rows if r["prune"] == "True" or r["archived"] == "True"]
print("retained:", len(retained), " archived:", len(archived))

def split(s):
    return [t.strip() for t in s.replace(";", ",").split(",") if t.strip()]

# what do the multi-value columns look like?
for col in ("animals", "afferents", "feedback"):
    vals = Counter()
    nonempty = 0
    for r in retained:
        v = r[col].strip()
        if v:
            nonempty += 1
            for t in split(v):
                vals[t] += 1
    print(f"\n{col}: nonempty retained = {nonempty}")
    for k, c in vals.most_common(15):
        print(f"   {k}: {c}")

# year coverage
years = [int(r["year"]) for r in retained if r["year"].strip().isdigit()]
print("\nretained year range:", min(years), "-", max(years), " n with year:", len(years))

# clusters
cl = json.load(open(EXPORT + r"\cluster_labels.json", encoding="utf-8"))
print("\ncluster_labels type:", type(cl).__name__, "len:", len(cl))
if isinstance(cl, dict):
    items = list(cl.items())[:3]
else:
    items = cl[:3]
print("sample:", json.dumps(items, ensure_ascii=False)[:400])

cp = json.load(open(EXPORT + r"\cluster_profiles.json", encoding="utf-8"))
print("\ncluster_profiles type:", type(cp).__name__, "len:", len(cp))
s = json.dumps(cp if not isinstance(cp, dict) else list(cp.items())[:2], ensure_ascii=False)
print("sample:", s[:600])

lay = json.load(open(EXPORT + r"\sadb_layout.json", encoding="utf-8"))
print("\nlayout type:", type(lay).__name__, "len:", len(lay))
