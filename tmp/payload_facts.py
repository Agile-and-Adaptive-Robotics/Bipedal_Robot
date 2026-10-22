import json

d = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
arch = [r for r in d if r.get("archived") is True]
keep = [r for r in d if r.get("archived") is not True]
print("total paper records:", len(d))
print("retained (not archived):", len(keep), " archived (pruned-from-Airtable):", len(arch))


def cnt(rows, f):
    return sum(1 for r in rows if r.get(f))


print()
print("WITH curation notes:  all", cnt(d, "has_notes"), "| retained", cnt(keep, "has_notes"), "| archived", cnt(arch, "has_notes"))
print("WITH pdf attachment flag: all", cnt(d, "has_pdf"), "| retained", cnt(keep, "has_pdf"), "| archived", cnt(arch, "has_pdf"))
print("WITH doi: all", sum(1 for r in d if r.get("doi")), "| retained", sum(1 for r in keep if r.get("doi")), "| archived", sum(1 for r in arch if r.get("doi")))
print("WITH authors (primary): all", sum(1 for r in d if r.get("primary")), "| archived", sum(1 for r in arch if r.get("primary")))
print("WITH year: all", sum(1 for r in d if r.get("year")), "| archived", sum(1 for r in arch if r.get("year")))

# what a pruned record actually looks like (first archived one, trimmed)
import copy

ex = copy.deepcopy(arch[0])
for k, v in list(ex.items()):
    if isinstance(v, str) and len(v) > 60:
        ex[k] = v[:57] + "..."
    if isinstance(v, list) and len(v) > 3:
        ex[k] = v[:3] + ["..."]
print()
print("example archived record fields:")
for k in sorted(ex):
    print(" ", k, "=", ex[k])
