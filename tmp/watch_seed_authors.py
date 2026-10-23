# Seed derivation for SADb_audit/watch/watch_config.json (one-off, 2026-10-22).
# Top ~20 primary authors by RETAINED-paper count from the export json.
# Papers only (recordType == 'paper'); archived papers NOT counted ("retained-paper count").
import json, collections, os

d = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\export\sadb_export.json", encoding="utf-8"))
papers = [r for r in d if r.get("recordType") == "paper"]
retained = [r for r in papers if not r.get("archived")]
print("papers", len(papers), "| retained", len(retained))

cnt = collections.Counter()
for r in retained:
    a = (r.get("primary") or "").strip()
    if a:
        cnt[a] += 1

print("\nTop 30 primary authors by retained-paper count:")
for a, n in cnt.most_common(30):
    print("%3d  %s" % (n, a))
