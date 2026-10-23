# Proof-3 analysis: diff the two scratch shortlists; check real state isolation.
import json

def dois(path):
    rows = json.load(open(path, encoding="utf-8"))
    return {r["doi"]: r for r in rows}, rows

base, base_rows = dois(r"D:\Github\Bipedal_Robot\tmp\watch_base\reports\2026-10_shortlist.json")
mod, mod_rows = dois(r"D:\Github\Bipedal_Robot\tmp\watch_mod\reports\2026-10_shortlist.json")

print("baseline rows:", len(base_rows), "| modified rows:", len(mod_rows))
print("baseline subset of modified:", set(base) <= set(mod))
print("rows in baseline but NOT in modified:", sorted(set(base) - set(mod)))
added = sorted(set(mod) - set(base))
print("rows ADDED by the config change:", len(added))
for d in added:
    r = mod[d]
    print("  + %-30s | trigger=%-16s | why=%s" % (d, r["trigger"], r["why"][:80]))

st = json.load(open(r"D:\Github\Bipedal_Robot\SADb_audit\watch\state.json", encoding="utf-8"))
print("\nreal state: routed DOIs:", len(st["routed_dois"]),
      "| run_log entries:", len(st["run_log"]),
      "| months:", sorted(st["months"]),
      "| stored 2026-10 entries:", len(st["months"]["2026-10"]["entries"]),
      "| stored window:", st["months"]["2026-10"]["window"])
