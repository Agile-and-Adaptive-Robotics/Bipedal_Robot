# Diff the Audet-supplementary scratch shortlist against the seeded baseline.
import json

def dois(path):
    rows = json.load(open(path, encoding="utf-8"))
    return {r["doi"]: r for r in rows}

base = dois(r"D:\Github\Bipedal_Robot\tmp\watch_base\reports\2026-10_shortlist.json")
aud = dois(r"D:\Github\Bipedal_Robot\tmp\watch_audet\reports\2026-10_shortlist.json")
print("baseline rows:", len(base), "| audet-config rows:", len(aud))
print("baseline subset:", set(base) <= set(aud))
print("ADDED rows:")
for d in sorted(set(aud) - set(base)):
    r = aud[d]
    print("  + %-32s | trigger=%-16s | why=%s" % (d, r["trigger"], r["why"][:90]))
