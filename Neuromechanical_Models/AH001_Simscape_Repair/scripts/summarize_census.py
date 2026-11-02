import json
import os
import collections

REPAIR = r"C:\Users\Ben\Documents\GitHub\Bipedal_Robot\Neuromechanical_Models\AH001_Simscape_Repair"
d = json.load(open(os.path.join(REPAIR, "logs", "sw_census_20261101.json")))

roots = collections.Counter()
exts = collections.Counter()
for c in d["components"]:
    p = c["path"].replace("\\", "/")
    idx = p.find("/Solid_Models/")
    sub = p[idx:] if idx >= 0 else p
    parts = sub.split("/")
    roots["/".join(parts[:4])] += 1
    exts[os.path.splitext(p)[1].lower()] += 1
for k, v in sorted(roots.items()):
    print(v, k)
print("extensions:", dict(exts))
print("registry mentions:")
for line in d.get("registry_link_mentions", []):
    print("  ", line)
