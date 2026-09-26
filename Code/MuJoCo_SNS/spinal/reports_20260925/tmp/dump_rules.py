"""Dump the structure of Ben's master circuit-rules connectome JSON."""
import json
from collections import Counter

p = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\ben_rules_20260924.json"
doc = json.load(open(p, encoding="utf-8"))
nodes = doc["nodes"]
edges = doc["edges"]
print(f"nodes={len(nodes)} edges={len(edges)}")

by_type = Counter(n["type"] for n in nodes)
print("\n== node types ==")
for t, c in by_type.most_common():
    print(f"  {t}: {c}")

print("\n== nodes by group (label prefix before ':') ==")
grp = Counter()
for n in nodes:
    lab = n["label"]
    g = lab.split(":")[0] if ":" in lab else lab
    grp[g] += 1
for g, c in sorted(grp.items()):
    print(f"  {g}: {c}")

print("\n== all edges (from -> to | sign | gain | tag) ==")
for e in edges:
    print(f"  {e['from']:35s} -> {e['to']:35s} {e.get('sign','?'):3s} "
          f"g={e.get('gain','?')} tag={e.get('tag','')}")
