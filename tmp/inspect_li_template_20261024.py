# Inspect the `li` template: node/edge schema, per-synapse gains present?
import json

p = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json"
d = json.load(open(p, encoding="utf-8"))
li = d["li"]
print("top keys:", sorted(li.keys()))
nodes = li["nodes"]
edges = li["edges"]
print("nodes:", len(nodes), "edges:", len(edges))
print("node example:", json.dumps(nodes[0], indent=1)[:300])
print("edge example:", json.dumps(edges[0], indent=1)[:400])
# what fields do edges carry?
keys = set()
for e in edges:
    keys.update(e.keys())
print("edge fields union:", sorted(keys))
# gain distribution if present
gains = [e.get("gain") for e in edges if e.get("gain") is not None]
print("edges with gain:", len(gains), "unique:", sorted(set(gains))[:20])
types = {}
for e in edges:
    t = e.get("type", e.get("kind", "?"))
    types[t] = types.get(t, 0) + 1
print("edge types:", types)
