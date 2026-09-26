import json, io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
T = json.load(open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json", encoding="utf-8"))
w = T["w2laproj"]
print("KEYS:", list(w.keys()))
print("NOTE:", w.get("_note", "")[:400])
print("n_nodes", len(w["nodes"]), "n_edges", len(w["edges"]))
print("TYPES:", sorted(set(n["type"] for n in w["nodes"])))
out = []
for n in w["nodes"]:
    out.append("N %-4s | %-10s | grp=%s | extra=%s" % (n["label"], n["type"], n.get("grp", ""), {k: v for k, v in n.items() if k not in ("label", "type", "grp", "x", "y")}))
out.append("")
from collections import defaultdict
for e in w["edges"]:
    out.append("E %-30s -> %-30s | g=%s | sign=%s | tag=%s | extra=%s" % (
        e["from"], e["to"], e.get("g", ""), e.get("sign", ""), e.get("tag", ""),
        {k: v for k, v in e.items() if k not in ("from", "to", "g", "sign", "tag", "pts")}))


# summarize edge tags
tags = defaultdict(int)
for e in w["edges"]:
    tags[e.get("tag", "?")] += 1
out.append("")
out.append("TAG COUNTS: " + json.dumps(dict(sorted(tags.items())), indent=None))
open(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal\reports_20260925\tmp\w2laproj_dump.txt", "w", encoding="utf-8").write("\n".join(out))
print("wrote dump, lines:", len(out))
