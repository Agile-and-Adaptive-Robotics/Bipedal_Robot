"""Dump current bilateralrg template edges (from,to,tag,sign) + node types."""
import json

TPL = json.load(open(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json"))
tpl = TPL["bilateralrg"]
nodes = {n["label"]: n["type"] for n in tpl["nodes"]}
print("== template node types ==")
from collections import Counter
for k, v in Counter(nodes.values()).most_common():
    print(f"   {v:4d}  {k}")
print("\n== template edges ==")
for e in sorted(tpl["edges"], key=lambda x: (x["tag"], x["from"], x["to"])):
    print(f"   {e['from']:<24} -> {e['to']:<24} {e['tag']:<20} {e['sign']}")
