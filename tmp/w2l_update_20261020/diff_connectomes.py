"""Diff the 2026-10-18 aproj mine vs the 2026-10-03 canonical export vs the
bilateralrg template (what the runnable build uses)."""
import json

M = json.load(open(r"D:\GitHub\Bipedal_Robot\tmp\w2l_update_20261020\w2l_aproj_mine_20261018.json"))
CAN = json.load(open(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\ben_connectome_canonical.json"))
TPL = json.load(open(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json"))
tpl = TPL["bilateralrg"]

nodes = M["nodes"]
print("== aproj node names (96 NonSpiking neurons) ==")
ns = sorted((n["name"] or "<unnamed>") for n in nodes if n["class"] == "NonSpiking")
for n in ns:
    print("  ", n)

print("\n== aproj adapters (PhysicalToNode) ==")
for n in sorted((n for n in nodes if n["class"] == "PhysicalToNodeAdapter"), key=lambda x: x["name"] or ""):
    print(f"   {n.get('name'):<28} src_data={n.get('SourceDataType')} origin={n.get('origin','')}")

print("\n== aproj named edges (all 204) ==")
for e in sorted(M["edges"], key=lambda x: (x["type"] or "", x["dst"] or "")):
    src = str(e["src"]); dst = str(e["dst"])
    print(f"   {src:<26} -> {dst:<26} {e['type']}  g={e['g']}")

# canonical export edge set (by from/to labels)
can_edges = {}
for e in CAN["edges"]:
    k = (e["from"], e["to"])
    can_edges.setdefault(k, []).append(e)

print("\n== template bilateralrg census ==")
print("  nodes:", len(tpl.get("nodes", [])), " edges:", len(tpl.get("edges", [])))
tc = {}
for e in tpl.get("edges", []):
    tc[e.get("tag", "?")] = tc.get(e.get("tag", "?"), 0) + 1
for k, v in sorted(tc.items()):
    print(f"   {v:4d}  {k}")
