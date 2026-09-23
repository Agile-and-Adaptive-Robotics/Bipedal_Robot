"""W2L BilateralRG .aproj -> block-editor JSON extractor.

The aproj stores the behavior diagram as an AddFlow graph: <Link
Org="i" Dst="j"> are INDICES into the diagram's node list, and the
behavior objects (NonSpiking neurons, Synapses, Muscles,
StretchReceptors) carry <ID> GUIDs that the graphical <Node>/<Link>
drawing elements reference. This extractor:
  1. pulls the AddFlow node array (in document order) and the
     behavior objects (ID, Name, ClassName),
  2. maps Link Org/Dst indices -> behavior names,
  3. classifies each synapse's sign from its chemical type/equilibrium,
  4. emits connectome_block.json (editor Load format).
Validation: every Org/Dst must land inside the node array.
"""
import json
import re
import sys

APROJ = (r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
         r'\Walker_2_Layer_CPG_BilateralRG'
         r'\Walker_2_Layer_CPG_BilateralRG.aproj')
xml = open(APROJ, encoding="utf-8", errors="ignore").read()

# --- behavior objects: <Node> blocks whose child <ClassName> is a
# behavior class; they carry <ID> and <Name>. Non-greedy block scan.
blocks = re.findall(
    r'<ClassName>([^<]+)</ClassName>.{0,4000}?<ID>([^<]+)</ID>'
    r'.{0,4000}?<Name>([^<]+)</Name>', xml, re.S)
objs = {}
for cls, oid, name in blocks:
    if oid not in objs:
        objs[oid] = {"cls": cls, "name": name}
print(f"behavior objects: {len(objs)}")

# --- synapse type details: find each Synapse object's type + equil
# (SynapseTypes.NonSpikingChemical blocks carry <Equilibrium>).
# Map ID -> ('inh'|'exc'|'elec')
signs = {}
for m in re.finditer(
        r'<ClassName>(IntegrateFireGUI\.DataObjects\.Behavior\.'
        r'SynapseTypes\.(\w+))</ClassName>.{0,3000}?<ID>([^<]+)</ID>'
        r'.{0,3000}?(?:<Equilibrium>([-\d.e]+)</Equilibrium>)?', xml,
        re.S):
    cls, kind, sid, equil = m.group(1), m.group(2), m.group(3), m.group(4)
    if kind == "NonSpikingChemical":
        try:
            signs[sid] = ("inh" if equil and float(equil) < 0
                          else "exc")
        except (TypeError, ValueError):
            signs[sid] = "exc"
    elif kind == "SpikingChemical":
        signs[sid] = "exc"     # spiking chemical: sign set by type field
    else:
        signs[sid] = "elec"
print(f"synapse signs resolved: {len(signs)}")

# --- AddFlow diagram: nodes in document order (the <Node Left=...>
# graphics entries) each reference a behavior object by ID in their
# child <ID>; links Org/Dst are 1-based? indices into that order.
af_nodes = re.findall(r'<Node (Left="[^"]*" Top="[^"]*"[^>]*)>'
                      r'(.{0,2000}?)</Node>', xml, re.S)
order = []          # list of (behavior_id | None)
for attrs, body in af_nodes:
    mid = re.search(r'<ID>([^<]+)</ID>', body)
    order.append(mid.group(1) if mid else None)
print(f"AddFlow graphic nodes: {len(order)}; "
      f"with behavior ID: {sum(1 for o in order if o)}")

links = re.findall(r'<Link Org="(\d+)" Dst="(\d+)"', xml)
print(f"AddFlow links: {len(links)}")
maxi = max((max(int(a), int(b)) for a, b in links), default=0)
print(f"max link index {maxi} vs node array {len(order)}")

nodes_out = []
byidx = []
for i, bid in enumerate(order):
    if bid and bid in objs:
        o = objs[bid]
        label = o["name"]
        if "Neurons.NonSpiking" in o["cls"]:
            t = "HC-RG-E" if "_E" in label.upper() else \
                "HC-RG-F" if "_F" in label.upper() else "IN-V0V"
        elif "Synapse" in o["cls"]:
            continue                     # synapses are edges, not nodes
        elif "Muscle" in o["cls"]:
            t = "MUSCLE"
        elif "StretchReceptor" in o["cls"]:
            t = "SN-Ia"
        else:
            t = "IN-V2a"
        byidx.append((i, len(nodes_out)))
        nodes_out.append({"id": f"n{i}", "type": t, "label": label})

idx2nid = {i: f"n{i}" for i, _ in byidx}
syn_ids = set(signs)
syn_out = []
for a, b in links:
    ia, ib = int(a), int(b)
    if ia < len(order) and ib < len(order):
        ba, bb = order[ia], order[ib]
        sa = objs.get(ba, {}).get("cls", "")
        sb = objs.get(bb, {}).get("cls", "")
        if "Synapse" in sa:
            # edge: presynaptic object -> the synapse object
            if ba in signs and bb in objs:
                syn_out.append({"from": objs[ba]["name"],
                                "to": objs[bb]["name"],
                                "sign": signs[ba], "gain": 0.5,
                                "tag": "w2l"})
        elif "Synapse" in sb:
            if bb in signs and ba in objs:
                syn_out.append({"from": objs[ba]["name"],
                                "to": objs[bb]["name"],
                                "sign": signs[bb], "gain": 0.5,
                                "tag": "w2l"})
# dedupe
seen, edges = set(), []
for e in syn_out:
    key = (e["from"], e["to"], e["sign"])
    if key not in seen:
        seen.add(key)
        edges.append(e)

spec = {"nodes": nodes_out, "synapses": edges}
json.dump(spec, open("w2l_bilateralrg_block.json", "w",
                     encoding="utf-8"), indent=1)
print(f"wrote w2l_bilateralrg_block.json: {len(nodes_out)} nodes, "
      f"{len(edges)} synapses")
