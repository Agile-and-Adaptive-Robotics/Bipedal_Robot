"""Full miner: nodes + synapse instances + adapters from the 2026-10-18 aproj.
Outputs w2l_aproj_mine_20261018.json with named edge list."""
import json
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict

APROJ = r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj"
OUT = r"D:\GitHub\Bipedal_Robot\tmp\w2l_update_20261020\w2l_aproj_mine_20261018.json"


def txt(el, tag, default=None):
    ch = el.find(tag)
    if ch is None:
        return default
    return ch.text if ch.text is not None else default


def val(el, tag):
    """AnimatLab param elements carry Value/Actual attrs; prefer Actual (SI)."""
    ch = el.find(tag)
    if ch is None:
        return None
    if "Actual" in ch.attrib:
        return float(ch.attrib["Actual"])
    if "Value" in ch.attrib:
        return float(ch.attrib["Value"])
    return ch.text


tree = ET.parse(APROJ)
root = tree.getroot()

# --- synapse TYPE table (per neural module; collect all) ---
syn_types = {}  # guid -> {name, class, erev, gmax, presyn_thr, presyn_sat}
for st in root.iter("SynapseTypes"):
    for lnk in st:
        gid = txt(lnk, "ID")
        syn_types[gid] = {
            "name": txt(lnk, "Name"),
            "class": (txt(lnk, "ClassName") or "").split(".")[-1],
            "eref": val(lnk, "EquilibriumPotential"),
            "gmax": val(lnk, "MaxSynapticConductance"),
            "syn_g": val(lnk, "SynapticConductance"),
            "pre_thr": val(lnk, "PreSynapticThreshold"),
            "pre_sat": val(lnk, "PreSynapticSaturationLevel"),
        }

# --- nodes ---
nodes = {}       # guid -> dict
edges = []       # {src, dst, type_guid, type_name, sign, g, delay}
adapter_links = []

for node in root.iter("Node"):
    cls = (txt(node, "ClassName") or "")
    if "NeuralModule" in cls:
        continue
    gid = txt(node, "ID")
    name = txt(node, "Text") or txt(node, "Name")
    short = cls.split(".")[-1]
    d = {"guid": gid, "name": name, "class": short}
    for p in ("FiringRateVMidpoint", "FiringRateVSlope", "TimeConstant",
              "TonicStimulus", "InitialVmem", "Adaptation", "AdaptationPM",
              "RelativeAccommodation", "AccommodationTimeConstant",
              "Gain", "SourceDataType", "Enabled", "FieldName", "DataTypeID",
              "TargetDataType", "StimulusOnDuration", "StimulusOffDuration",
              "StimulusStartTime", "StimulusIncrement"):
        v = val(node, p)
        if v is not None:
            d[p] = v
    # body binding of adapters lives in <OriginID> child (a body GUID)
    o = txt(node, "OriginID")
    if o:
        d["origin"] = o
    nodes[gid] = d

# --- synapse instances: Link elements with ClassName *Synapse, under any Links ---
syn_links = []
for lk in root.iter("Link"):
    cls = (txt(lk, "ClassName") or "")
    if not cls.endswith("Synapse"):
        continue
    stid = txt(lk, "SynapticTypeID")
    t = syn_types.get(stid, {})
    syn_links.append({
        "src_guid": txt(lk, "OriginID"),
        "dst_guid": txt(lk, "DestinationID"),
        "type": t.get("name"),
        "class": t.get("class"),
        "eref": t.get("eref"),
        "type_gmax": t.get("gmax"),
        "type_syn_g": t.get("syn_g"),
        "g": val(lk, "SynapticConductance"),
        "delay": val(lk, "ConductionDelay"),
    })

# --- adapter Links (data-flow: body/organ -> adapter node) ---
for lk in root.iter("Link"):
    cls = (txt(lk, "ClassName") or "")
    if "Adapter" in cls or "Gain" in cls or txt(lk, "DestinationID") or txt(lk, "OriginID"):
        # only keep data links under Links sections (not SynapseTypes/drawings)
        pass

# body name map (for adapter origins)
body_names = {}
for b in root.iter("RigidBody"):
    bid = txt(b, "ID")
    nm = txt(b, "Name")
    if bid and nm:
        body_names[bid] = nm

# resolve source/dest names
for e in syn_links:
    e["src"] = nodes.get(e["src_guid"], {}).get("name", body_names.get(e["src_guid"], e["src_guid"]))
    e["dst"] = nodes.get(e["dst_guid"], {}).get("name", body_names.get(e["dst_guid"], e["dst_guid"]))
edges = syn_links

out = {
    "source": APROJ,
    "syn_types": {v["name"]: v for v in syn_types.values()},
    "nodes": list(nodes.values()),
    "bodies": body_names,
    "edges": edges,
}
with open(OUT, "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)

print(f"nodes={len(nodes)} synapse_types={len(syn_types)} edges={len(edges)} bodies={len(body_names)}")
print("\n== node class census ==")
for k, v in Counter(n["class"] for n in nodes.values()).most_common():
    print(f"{v:5d}  {k}")
print("\n== edge type census ==")
for k, v in Counter(e["type"] for e in edges).most_common():
    print(f"{v:5d}  {k}")
print("\n== contact-related nodes ==")
for n in nodes.values():
    if n["class"] == "Adapter" or "contact" in (n["name"] or "").lower():
        print(json.dumps(n))
print("\n== cafe-type edges (the new build-script wiring) ==")
for e in edges:
    if e["type"] and ("Contra" in e["type"] or "AnkMNext" in e["type"] or "Ia_HipFlx" in e["type"]):
        print(f"{e['src']:>28} -> {e['dst']:<28} {e['type']}  g={e['g']}")
