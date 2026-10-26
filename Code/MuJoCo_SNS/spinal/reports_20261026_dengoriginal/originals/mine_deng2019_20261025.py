# Mine Deng 2019's actual AnimatLab project (ModelDB 264844) for the
# per-synapse wiring + conductances. Reuses the bilateralrg2 miner's
# aproj-parsing pattern (spinal/w2l_cpg/make_bilateralrg2.py), extended to
# record PER-INSTANCE synapse fields per Ben's 2026-10-24 per-target
# ruling (and to dump the Link child-tag inventory so we know what
# per-instance fields AnimatLab actually stores).
import io
import json
import sys
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

DEV = Path(r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev")
APROJ = DEV / "deng2019_modeldb" / "K_A wirh RG  march 21" / "K_A synergy.aproj"
OUT = DEV / "deng2019_mine_20261025.json"


def _val(el, tag):
    ch = el.find(tag)
    if ch is None:
        return None
    for k in ("Actual", "Value"):
        if k in ch.attrib:
            return float(ch.attrib[k])
    if ch.text:
        try:
            return float(ch.text)
        except ValueError:
            return ch.text
    return None


def _txt(el, tag, default=None):
    ch = el.find(tag)
    return (ch.text if ch is not None and ch.text else default)


root = ET.parse(APROJ).getroot()

# synapse TYPES (class-level gmax/Eref)
syn_types = {}
for st in root.iter("SynapseTypes"):
    for lnk in st:
        syn_types[_txt(lnk, "ID")] = {
            "name": _txt(lnk, "Name"),
            "class": (_txt(lnk, "ClassName") or "").split(".")[-1],
            "eref": _val(lnk, "EquilibriumPotential"),
            "gmax": _val(lnk, "MaxSynapticConductance"),
        }

nodes = {}
for node in root.iter("Node"):
    cls = (_txt(node, "ClassName") or "")
    if "NeuralModule" in cls:
        continue
    gid = _txt(node, "ID")
    nodes[gid] = {
        "guid": gid,
        "name": _txt(node, "Text") or _txt(node, "Name"),
        "class": cls.split(".")[-1],
        "tonic": _val(node, "TonicStimulus"),
        "vrest": _val(node, "RestingPotential"),
        "g leakage": _val(node, "LeakageConductance"),
        "cm": _val(node, "MembraneCapacitance"),
        "gna": _val(node, "FastSodiumConductance") or _val(node, "PersistentSodiumConductance"),
        "tau_h_max": _val(node, "HASwoNaChannelTimeConstant") or _val(node, "TimeConstant"),
    }

link_tags = Counter()
edges = []
for lk in root.iter("Link"):
    if not (_txt(lk, "ClassName") or "").endswith("Synapse"):
        continue
    for ch in lk:
        link_tags[ch.tag] += 1
    t = syn_types.get(_txt(lk, "SynapticTypeID"), {})
    edges.append({
        "src": nodes.get(_txt(lk, "OriginID"), {}).get("name"),
        "dst": nodes.get(_txt(lk, "DestinationID"), {}).get("name"),
        "type": t.get("name"),
        "eref": t.get("eref"),
        "type_gmax": t.get("gmax"),
        "instance_amp": _val(lk, "SynAmp") or _val(lk, "SynapticConductance"),
        "delay": _val(lk, "Delay"),
        "src_guid": _txt(lk, "OriginID"),
        "dst_guid": _txt(lk, "DestinationID"),
    })

print("APROJ:", APROJ.name)
print("nodes:", len(nodes), " synapse edges:", len(edges), " synapse types:", len(syn_types))
print("Link child-tag inventory:", dict(link_tags))
types_used = Counter(e["type"] for e in edges)
print("types used:", dict(types_used))
# per-type conductance table
for tid, t in syn_types.items():
    if t["name"] in types_used:
        print("  TYPE %-28s class=%-22s gmax=%s eref=%s"
              % (t["name"], t["class"], t["gmax"], t["eref"]))
# per-instance deviation check
inst = [e for e in edges if e["instance_amp"] is not None]
print("edges carrying a per-instance amplitude field:", len(inst))

out = {"_source": str(APROJ), "_note": "Deng 2019 ModelDB 264844 K_A synergy project, mined 2026-10-25 with per-instance fields recorded per Ben's per-target ruling.",
       "nodes": list(nodes.values()), "edges": edges,
       "synapse_types": list(syn_types.values())}
OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")
print("wrote", OUT)

# quick neuron census
census = Counter(n["class"] for n in nodes.values())
print("neuron census:", dict(census))
