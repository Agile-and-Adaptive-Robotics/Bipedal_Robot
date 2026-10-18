"""Resolve the 39 'Afferent HC Excite' instances: source neuron, target neuron,
conductance attributes — in document order (for later text-surgery)."""
import json
import re
import xml.etree.ElementTree as ET

PATH = "W2LBRG.aproj"
OLD_TYPE_GUID = "cafe0043-0000-4000-8000-900000000043"

root = ET.parse(PATH).getroot()

def short(t):
    return t.split("}")[-1] if "}" in t else t

def child_map(el):
    m = {}
    for c in el:
        k = short(c.tag)
        if k not in m:
            m[k] = c
    return m

def label_of(el):
    cm = child_map(el)
    for key in ("Text", "UserText", "Name"):
        node = cm.get(key)
        if node is not None and (node.text or "").strip():
            return (node.text or "").strip()
    return "?"

# GUID -> (tag, label) for all ID-bearing elements (FIRST <ID> child wins)
by_guid = {}
for el in root.iter():
    cm = child_map(el)
    idel = cm.get("ID")
    if idel is not None and (idel.text or "").strip():
        by_guid[(idel.text or "").strip()] = (short(el.tag), label_of(el))

rows = []
for el in root.iter():
    cm = child_map(el)
    tid = cm.get("SynapticTypeID")
    if tid is None or (tid.text or "").strip() != OLD_TYPE_GUID:
        continue
    org = (cm.get("OriginID").text or "").strip() if cm.get("OriginID") is not None else "?"
    dst = (cm.get("DestinationID").text or "").strip() if cm.get("DestinationID") is not None else "?"
    cond = cm.get("SynapticConductance")
    condv = cond.attrib.get("Value") if cond is not None else "?"
    conda = cond.attrib.get("Actual") if cond is not None else "?"
    src = by_guid.get(org, ("?", f"?{org[:8]}"))
    tgt = by_guid.get(dst, ("?", f"?{dst[:8]}"))
    rows.append({
        "origin_guid": org, "dest_guid": dst,
        "src": src[1], "src_tag": src[0],
        "tgt": tgt[1], "tgt_tag": tgt[0],
        "cond_value": condv, "cond_actual": conda,
    })

for i, r in enumerate(rows):
    print(f"{i:2d} {r['src']:<28} -> {r['tgt']:<28} cond={r['cond_value']}u ({r['cond_actual']})")

raw = open(PATH, encoding="utf-8", errors="surrogateescape").read()
occ = [m.start() for m in re.finditer(re.escape(OLD_TYPE_GUID), raw)]
print(f"\nraw-text occurrences of old type GUID: {len(occ)} (expect 40 = 1 type ID + 39 instances)")
with open("aff_instances.json", "w") as f:
    json.dump(rows, f, indent=1)
