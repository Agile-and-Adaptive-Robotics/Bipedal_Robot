"""Dump full structure: one complete Node, the SynapseTypes table, node type census."""
import json
import xml.etree.ElementTree as ET
from collections import Counter

APROJ = r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj"

tree = ET.parse(APROJ)
root = tree.getroot()

# Full census top 60
c = Counter(el.tag for el in root.iter())
print("== FULL census top 60 ==")
for tag, n in c.most_common(60):
    print(f"{n:6d}  {tag}")

# One full Node dump (pretty)
print("\n== FULL first Node ==")
for el in root.iter("Node"):
    print(ET.tostring(el, encoding="unicode")[:3000])
    break

# Node Type census
print("\n== Node TypeName census ==")
tc = Counter()
for el in root.iter("Node"):
    tn = el.find("TypeName")
    if tn is not None:
        tc[tn.text] += 1
    else:
        tc["<" + ">".join(list(ch.tag for ch in el)[:6]) + ">"] += 1
for k, v in tc.most_common():
    print(f"{v:5d}  {k}")

# SynapseTypes table
print("\n== SynapseTypes ==")
for st in root.iter("SynapseTypes"):
    for lnk in st:
        print(json.dumps(lnk.attrib))
        pl = {ch.tag: (ch.text or "").strip() for ch in lnk if len(ch) == 0}
        if pl:
            print("   ", json.dumps(pl))
