# -*- coding: utf-8 -*-
"""Survey aproj: dump RigidBody names/classes, Joint entries."""
import xml.etree.ElementTree as ET
from collections import Counter

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

def txt(el, tag):
    c = el.find(tag)
    return (c.text or "").strip() if c is not None else None

print("RigidBody entries (name | class | id | parent):")
bodies = root.iter("RigidBody")
rows = []
for rb in root.iter("RigidBody"):
    name = txt(rb, "Name"); cid = txt(rb, "ID")
    cn = txt(rb, "ClassName")
    rows.append((name, cn, cid))
for name, cn, cid in rows:
    print(f"  {name} | {cn} | {cid}")

print("\nJoints:")
for jt in root.iter("Joint"):
    print("  ", txt(jt, "Name"), "|", txt(jt, "ClassName"), "|", txt(jt, "ID"))
    for c in jt:
        if c.tag not in ("Name", "ID", "ClassName", "Enabled"):
            pass

# Muscle entries in behavior nodes
print("\nBehavior Muscle nodes:")
for m in root.iter("Muscle"):
    print("  ", txt(m, "Name"), "|", txt(m, "ID"))

# Structure children under Organism
org = root.find(".//Organism")
print("\nOrganism direct children tags:", [c.tag for c in org])
body = org.find("BodyPart") or org
