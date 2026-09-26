# -*- coding: utf-8 -*-
"""Survey aproj physical structure: RigidBody tree, joints, muscles."""
import xml.etree.ElementTree as ET
from collections import Counter

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

env = root.find(".//Environment")
mu = env.find("MassUnits"); du = env.find("DistanceUnits")
print("MassUnits:", dict(mu.attrib) if mu is not None else None)
print("DistanceUnits:", dict(du.attrib) if du is not None else None)

# Organism / Structures
orgs = root.find(".//Organisms")
print("\nOrganisms children:")
for o in orgs:
    print("  ", o.tag, dict(o.attrib))
    for c in o:
        if c.tag in ("Name", "ID"):
            print("     ", c.tag, c.text)

cnt = Counter()
for el in root.iter():
    cnt[el.tag] += 1
print("\nTag counts of interest:")
for t in ("RigidBody", "Attachment", "Muscle", "Joint", "Hinge", "LinearHillMuscle",
          "Structure", "Organism", "BodyPart", "Box", "MaterialType", "Material",
          "Spring", "StretchReceptor", "Biped"):
    print(f"  {t}: {cnt.get(t, 0)}")
