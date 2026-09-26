# -*- coding: utf-8 -*-
"""Inspect CollisionExclusionPairs in the aproj."""
import xml.etree.ElementTree as ET
P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()
org = root.find(".//Organism")
cep = org.find("CollisionExclusionPairs")
if cep is None or len(list(cep)) == 0:
    print("CollisionExclusionPairs: EMPTY or absent -> all collisions active (children tags:",
          [c.tag for c in cep] if cep is not None else "absent", ")")
else:
    ids = {txt.text.strip() for txt in cep.iter() if txt.tag in ("ID", "BodyID") and txt.text}
    for c in cep:
        print(c.tag, {k: v for k, v in c.attrib.items()}, (c.text or "").strip()[:60])
    # resolve GUIDs to names
    def name_of(gid):
        for rb in org.iter("RigidBody"):
            i = rb.find("ID")
            if i is not None and i.text.strip() == gid:
                return rb.find("Name").text.strip()
        return gid
    for c in cep.iter("CollisionExclusionPair"):
        pass
    for gid_el in cep.iter():
        if gid_el.tag in ("ID",) and gid_el.text and len(gid_el.text) == 36:
            print("  ", name_of(gid_el.text.strip()))
