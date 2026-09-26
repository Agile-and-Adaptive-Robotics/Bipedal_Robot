# -*- coding: utf-8 -*-
"""Check Box bodies for LWH children; inventory femur_L element."""
import xml.etree.ElementTree as ET

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

def txt(el, tag):
    c = el.find(tag)
    return (c.text or "").strip() if c is not None else None

org = root.find(".//Organism")

boxes = []
def find_boxes(el):
    for rb in el.findall("RigidBody"):
        if txt(rb, "Type") == "Box":
            boxes.append(rb)
        cb = rb.find("ChildBodies")
        if cb is not None:
            find_boxes(cb)
find_boxes(org)

print(f"{len(boxes)} Box bodies found:")
for rb in boxes:
    nm = txt(rb, "Name")
    lwh = {}
    for ax in ("Length", "Width", "Height"):
        c = rb.find(ax)
        lwh[ax] = (c.get("Value"), c.get("Scale"), c.get("Actual")) if c is not None else None
    print(f"  {nm}: L={lwh['Length']} W={lwh['Width']} H={lwh['Height']}")

# inventory femur_L
print("\nfemur_L direct children (tags in order):")
for rb in org.iter("RigidBody"):
    if txt(rb, "Name") == "femur_L" and txt(rb, "ID") == "20eb335a-ab76-41dd-be7d-885c53010b81":
        print("  ", [c.tag for c in rb])
        break

# also check GroundPlane/WalkingPath types
print("\nNon-organism structures:")
env = root.find(".//Environment")
st = env.find("Structures")
if st is not None:
    for s in st:
        print("  ", s.tag, txt(s, "Name"), txt(s, "Type") or txt(s, "PartType"))
