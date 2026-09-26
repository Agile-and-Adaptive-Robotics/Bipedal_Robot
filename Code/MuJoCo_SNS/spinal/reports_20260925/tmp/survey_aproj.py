# -*- coding: utf-8 -*-
"""Survey Walker_2_Layer_CPG.aproj structure (read-only)."""
import xml.etree.ElementTree as ET
from collections import Counter

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()
print("ROOT:", root.tag, dict(root.attrib))

# Environment / units
for child in root:
    print("  <-", child.tag, {k: v for k, v in child.attrib.items() if k in ("Name", "ID")})

# find Environment node with units
def walk(el, path="", depth=0, maxdepth=3):
    if depth > maxdepth:
        return
    for c in el:
        a = {k: v for k, v in c.attrib.items()}
        print("  " * depth + f"<{c.tag}> {a}")
        walk(c, path + "/" + c.tag, depth + 1, maxdepth)

# locate Environment
for el in root.iter("Environment"):
    print("\nENVIRONMENT:", dict(el.attrib))
    break

# Structure: Organism > structures
org = root.find(".//Organism")
if org is not None:
    print("\nORGANISM:", dict(org.attrib))

# count element types in whole file
cnt = Counter()
for el in root.iter():
    cnt[el.tag] += 1
print("\nELEMENT COUNTS (top 40):")
for tag, n in cnt.most_common(40):
    print(f"  {tag}: {n}")
