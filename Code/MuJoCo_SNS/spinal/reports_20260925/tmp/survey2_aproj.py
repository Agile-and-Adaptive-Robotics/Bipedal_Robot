# -*- coding: utf-8 -*-
"""Survey aproj: find physical structure classes."""
import xml.etree.ElementTree as ET

P = r"D:\Github\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG\Walker_2_Layer_CPG.aproj"
tree = ET.parse(P)
root = tree.getroot()

# Find all ClassName values
from collections import Counter
cls = Counter()
for el in root.iter("ClassName"):
    cls[(el.text or "").strip()] += 1
print("CLASSNAME counts:")
for k, v in cls.most_common():
    print(f"  {k}: {v}")

# Environment: units
env = root.find(".//Environment")
print("\nEnvironment direct children:")
for c in env:
    print("  ", c.tag, dict(c.attrib))
