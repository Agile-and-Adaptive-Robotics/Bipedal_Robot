"""Dump the element structure of an AnimatLab .aproj: top 3 levels, tag census,
and a sample of SynapseTypes children + a sample synapse/neuron object."""
import sys
import xml.etree.ElementTree as ET
from collections import Counter

root = ET.parse(sys.argv[1]).getroot()
print("root:", root.tag, dict(list(root.attrib.items())[:5]))

def short(tag):
    return tag.split("}")[-1] if "}" in tag else tag

census = Counter(short(el.tag) for el in root.iter())
print("tag census (top 40):")
for tag, n in census.most_common(40):
    print(f"  {tag}: {n}")

print("\nLevel-1 children:")
for c in root:
    print(" ", short(c.tag), c.attrib.get("Name", ""), len(list(c)))
    for cc in list(c)[:6]:
        print("    ", short(cc.tag), cc.attrib.get("Name", ""), len(list(cc)))

# sample: the SynapseTypes subtree and any element named like a synapse instance
for el in root.iter():
    if short(el.tag) == "SynapseTypes":
        print("\nSynapseTypes children:")
        for st in el:
            print("  ", short(st.tag), st.attrib)
            for p in st:
                print("      ", short(p.tag), p.attrib, (p.text or "").strip()[:40])
        break
