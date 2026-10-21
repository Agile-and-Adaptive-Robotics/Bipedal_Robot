"""Find Name elements' parents + full NonSpiking node scan for name-ish fields."""
import xml.etree.ElementTree as ET
from collections import Counter

APROJ = r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj"
tree = ET.parse(APROJ)
root = tree.getroot()
parent = {c: p for p in root.iter() for c in p}

pc = Counter()
samples = {}
for el in root.iter("Name"):
    p = parent[el]
    gp = parent.get(p)
    key = f"{p.tag}/{gp.tag if gp is not None else '-'}"
    pc[key] += 1
    samples.setdefault(key, el.text)
for k, v in pc.most_common():
    print(f"{v:5d}  {k}   e.g. {samples[k]!r}")

# Also: does a NonSpiking Node contain a Name anywhere in its subtree?
for node in root.iter("Node"):
    cls = (node.findtext("ClassName") or "")
    if cls.endswith("Neurons.NonSpiking"):
        names = [el.text for el in node.iter("Name")]
        texts = [el.text for el in node.iter("Text")]
        print("\nNonSpiking subtree Name elements:", names[:5], " Text elements:", texts[:5])
        # what does InLinks hold?
        inl = node.find("InLinks")
        if inl is not None:
            print("InLinks children:", [c.tag for c in inl][:10])
        break

# NodeToPhysicalAdapter dump (motor path) + where muscle names are
for node in root.iter("Node"):
    cls = (node.findtext("ClassName") or "")
    if "Muscle" in cls:
        s = ET.tostring(node, encoding="unicode")
        import re
        m = re.findall(r"<Name>([^<]*)</Name>", s)
        print("\nMuscle node Names:", m[:3])
        break
