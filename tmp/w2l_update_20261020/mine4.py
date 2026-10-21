"""Find where synapse instances live + dump one."""
import xml.etree.ElementTree as ET

APROJ = r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj"
tree = ET.parse(APROJ)
root = tree.getroot()
parent = {c: p for p in root.iter() for c in p}

for el in root.iter():
    if el.find("SynapticTypeID") is not None:
        chain = []
        cur = el
        while cur in parent and len(chain) < 8:
            chain.append(cur.tag)
            cur = parent[cur]
        print("synapse element tag:", el.tag)
        print("chain:", " <- ".join(chain))
        print(ET.tostring(el, encoding="unicode")[:2500])
        break
