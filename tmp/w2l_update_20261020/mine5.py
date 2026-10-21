"""Dump raw XML of one NonSpiking neuron + one PhysicalToNodeAdapter."""
import xml.etree.ElementTree as ET

APROJ = r"D:\GitHub\Bipedal_Robot\Neuromechanical_Models\Walker_2_Layer_CPG_BilateralRG\Walker_2_Layer_CPG_BilateralRG.aproj"
tree = ET.parse(APROJ)
root = tree.getroot()

done = set()
for node in root.iter("Node"):
    cls = (node.findtext("ClassName") or "")
    short = cls.split(".")[-1]
    if short in ("NonSpiking", "PhysicalToNodeAdapter") and short not in done:
        done.add(short)
        s = ET.tostring(node, encoding="unicode")
        # strip the drawing-format boilerplate children for readability
        print(f"===== {short} =====")
        print(s[:4000])
        print()
    if len(done) == 2:
        break
