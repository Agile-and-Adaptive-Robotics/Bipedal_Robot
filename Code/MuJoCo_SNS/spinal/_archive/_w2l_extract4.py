"""Fourth probe: sample full <Node> and <Link> blocks to design the
aproj -> block-editor JSON extractor."""
import re

xml = open(r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
           r'\Walker_2_Layer_CPG_BilateralRG'
           r'\Walker_2_Layer_CPG_BilateralRG.aproj',
           encoding='utf-8', errors='ignore').read()
out = []
nodes = re.findall(r'<Node [^>]*>', xml)
out.append(f"{len(nodes)} <Node> tags; first 8 verbatim:")
for n in nodes[:8]:
    out.append("  " + n[:200])
links = re.findall(r'<Link [^>]*>', xml)
out.append(f"{len(links)} <Link> tags; first 5 verbatim:")
for l in links[:5]:
    out.append("  " + l[:220])
# ClassNames inventory
cn = re.findall(r'<Node [^>]*ClassName="([^"]+)"', xml)
from collections import Counter
out.append("ClassNames: " + repr(Counter(cn).most_common(20)))
# where are synapse types?
sn = re.findall(r'<Node [^>]*ClassName="([^"]*[Ss]ynapse[^"]*)"', xml)
out.append(f"synapse-class nodes: {Counter(sn)}")
open("_w2l_probe4.txt", "w", encoding="utf-8").write("\n".join(out))
print("wrote _w2l_probe4.txt")
