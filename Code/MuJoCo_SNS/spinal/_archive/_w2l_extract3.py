"""Third probe: list the distinct XML tags in the aproj and find where
neurons/synapses actually live."""
import re
from collections import Counter

xml = open(r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
           r'\Walker_2_Layer_CPG_BilateralRG'
           r'\Walker_2_Layer_CPG_BilateralRG.aproj',
           encoding='utf-8', errors='ignore').read()
tags = Counter(t for t in re.findall(r'<([A-Za-z_][A-Za-z0-9_]*)', xml))
lines = [f"{n:6d}  {t}" for t, n in tags.most_common(50)]
lines.append("")
# case-insensitive hunt for synapse-ish content
for kw in ("ynapse", "Neuron", "NeuralModule", "Link", "Connex"):
    idxs = [m.start() for m in re.finditer(kw, xml)][:3]
    lines.append(f"--- '{kw}': {len(re.findall(kw, xml))} hits")
    for ix in idxs:
        lines.append("   ..." + xml[max(0, ix-80):ix+160]
                     .replace("\n", " ")[:220])
open("_w2l_probe3.txt", "w", encoding="utf-8").write("\n".join(lines))
print("wrote _w2l_probe3.txt")
