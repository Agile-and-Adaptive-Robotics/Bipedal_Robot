"""Second probe: AnimatLab stores nodes as <Node Type="...">. Find the
neuron/synapse type strings + name/id/link structure."""
import re

xml = open(r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
           r'\Walker_2_Layer_CPG_BilateralRG'
           r'\Walker_2_Layer_CPG_BilateralRG.aproj',
           encoding='utf-8', errors='ignore').read()
out = []
types = re.findall(r'<Node Type="([^"]+)"', xml)
from collections import Counter
c = Counter(types)
for t, n in c.most_common(30):
    out.append(f"{n:5d}  {t}")
# sample a neuron-ish node block
m = re.search(r'<Node Type="[^"]*Neuron[^"]*"[^>]*>.{0,700}', xml, re.S)
if m:
    out.append("\nSAMPLE Neuron node:\n" + m.group(0)[:700])
m2 = re.search(r'<Node Type="[^"]*Synapse[^"]*"[^>]*>.{0,700}', xml, re.S)
if m2:
    out.append("\nSAMPLE Synapse node:\n" + m2.group(0)[:700])
open("_w2l_probe2.txt", "w", encoding="utf-8").write("\n".join(out))
print("wrote _w2l_probe2.txt")
