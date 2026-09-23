"""Fifth probe: ClassName VALUES (child elements) + hunt for SNS
neuron-class element names."""
import re
from collections import Counter

xml = open(r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
           r'\Walker_2_Layer_CPG_BilateralRG'
           r'\Walker_2_Layer_CPG_BilateralRG.aproj',
           encoding='utf-8', errors='ignore').read()
out = []
cn = re.findall(r'<ClassName>([^<]+)</ClassName>', xml)
out.append("ClassName values: " + repr(Counter(cn).most_common(25)))
# element names that look neuron-ish
tags = Counter(t for t in re.findall(r'<([A-Za-z_][A-Za-z0-9_.]*)', xml)
               if any(k in t.lower() for k in
                      ("neuron", "spik", "synaps", "node", "ion", "cell")))
out.append("neuron-ish tags: " + repr(tags.most_common(30)))
# sample
for pat in (r'<[A-Za-z_.]*Neuron[A-Za-z_.]*[^>]*>.{0,500}',
            r'<Spiking[^>]*>.{0,400}'):
    m = re.search(pat, xml)
    if m:
        out.append("SAMPLE:\n" + m.group(0)[:500])
        break
open("_w2l_probe5.txt", "w", encoding="utf-8").write("\n".join(out))
print("wrote _w2l_probe5.txt")
