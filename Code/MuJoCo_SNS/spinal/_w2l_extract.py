"""Probe the W2L BilateralRG .aproj: count neurons/synapses, list names.
Writes findings to _w2l_probe.txt (avoids inline-quoting issues)."""
import re

out = []
xml = open(r'D:\Github\Bipedal_Robot\Neuromechanical_Models'
           r'\Walker_2_Layer_CPG_BilateralRG'
           r'\Walker_2_Layer_CPG_BilateralRG.aproj',
           encoding='utf-8', errors='ignore').read()
out.append(f"size {len(xml)}")
out.append(f"Neuron tags: {len(re.findall(r'<Neuron ', xml))}")
out.append(f"Synapse tags: {len(re.findall(r'<Synapse ', xml))}")
neurons = re.findall(r'<Neuron [^>]*?Name="([^"]+)"', xml)
out.append(f"neuron names ({len(neurons)}): {neurons[:40]}")
syns = re.findall(r'<Synapse [^>]*?Name="([^"]+)"', xml)
out.append(f"synapse names ({len(syns)}): {syns[:25]}")
# AnimatLab links connect via IDs: findOrganism/adapter structure;
# sample one Synapse block
m = re.search(r'<Synapse .{0,900}', xml, re.S)
if m:
    out.append("sample synapse block:\n" + m.group(0)[:900])
open("_w2l_probe.txt", "w", encoding="utf-8").write("\n".join(out))
print("wrote _w2l_probe.txt")
