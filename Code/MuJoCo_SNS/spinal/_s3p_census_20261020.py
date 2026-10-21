"""Defaults regression census after the 2026-10-20 s3p edits (coact knob,
mechanism offsets, platform hook, ctrl log). GATE: (410, 376, 1186) at
defaults; coact>0 must ADD synapses (antagonist posture edges)."""
import io, sys
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import muscle_map as MM
acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
import build_network as BN
import params as P

net = BN.build(acts)
n = net.net
c0 = (n.get_num_neurons(), n.get_num_inputs_actual(),
      n.get_num_connections())
print("defaults census:", c0)
print("GATE:", "PASS" if c0 == (410, 376, 1186) else "FAIL")
P.G["coact"] = 0.8
net2 = BN.build(acts)
n2 = net2.net
c1 = (n2.get_num_neurons(), n2.get_num_inputs_actual(),
      n2.get_num_connections())
print("coact=0.8 census:", c1, "delta synapses:", c1[2] - c0[2])
P.G["coact"] = 0.0
net3 = BN.build(acts)
n3 = net3.net
c2 = (n3.get_num_neurons(), n3.get_num_inputs_actual(),
      n3.get_num_connections())
print("coact back to 0 census:", c2,
      "RESTORE:", "PASS" if c2 == (410, 376, 1186) else "FAIL")
