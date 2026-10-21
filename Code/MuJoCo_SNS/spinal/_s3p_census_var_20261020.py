"""Variant census regression after the coact patch (defaults must hold:
w2lvar 888/382/7866, syn6 794/382/3033; coact>0 must add edges)."""
import io, sys, os
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
import muscle_map as MM
acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")
import params as P

os.environ["AARL_NET"] = "w2lvar"
import build_network_w2lvar as BW
n = BW.build(acts).net
c0 = (n.get_num_neurons(), n.get_num_inputs_actual(), n.get_num_connections())
print("w2lvar defaults:", c0, "GATE:", "PASS" if c0 == (888, 382, 7866) else "FAIL")
P.G["coact"] = 0.8
n = BW.build(acts).net
c1 = (n.get_num_neurons(), n.get_num_inputs_actual(), n.get_num_connections())
print("w2lvar coact=0.8:", c1, "delta:", c1[2] - c0[2])
P.G["coact"] = 0.0

os.environ["AARL_NET"] = "syn6"
P.G["syn6"] = 1.0; P.G["syn6_brainstem"] = 0.0
import build_network_syn6 as BS
n = BS.build(acts).net
d0 = (n.get_num_neurons(), n.get_num_inputs_actual(), n.get_num_connections())
print("syn6 defaults:", d0, "GATE:", "PASS" if d0 == (794, 382, 3033) else "FAIL")
P.G["coact"] = 0.8
n = BS.build(acts).net
d1 = (n.get_num_neurons(), n.get_num_inputs_actual(), n.get_num_connections())
print("syn6 coact=0.8:", d1, "delta:", d1[2] - d0[2])
P.G["coact"] = 0.0
