"""Diagnostic: build the spiking mirror in a FRESH process at the same
TUNED config and print the true counts + edge-class census, to determine
whether the 8092-synapse harvest inherited state mutated by the variant
curriculum loaders inside make_editor_templates.main()."""
import io, os, sys
from collections import Counter

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
os.chdir(SPINAL)
sys.path.insert(0, SPINAL)
assert "AARL_NET" not in os.environ

import params as P
import build_network_spiking as BSS
import make_editor_templates as MT

keep = dict(P.G)
P.G.update(MT.SPIKE_TUNED)
P.G["joint_pf"] = 0.0
acts = []
import muscle_map as MM
for b in MM._GROUPS_BY_NAME:
    acts.append(b + "_r")
    acts.append(b + "_l")
net = BSS.build(acts, interleg=True)
P.G.clear(); P.G.update(keep)

nobj = net.net
names = [q["name"] for q in nobj.populations]
print("FRESH-PROCESS build: %d neurons / %d inputs / %d synapses"
      % (nobj.get_num_neurons(), nobj.get_num_inputs_actual(),
         nobj.get_num_connections()))
print("RD taps:", sum(1 for n in names if n.startswith("RD_")))

census = Counter()
for c in nobj.connections:
    s, d = names[c["source"]], names[c["destination"]]
    ks, kd = MT._spk_kind(s)[0], MT._spk_kind(d)[0]
    sign = "exc" if c["params"].get("reversal_potential", 0) > -1e-6 else "inh"
    census[(ks, kd, sign)] += 1
for k in sorted(census):
    print("  %-28s %5d" % ("%s->%s %s" % k, census[k]))
print("total:", sum(census.values()))
