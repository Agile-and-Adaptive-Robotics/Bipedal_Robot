"""Audit probe 3: syn6 motif census (correct names) + antagonist
projections; s3k default-build motif census for comparison."""
import io
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, str(SPINAL))
os.chdir(SPINAL)


def dump(n):
    pop_names = [p["name"] for p in n.populations]
    edges = defaultdict(list)
    for c in n.connections:
        key = (pop_names[c["source"]], pop_names[c["destination"]])
        syn = c["params"]
        edges[key].append((float(syn.get("max_conductance", np.nan)),
                           "exc" if float(syn.get(
                               "reversal_potential", 0)) > -1e-6 else "inh"))
    return edges


import muscle_map as MM  # noqa: E402

acts = []
for base in MM._GROUPS_BY_NAME:
    acts.append(base + "_r")
    acts.append(base + "_l")

print("== syn6 motif census ==")
os.environ["AARL_NET"] = "syn6"
import importlib  # noqa: E402
import build_network as bn  # noqa: E402
net = bn.build(acts)
edges = dump(net.net)
missing_ia = [a for a in acts
              if (f"Ia_{a}", f"MN_{a}") not in edges]
missing_ibin = [a for a in acts
                if (f"IBIN_{a}", f"MN_{a}") not in edges]
missing_iix = [a for a in acts
               if (f"IIX_{a}", f"MN_{a}") not in edges]
print("muscles missing Ia->MN:", missing_ia)
print("muscles missing IBIN->MN:", missing_ibin)
print("muscles missing IIX->MN:", missing_iix)
for a in ("iliacus_r", "psoas_r"):
    print(a, "Ia->MN:", edges.get((f"Ia_{a}", f"MN_{a}")))
# antagonist projection sample: knee ext vs flex on r
print("IaIN_vas_lat_r -> MN_semiten_r:",
      edges.get(("IaIN_vas_lat_r", "MN_semiten_r")))
print("IIIN_vas_lat_r -> MN_bifemsh_r:",
      edges.get(("IIIN_vas_lat_r", "MN_bifemsh_r")))
print("IaIN_vas_lat_r <-> IaIN_vas_med_r:",
      edges.get(("IaIN_vas_lat_r", "IaIN_vas_med_r")),
      edges.get(("IaIN_vas_med_r", "IaIN_vas_lat_r")))
print("IBIN_vas_lat_r <-> IBIN_vas_med_r:",
      edges.get(("IBIN_vas_lat_r", "IBIN_vas_med_r")))
n_ia_ant = sum(1 for (s, d) in edges
               if s.startswith("IaIN_") and d.startswith("MN_"))
n_iiin_ant = sum(1 for (s, d) in edges
                 if s.startswith("IIIN_") and d.startswith("MN_"))
n_iain_mut = sum(1 for (s, d) in edges
                 if s.startswith("IaIN_") and d.startswith("IaIN_")
                 and s != d)
n_ibin_mut = sum(1 for (s, d) in edges
                 if s.startswith("IBIN_") and d.startswith("IBIN_")
                 and s != d)
print(f"IaIN->ant MN edges: {n_ia_ant}, IIIN->ant MN: {n_iiin_ant}, "
      f"IaIN<->IaIN: {n_iain_mut}, IBIN<->IBIN: {n_ibin_mut}")
del os.environ["AARL_NET"]
importlib.reload(bn)
