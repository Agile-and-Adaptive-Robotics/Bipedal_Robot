# Per-instance vs type-level conductance in the Deng 2019 paper project.
import json
from pathlib import Path

d = json.load(open(
    r"D:\Github\Bipedal_Robot\Code\Matlab\SNS_Simscape\dev\deng2019_mine_20261025.json",
    encoding="utf-8"))
dev = [e for e in d["edges"]
       if e["instance_amp"] is not None and e["type_gmax"] is not None
       and abs(e["instance_amp"] - e["type_gmax"]) / max(e["type_gmax"], 1e-12) > 0.01]
print("edges whose per-instance SynapticConductance deviates >1% from its "
      "type gmax:", len(dev), "of", len(d["edges"]))
for e in dev[:12]:
    print("  %-22s -> %-22s type=%-22s type_g=%-11s inst=%s"
          % (e["src"], e["dst"], e["type"], e["type_gmax"], e["instance_amp"]))
# the PF->MN edges with names, for the twin comparison
print("--- PF->MN excitation edges (paper values) ---")
for e in d["edges"]:
    if e["type"] and "MN Excite" in e["type"]:
        print("  %-28s %-28s g=%-11s eref=%s inst=%s"
              % (e["src"], e["dst"], e["type_gmax"], e["eref"], e["instance_amp"]))
