"""Build fingerprint for the default-off regression gate.

Writes a stable digest of the compiled network: population count, input
count, connection count, and a SHA-256 over the sorted connection
(src_name, dst_name, sign, weight-hex) list plus population names+tau.
Default parameters only (G["joint_pf"] must be absent or 0 here).
"""
import hashlib
import io
import json
import sys

import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import build_network as bn
import params as P

P.G["joint_pf"] = 0.0
net = bn.build(["vas_lat_r", "semimem_r", "vas_lat_l", "semimem_l"],
               interleg=True)
names = [p["name"] for p in net.net.populations]
conns = []
for c in net.net.connections:
    conns.append((names[c["source"]], names[c["destination"]],
                  repr(float(c["params"].get("reversal_potential", 0.0))),
                  repr(float(c["params"].get("max_conductance", 0.0)))))
h = hashlib.sha256()
for row in sorted(conns):
    h.update(("|".join(row) + "\n").encode("utf-8"))
for nm in sorted(names):
    p = net.net.populations[names.index(nm)]
    tau = p.get("tau", None)
    h.update((nm + "|" + (repr(float(tau)) if tau is not None
                          else "?") + "\n").encode("utf-8"))
fp = {"n_pop": len(names), "n_conn": len(conns),
      "n_inputs": len(net.inputs), "digest": h.hexdigest()}
out = sys.argv[1] if len(sys.argv) > 1 else "_reg_fingerprint.json"
open(out, "w", encoding="utf-8").write(json.dumps(fp, indent=1))
print(json.dumps(fp, indent=1))
