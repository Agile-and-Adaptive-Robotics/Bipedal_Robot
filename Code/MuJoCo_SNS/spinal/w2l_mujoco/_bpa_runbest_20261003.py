"""Run one BPA protocol config and save the npz (figure input).

Usage: python _bpa_runbest_20261003.py <ns|spk> <json_params> <npz_out>
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SPINAL = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, SPINAL)

import _curriculum_bpa as NS  # noqa: E402  (also chdirs to HERE)

which, pj, npz_out = sys.argv[1], sys.argv[2], sys.argv[3]
p = dict(NS.SEED)
if os.path.exists(pj):
    p.update(json.loads(open(pj, encoding="utf-8").read())["params"])
if which == "ns":
    import build_w2l_split_net as B
    B.NAP["tau_max_h"] = p["tau"]
    net = B.build(comm=p["comm"],
                  template_path=os.path.join(SPINAL,
                                             "connectome_templates.json"))
else:
    import build_w2l_split_net_spiking as BS
    net = BS.build(comm=p.get("comm", 1.0), tau=p.get("tau", 0.25),
                   thr_inc=p.get("thr_inc", 4.0),
                   template_path=os.path.join(SPINAL,
                                              "connectome_templates.json"))
met = NS.run_bpa(net, p, save_npz=npz_out, verbose=True)
print("score", NS.score_bpa(met))
