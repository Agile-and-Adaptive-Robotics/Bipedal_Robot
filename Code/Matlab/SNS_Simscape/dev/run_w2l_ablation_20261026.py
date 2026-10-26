# Order A (Ben 2026-10-26): feedback intact vs ablated, ONE frozen tuning,
# walking B2L model, desktop harness. Predictions P4-P6 pre-registered in
# SNS_Simscape\reports_20261024\PREREG_REPRO_ABLATION_20261026.md.
# Ablated arm zeroes every harness-exposed sensory channel: contact
# cascades (ch_hc, ch_ank, ch_contra), the contact encoder (ctn_amp=0),
# and the Ib force-sensor current (ibnA=0). The Li-style Ia/II chains are
# NOT knobbed in this harness and stay at template defaults - labeled.
from __future__ import annotations
import json
import os
import sys

HERE = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "..", "Code", "MuJoCo_SNS", "spinal",
    "reports_20261020", "w2l_contact2"))
WEAN = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "..", "..", "..", "Code", "MuJoCo_SNS", "spinal",
    "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

OUT = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "results"))
os.makedirs(OUT, exist_ok=True)

arm = sys.argv[1]          # intact | ablated
seed = dict(te=3.0, tf=4.0, ch_hc=0.08, ch_ank=1.0, ch_contra=0.08,
            ctn_amp=4.0, cref=50.0, ibnA=1.0, comm=1.0)   # Ben's defaults
if arm == "ablated":
    seed.update(ch_hc=0.0, ch_ank=0.0, ch_contra=0.0, ctn_amp=0.0, ibnA=0.0)

os.environ["AARL_NPZ"] = os.path.join(
    OUT, f"ablation_{arm}_20261026.npz").replace("\\", "/")

import _w2l_wean_eval as EV  # noqa: E402
EV.T.KV["net"] = "bilateralrg2"

kn = dict(EV.BASE, dscale=1.0,
          te=seed["te"], tf=seed["tf"],
          ch_hc=seed["ch_hc"], ch_ank=seed["ch_ank"],
          ch_contra=seed["ch_contra"], ctn_amp=seed["ctn_amp"],
          cref=seed["cref"], ibnA=seed["ibnA"], comm=seed["comm"],
          dur_walk=30.0)
import io  # noqa: E402
_real = sys.stdout
sys.stdout = io.StringIO()          # shield: the runner stack closes the
m = EV.evaluate(kn, tag=f"ablation_{arm}_20261026")  # original stdout
sys.stdout = _real
res = {k: (v if isinstance(v, (int, float, bool, str)) else str(v))
       for k, v in (m.items() if hasattr(m, "items") else vars(m).items())}
res["arm"] = arm
res["knobs"] = seed
with open(os.path.join(OUT, f"ablation_{arm}_20261026.json"), "w") as f:
    json.dump(res, f, indent=1)
print(f"ARM {arm}: " + json.dumps(res, default=str)[:600])
