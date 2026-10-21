"""Fresh-process winner verification for the w2l_contact2 campaign.

Usage: python _verify_finalist.py <knobs.json>
Prints 'RESULT {json}' with the _w2l_wean_eval METRIC_KEYS."""
from __future__ import annotations

import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
WEAN = os.path.normpath(os.path.join(
    HERE, "..", "..", "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import _w2l_wean_eval as EV  # noqa: E402

kn = json.load(open(sys.argv[1], encoding="utf-8"))
kn.pop("net", None)
EV.T.KV["net"] = "bilateralrg2"
m = EV.evaluate(kn, tag="fresh_verify",
                csv_path=os.path.join(HERE, "verify_trials.csv"))
print("RESULT " + json.dumps({k: EV._j(m[k]) for k in EV.METRIC_KEYS}))
