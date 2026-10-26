"""STRICT basin gate for the condsolve finalist (Batch BW Tier 2, 2026-10-24).
Contact2 campaign's _basin_contact2.py pattern: +/-1% uniform perturbation
of EVERY searched knob (conductance scales included), 12 FRESH-process
trials, survival must persist in >= 11/12 (the acceptance bar; the best
contact2 result was 9/12).

Usage: python _basin_condsolve.py [winner.json]
"""
from __future__ import annotations

import io
import json
import os
import subprocess
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
import numpy as np  # noqa: E402

src = sys.argv[1] if len(sys.argv) > 1 else "condsolve_best.json"
d = json.load(open(os.path.join(HERE, src), encoding="utf-8"))
P = d["params"]

rng = np.random.default_rng(20261024)
survived = 0
rows = []
for k in range(12):
    kn = {"dscale": 1.0, "dur_walk": 30.0,
          "legcoll": 0.0, "sink": 0.06, "rig": 0.25, "walker_yield": 1.5,
          "te": 3.0, "tf": 4.0}
    for key, v in P.items():
        f = 1.0 + rng.uniform(-0.01, 0.01)
        kn[key] = float(v) * float(f)
    tmp = os.path.join(HERE, "_basin_knobs.json")
    json.dump(kn, open(tmp, "w"))
    r = subprocess.run([PY, os.path.join(HERE, "_verify_finalist.py"), tmp],
                       capture_output=True, text=True, timeout=1800)
    res = None
    for line in r.stdout.splitlines():
        if line.startswith("RESULT "):
            res = json.loads(line[7:])
    ok = bool(res and res.get("survive"))
    survived += ok
    rows.append((k, ok, None if res is None else (res["cycL"], res["cycR"],
                                                  res["phase_med"])))
    print(f"  trial {k:2d}: survive={ok} "
          f"{rows[-1][2] if rows[-1][2] else ''}", flush=True)

verdict = "PASS" if survived >= 11 else "FAIL"
print(f"\nBASIN {verdict}: {survived}/12 survived (bar 11/12; contact2 best was 9/12)")
out = dict(verdict=verdict, survived=survived, bar="11/12",
           rows=[(k, ok, list(v) if v else None) for k, ok, v in rows])
json.dump(out, open(os.path.join(HERE, "basin_condsolve.json"), "w"), indent=1)
