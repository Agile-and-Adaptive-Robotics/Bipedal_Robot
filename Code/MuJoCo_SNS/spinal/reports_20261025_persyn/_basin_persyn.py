"""Strict basin gate for the persyn w10 winner (ruling: validation of the
per-conductance winner, never a reason to shrink the space).
+/-1% uniform perturbation of ALL 61 knobs, 12 FRESH-process trials at
wsup=0.10, bar 11/12. Usage: python _basin_persyn.py"""
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

d = json.load(open(os.path.join(HERE, "round_w25.json"), encoding="utf-8"))
P = d["params"]
WSUP = d["wsup"]

rng = np.random.default_rng(20261025)
survived, rows = 0, []
for k in range(12):
    kn = {"dscale": 1.0, "dur_walk": 30.0,
          "legcoll": 0.0, "sink": 0.06, "rig": 0.25, "walker_yield": 1.5,
          "te": 3.0, "tf": 4.0, "wsup": WSUP, "sagittal": 50.0}
    for key, v in P.items():
        kn[key] = float(v) * float(1.0 + rng.uniform(-0.01, 0.01))
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
print(f"\nBASIN {verdict}: {survived}/12 (bar 11/12)")
json.dump(dict(verdict=verdict, survived=survived,
               rows=[(k, ok, list(v) if v else None) for k, ok, v in rows]),
          open(os.path.join(HERE, "basin_w25.json"), "w"), indent=1)
