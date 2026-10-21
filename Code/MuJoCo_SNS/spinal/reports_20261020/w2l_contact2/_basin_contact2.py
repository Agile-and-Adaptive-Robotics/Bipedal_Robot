"""Basin gate for the w2l_contact2 finalist (the wean campaign's
_w2l_wean_basin.py pattern, FRESH-PROCESS trials): +/-1% uniform
perturbation of every SPACE knob, 12 trials, survival must persist in
>= 11/12 (the wean record used 12/12; 11/12 = one-blip tolerance stated
up front).

Usage: python _basin_contact2.py [round_best_rXXX.json]
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

src = sys.argv[1] if len(sys.argv) > 1 else "round_best_r030.json"
d = json.load(open(os.path.join(HERE, src), encoding="utf-8"))
P = d["params"]
DS = d["dscale"]

rng = np.random.default_rng(20261020)
survived = 0
rows = []
for k in range(12):
    kn = {"dscale": DS, "dur_walk": 30.0,
          "legcoll": 0.0, "sink": 0.06, "rig": 0.25, "walker_yield": 1.5}
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
print(f"BASIN {verdict} {survived}/12 (finalist {src}, dscale {DS})")
json.dump({"source": src, "survived": survived, "trials": rows},
          open(os.path.join(HERE, "basin_contact2.json"), "w"), indent=1)
sys.exit(0 if survived >= 11 else 1)
