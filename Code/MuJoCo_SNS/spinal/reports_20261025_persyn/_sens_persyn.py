"""Gain-sensitivity probe + remove-and-retest for the persyn w10 winner
(ruling: anti-lottery validation of the per-conductance winner).

1. Sensitivity: every knob individually +/-10% (in-process 20 s evals);
   report the objective delta per knob -> which conductances matter.
2. Remove-and-retest: the 8 LEAST influential knobs zeroed/anchored one at
   a time (retest) — a knob whose removal does not degrade the walk is a
   lottery-weight candidate (reported, never silently dropped).

Usage: python _sens_persyn.py"""
from __future__ import annotations

import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
WEAN = os.path.normpath(os.path.join(
    HERE, "..", "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
import _w2l_wean_eval as EV  # noqa: E402
import numpy as np           # noqa: E402

d = json.load(open(os.path.join(HERE, "round_w10.json"), encoding="utf-8"))
P, WSUP = d["params"], d["wsup"]
EV.T.KV["net"] = "bilateralrg2"


def run(p, tag):
    kn = dict(EV.BASE, dscale=1.0, te=3.0, tf=4.0, wsup=WSUP, sagittal=50.0,
              dur_walk=20.0, **p)
    return EV.evaluate(kn, tag=tag,
                       csv_path=os.path.join(HERE, "sens_trials.csv"))


base = run(P, "base")
print(f"base: obj={base['obj']:.1f} survive={base['survive']} "
      f"cyc {base['cycL']}/{base['cycR']}", flush=True)
deltas = {}
for k in sorted(P):
    for sgn, lab in ((1.10, "+10%"), (0.90, "-10%")):
        p2 = dict(P, **{k: P[k] * sgn})
        m = run(p2, f"sens_{k}_{lab}")
        deltas[(k, lab)] = float(m["obj"]) - float(base["obj"])
    print(f"  {k:16s} d(+10%)={deltas[(k,'+10%')]:+7.1f} "
          f"d(-10%)={deltas[(k,'-10%')]:+7.1f}", flush=True)

mag = {k: abs(deltas[(k, "+10%")]) + abs(deltas[(k, "-10%")])
       for k in P}
ranked = sorted(mag, key=mag.get)
print("\nleast influential (remove-and-retest candidates):", ranked[:8])
anchor = json.load(open(os.path.join(HERE, "round_w70.json"),
                        encoding="utf-8"))  # documented anchor values
retest = {}
for k in ranked[:8]:
    p2 = dict(P, **{k: anchor["params"][k]})
    m = run(p2, f"rm_{k}")
    retest[k] = float(m["obj"]) - float(base["obj"])
    print(f"  REMOVE {k:16s} -> obj delta {retest[k]:+7.1f} "
          f"survive={m['survive']}", flush=True)
json.dump(dict(base_obj=float(base["obj"]),
               deltas={f"{k}|{l}": v for (k, l), v in deltas.items()},
               magnitude=mag, retest=retest),
          open(os.path.join(HERE, "sensitivity_w10.json"), "w"), indent=1)
print("wrote sensitivity_w10.json")
