"""CLOSED-LOOP basin gate for the wean finalist (basin_gate.py pattern
adapted): the network-only constant-DRIVE gate is inapplicable because
this config family's rhythm is schedule/closed-loop sustained (network-
only baseline decays to zero — basin_round8_best3.json.json, and the
2026-09-20 s3b finding). Here: perturb every scalar param +/-1% (basin_
gate.perturb), run the FULL ground eval (the loop the winner actually
operates in), and require the gait to persist (bilateral >=3+3 cycles).
Usage: python basin_closed_20261016.py <winner.json> [--trials 8]
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import wean_tune_20261016 as WT
import basin_gate as BG

ap = argparse.ArgumentParser()
ap.add_argument("winner")
ap.add_argument("--trials", type=int, default=8)
ap.add_argument("--sigma", type=float, default=0.01)
ap.add_argument("--out", default=None)
args = ap.parse_args()

win = json.loads(open(args.winner, encoding="utf-8").read())
scale = float(win.get("scale", 1.0))
knobs = win["knobs"]
fb = {k: float(win["params"].get(k, 0.0)) for k in WT.FB_KEYS}

# snapshot the applied state (params tables) so each trial restores
# EXACTLY the winner before perturbing
WT.apply_config(knobs, scale, fb=fb)
P, R = WT.MODS["P"], WT.MODS["R"]
import params

snapshot = {}
for name in ("G", "TAU", "W_POSTURE", "BAL", "PF_SHAPE", "W_PF_MN",
             "AFF", "MOD"):
    snapshot[name] = json.dumps(getattr(params, name),
                                default=lambda o: list(o)
                                if isinstance(o, tuple) else o)


def restore():
    for name in snapshot:
        tbl = getattr(params, name)
        tbl.clear()
        d = json.loads(snapshot[name])
        for k, v in d.items():
            if isinstance(v, list):
                v = tuple(v)
            elif isinstance(v, dict):
                v = {g: (tuple(w) if isinstance(w, list) else w)
                     for g, w in v.items()}
            tbl[k] = v


def eval_once():
    m = R.main(["--eval", "--drive", repr(WT.BASE_DRIVE)])
    kk = m.get("kine") or {}
    return dict(kine_score=m.get("kine_score"),
                nan=bool(m.get("nan")),
                cyc_r=int(kk.get("n_cycles_r", 0) or 0),
                cyc_l=int(kk.get("n_cycles_l", 0) or 0),
                duty_r=kk.get("duty_r"), duty_l=kk.get("duty_l"),
                kz=m.get("kz"))


base = eval_once()
print(f"baseline (unperturbed): {json.dumps(base)}", flush=True)

rows = []
for t in range(args.trials):
    restore()
    BG.perturb(args.sigma, np.random.default_rng(3000 + t))
    res = eval_once()
    ok = (not res["nan"] and res["cyc_r"] >= 3 and res["cyc_l"] >= 3)
    res["pass"] = bool(ok)
    rows.append(res)
    print(f"  trial {t}: {json.dumps(res)}", flush=True)

report = dict(winner=os.path.abspath(args.winner), scale=scale,
              sigma=args.sigma, trials=args.trials, baseline=base,
              pass_fraction=sum(r["pass"] for r in rows) / len(rows),
              rows=rows)
out = args.out or os.path.join(
    HERE, "reports_20261016", "wean_human",
    f"basin_closed_{os.path.basename(args.winner)}.json")
json.dump(report, open(out, "w", encoding="utf-8"), indent=1)
print(f"--> closed-loop basin pass {sum(r['pass'] for r in rows)}/"
      f"{len(rows)}  wrote {out}", flush=True)
