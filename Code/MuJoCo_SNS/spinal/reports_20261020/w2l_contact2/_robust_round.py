"""Robust-objective round for w2l_contact2 (basin-gate-driven): TPE over
the same SPACE, objective = MEAN obj of the config + two +/-1%
perturbations (worst-case mean proxy; the basin gate requires >=11/12
survivals, the marginal-basin fix). All runs in THIS process (the study
objective); the finalist then gets the fresh-process basin gate.

Usage: python _robust_round.py [n_trials=36]
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WEAN = os.path.normpath(os.path.join(HERE, "..", "..",
                                     "reports_20261016", "w2l_wean"))
for p in (WEAN, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import numpy as np                      # noqa: E402
import optuna                           # noqa: E402
import _w2l_wean_eval as EV             # noqa: E402
from _w2l_contact2_campaign import (DB, SPACE, SEED, knobs_of)  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)
EV.T.KV["net"] = "bilateralrg2"

RNG = np.random.default_rng(SEED)
STUDY = "curr_w2l_contact2_robust"


def robust_obj(p):
    """Mean obj over base + 2 perturbed copies (3 x 20 s runs)."""
    vals = []
    for pert in (None, -0.01, 0.01):
        pp = dict(p)
        if pert is not None:
            lo, hi = (pert, pert * 0.6) if pert < 0 else (pert * 0.6, pert)
            for k in pp:
                pp[k] = float(pp[k]) * (1.0 + RNG.uniform(lo, hi))
        m = EV.evaluate(knobs_of(1.0, pp, 20.0), tag="rob")
        vals.append(float(m["obj"]))
    return float(np.mean(vals))


def main() -> int:
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 36
    study = optuna.create_study(study_name=STUDY, storage=DB,
                                direction="maximize", load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(seed=SEED))
    have = len(study.trials)
    if have == 0:
        for f in ("round_best_r030.json", "round_best_bendefaults.json",
                  "round_best_r100.json"):
            fp = os.path.join(HERE, f)
            if os.path.exists(fp):
                study.enqueue_trial(json.load(open(fp))["params"])

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi)
             for k, (lo, hi) in SPACE.items()}
        v = robust_obj(p)
        trial.set_user_attr("obj_mean", v)
        return v

    if n - have > 0:
        study.optimize(obj, n_trials=n - have)
    best = study.best_trial
    out = dict(study=STUDY, n_trials=len(study.trials),
               params={k: round(float(v), 5) for k, v in best.params.items()},
               obj_mean=float(best.value))
    json.dump(out, open(os.path.join(HERE, "robust_best.json"), "w"),
              indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
