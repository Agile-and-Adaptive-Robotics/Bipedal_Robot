r"""Di Russo stage-2 optimization (her eq. 12 -> dynamic stability; our
Optuna stands in for her CMA-ES, documented). D2, 2026-10-25 night.

Per the 2026-10-25 per-synapse ruling: EVERY searched weight is individual
(per muscle), seeded from stage-1 / the paper's pre-optimization bounding.
L/R shared by symmetry (Ben's reduction).

Search (83 dims): per-muscle feedforward scale (9, seed 1.0); per-muscle
individual reflex weights ia_mono, ia_recip, ii_exc, ii_inh, ib_auto,
ib_rev, rc_mn, cut_exc (9x8=72, seed 0.1 = the pre-opt bound); omega
(seed 5.71), gamma (seed 1.0).

Objective (her stage-2, adapted and labeled): survival gate (finite, no
fall = z>0.75 m and |tilt|<45 deg, >=6 heel-strike cycles/leg in 20 s,
inter-leg phase 0.35-0.65) then maximize speed toward 1.2 m/s + cycle
count + antiphase quality - tilt penalty. Harness wsup=0.7 fixed this
stage (her sim is unassisted; wean comes after stability exists).

Usage: python dirusso_stage2.py [--n=120]
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import controller as C          # noqa: E402
import dirusso_run as R         # noqa: E402

DB = "sqlite:///" + os.path.join(HERE, "optuna_dirusso_s2.db").replace("\\", "/")
STUDY = "dirusso_s2_stability"
DUR = 20.0
Z_FALL, TILT_FALL = 0.75, 45.0
SPEED_TARGET = 1.2

SPACE = {}
for m in C.MUSCLES:
    SPACE[f"ff_{m}"] = (0.3, 3.0)
    for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
              "ib_rev", "rc_mn", "cut_exc"):
        SPACE[f"{g}_{m}"] = (0.0, 2.0)
SPACE["omega"] = (4.0, 8.0)
SPACE["gamma"] = (0.0, 2.0)


def seed_params():
    p = {"omega": 2 * np.pi / 1.1, "gamma": 1.0}
    for m in C.MUSCLES:
        p[f"ff_{m}"] = 1.0
        for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
                  "ib_rev", "rc_mn", "cut_exc"):
            p[f"{g}_{m}"] = 0.1
    return p


def evaluate(p, tag):
    W0 = np.load(os.path.join(HERE, "dirusso_W_stage1.npy"))
    W = W0.copy()
    gains = {}
    for m in C.MUSCLES:
        i = C.MIDX[m]
        W[0:9][i] *= p[f"ff_{m}"]        # L rows
        W[9:18][i] *= p[f"ff_{m}"]       # R rows (symmetry)
        for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
                  "ib_rev", "rc_mn", "cut_exc"):
            if g not in gains:
                gains[g] = np.zeros(C.NM)
            gains[g][i] = p[f"{g}_{m}"]
    log, strikes = R.run(W, dur=DUR, wsup=0.7, omega=p["omega"],
                         gamma=p["gamma"], gains_in=gains,
                         out=f"_s2_{tag}.npz")
    z = log["z"]; t = log["t"]; x = log["x"]
    seg = slice(len(z) // 3, None)
    slope = np.polyfit(t[seg], x[seg], 1)[0]
    cycL = log["cycleL"][-1]; cycR = log["cycleR"][-1]
    tilt = np.degrees(log["tilt"])
    finite = np.isfinite(z).all()
    fell = (not finite) or z.min() < Z_FALL or np.abs(tilt).max() > TILT_FALL
    ph = None
    # inter-leg phase: median L-strike phase within consecutive R strikes
    if cycL >= 6 and cycR >= 6:
        rT = np.array(strikes["R"]); lT = np.array(strikes["L"])
        phs = []
        for i in range(len(rT) - 1):
            inside = lT[(lT >= rT[i]) & (lT < rT[i + 1])]
            if len(inside) == 1:
                phs.append((inside[0] - rT[i]) / (rT[i + 1] - rT[i]))
        if phs:
            ph = float(np.median(phs))
    survive = (finite and not fell and min(cycL, cycR) >= 6
               and ph is not None and 0.35 <= ph <= 0.65)
    obj = -100.0 if (not finite or fell) else \
        5 * min(min(cycL, cycR), 12) \
        + 60 * (1 - min(abs(slope - SPEED_TARGET) / SPEED_TARGET, 1.0)) \
        + (20 * (1 - min(abs(ph - 0.5) / 0.5, 1.0)) if ph is not None else 0) \
        - 0.5 * max(np.abs(tilt).max() - 16.0, 0.0)
    return dict(obj=obj, survive=survive, cycL=int(cycL), cycR=int(cycR),
                phase=ph, speed=float(slope), zmin=float(z.min()) if finite else None,
                tiltmax=float(np.abs(tilt).max()) if finite else None)


def main():
    n = 120
    for a in sys.argv[1:]:
        if a.startswith("--n="):
            n = int(a[4:])
    study = optuna.create_study(study_name=STUDY, storage=DB,
                                direction="maximize", load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(seed=20261025))
    if not study.trials:
        study.enqueue_trial(seed_params())

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi) for k, (lo, hi) in SPACE.items()}
        m = evaluate(p, tag=str(trial.number))
        for k, v in m.items():
            if v is not None:
                trial.set_user_attr(k, v)
        print(f"[{STUDY}#{trial.number}] survive={m['survive']} obj={m['obj']:.1f} "
              f"cyc {m['cycL']}/{m['cycR']} phase {m['phase']} "
              f"speed {m['speed']:+.3f} zmin {m['zmin']}", flush=True)
        return float(m["obj"])

    study.optimize(obj, n_trials=n)
    b = study.best_trial
    json.dump(dict(params=b.params, value=b.value,
                   user_attrs=b.user_attrs),
              open(os.path.join(HERE, "dirusso_s2_best.json"), "w"), indent=1)
    print("best obj", b.value, "-> dirusso_s2_best.json")


if __name__ == "__main__":
    main()
