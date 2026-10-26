"""Di Russo stage-2 WITH replay bootstrap + proper budget (2026-10-26).

Changes from the first stage-2 (dirusso_stage2.py):
1. REPLAY BOOTSTRAP: first 2 s impose reference kinematics (her protocol)
2. 500 TRIALS (was 120, under-funded)
3. Uses the parallel_eval pool where possible (the Di Russo model is
   larger than the W2L so each trial is slower, but parallelism still helps)

Objective: survival gate + speed toward 1.2 m/s + cycle count + antiphase
(her stage-2 dynamic stability goal, same as before but now starting
from a walking-initialized state).
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import controller as C          # noqa: E402
import dirusso_run as R         # noqa: E402

DB = "sqlite:///" + os.path.join(HERE, "optuna_dirusso_s2r.db").replace("\\", "/")
STUDY = "dirusso_s2r_replay"
DUR = 15.0   # 2 s replay + 13 s free
Z_FALL, TILT_FALL = 0.75, 45.0
SPEED_TARGET = 1.2

SPACE = {}
for m in C.MUSCLES:
    SPACE[f"ff_{m}"] = (0.3, 3.0)
    for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
              "ib_rev", "rc_mn", "cut_exc"):
        SPACE[f"{g}_{m}"] = (0.0, 2.0)
SPACE["omega"] = (3.5, 7.0)     # widened around the ref cycle 5.11 rad/s
SPACE["gamma"] = (0.0, 2.0)
SPACE["wsup"] = (0.3, 1.0)      # harness support as a search dim


def seed_params():
    p = {"omega": 2 * np.pi / 1.23, "gamma": 1.0, "wsup": 0.7}
    for m in C.MUSCLES:
        p[f"ff_{m}"] = 1.0
        for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
                  "ib_rev", "rc_mn", "cut_exc"):
            p[f"{g}_{m}"] = 0.1
    return p


def evaluate(p, tag=""):
    W0 = np.load(os.path.join(HERE, "dirusso_W_stage1.npy"))
    W = W0.copy()
    gains = {}
    for m in C.MUSCLES:
        i = C.MIDX[m]
        W[i, :] *= p[f"ff_{m}"]           # L
        W[9 + i, :] *= p[f"ff_{m}"]       # R
        for g in ("ia_mono", "ia_recip", "ii_exc", "ii_inh", "ib_auto",
                  "ib_rev", "rc_mn", "cut_exc"):
            if g not in gains:
                gains[g] = np.zeros(C.NM)
            gains[g][i] = p[f"{g}_{m}"]
    log, strikes = R.run(W, dur=DUR, wsup=p["wsup"],
                         omega=p["omega"], gamma=p["gamma"],
                         gains_in=gains, replay_s=2.0,
                         out=f"_s2r_{tag}.npz")
    z = log["z"]; t = log["t"]; x = log["x"]
    free = t >= 3.0   # metrics after replay + settle
    zf, tf, xf = z[free], t[free], x[free]
    if len(tf) < 2:
        return dict(obj=-100.0, survive=False, cycles=0, speed=0, zmin=0)
    slope = np.polyfit(tf, xf, 1)[0]
    cycL = log["cycleL"][-1]; cycR = log["cycleR"][-1]
    tilt = np.degrees(log["tilt"])
    finite = np.isfinite(z).all()
    fell = (not finite) or z.min() < Z_FALL or np.abs(tilt).max() > TILT_FALL
    ph = None
    if cycL >= 3 and cycR >= 3:
        rT = np.array(strikes["R"]); lT = np.array(strikes["L"])
        phs = []
        for i in range(len(rT) - 1):
            ins = lT[(lT >= rT[i]) & (lT < rT[i + 1])]
            if len(ins) == 1:
                phs.append((ins[0] - rT[i]) / (rT[i + 1] - rT[i]))
        if phs:
            ph = float(np.median(phs))
    survive = (finite and not fell and min(cycL, cycR) >= 3
               and ph is not None and 0.35 <= ph <= 0.65)
    obj = -100.0 if (not finite or fell) else \
        5 * min(min(cycL, cycR), 10) \
        + 60 * (1 - min(abs(slope - SPEED_TARGET) / SPEED_TARGET, 1.0)) \
        + (20 * (1 - min(abs(ph - 0.5) / 0.5, 1.0)) if ph is not None else 0) \
        - 0.5 * max(np.abs(tilt).max() - 16.0, 0.0)
    return dict(obj=obj, survive=survive, cycles=min(cycL, cycR),
                speed=float(slope), zmin=float(z.min()) if finite else None,
                tiltmax=float(np.abs(tilt).max()) if finite else None,
                phase=ph)


def main():
    n = 500
    for a in sys.argv[1:]:
        if a.startswith("--n="):
            n = int(a[4:])
    study = optuna.create_study(study_name=STUDY, storage=DB,
                                direction="maximize", load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(seed=20261026))
    if not study.trials:
        study.enqueue_trial(seed_params())

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi)
             for k, (lo, hi) in SPACE.items()}
        m = evaluate(p, tag=str(trial.number))
        for k, v in m.items():
            if v is not None:
                trial.set_user_attr(k, v)
        print(f"[{STUDY}#{trial.number}] obj={m['obj']:.1f} "
              f"survive={m['survive']} cyc={m['cycles']} "
              f"speed={m['speed']:+.3f} z={m['zmin']}", flush=True)
        return float(m["obj"])

    t0 = time.time()
    study.optimize(obj, n_trials=n)
    b = study.best_trial
    json.dump(dict(params=b.params, value=b.value, user_attrs=b.user_attrs),
              open(os.path.join(HERE, "dirusso_s2r_best.json"), "w"), indent=1)
    print(f"done in {time.time()-t0:.0f} s | best {b.value:.1f}")


if __name__ == "__main__":
    main()
