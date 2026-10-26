"""High-budget per-synapse search (2026-10-26, Ben's "high gear" directive).

Uses the parallel_eval pool (10 workers, wraps the proven run_phase path)
for 7.5x throughput. All 61 conductance keys are genuinely live (the
2026-10-26 forwarding fix is active). Budget: 500 trials at stage 1,
then the weight-support wean ladder at 100 trials per rung (properly
funded this time; the 140-trial redo was under-budget and its winner
failed fresh verification).

Study: curr_w2l_persyn_v2 + wean rungs in optuna_w2l_persyn_v2.db.
"""
from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "reports_20261016", "w2l_wean"))
sys.path.insert(0, HERE)

import _w2l_wean_eval as EV       # noqa: E402
import test_w2l_ground as T       # noqa: E402
import optuna                     # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)
EV.T.KV["net"] = "bilateralrg2"

# Forwarding fix: all searched keys KV-resident
from _w2l_persyn_campaign import SPACE, ANCHORS, SEEDS  # noqa: E402
for _k in SPACE:
    T.KV.setdefault(_k, float(ANCHORS[_k]))

DB = "sqlite:///" + os.path.join(HERE, "optuna_w2l_persyn_v2.db").replace("\\", "/")
PY = sys.executable
LADDER = [0.55, 0.40, 0.25, 0.10]
TE_PIN, TF_PIN = 3.0, 4.0
SAG = 50.0


def knobs_of(wsup, p, dur):
    return dict(EV.BASE, dscale=1.0, te=TE_PIN, tf=TF_PIN,
                wsup=wsup, sagittal=SAG,
                **{k: float(p[k]) for k in SPACE}, dur_walk=dur)


def fresh_verify(wsup, params):
    import subprocess
    kn = knobs_of(wsup, params, 30.0)
    tmp = os.path.join(HERE, "_finalist_v2.json")
    json.dump(kn, open(tmp, "w"), indent=1)
    r = subprocess.run(
        [PY, os.path.join(HERE, "_verify_finalist.py"), tmp],
        capture_output=True, text=True, timeout=1800)
    for line in r.stdout.splitlines():
        if line.startswith("RESULT "):
            return json.loads(line[7:])
    raise RuntimeError(f"fresh verify failed:\n{r.stdout[-500:]}")


def run_round(wsup, seed_params, ntrials):
    name = ("curr_w2l_persyn_v2" if abs(wsup - 0.70) < 1e-9 else
            f"curr_w2l_persyn_v2_w{int(round(wsup*100)):02d}")
    study = optuna.create_study(
        study_name=name, storage=DB, direction="maximize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=20261026))
    have = len(study.trials)
    if have == 0:
        for sp in ([seed_params] if seed_params else []) + SEEDS:
            study.enqueue_trial(sp)
    need = max(ntrials - have, 0)

    def obj(trial):
        p = {k: trial.suggest_float(k, lo, hi)
             for k, (lo, hi) in SPACE.items()}
        m = EV.evaluate(knobs_of(wsup, p, 20.0),
                        tag=f"{name}#{trial.number}",
                        csv_path=os.path.join(HERE, "v2_trials.csv"))
        for k, v in m.items():
            if not k.startswith("_") and v is not None:
                try:
                    trial.set_user_attr(k, v)
                except (TypeError, ValueError):
                    pass
        return float(m["obj"])

    if need > 0:
        study.optimize(obj, n_trials=need)
    best = study.best_trial
    out = dict(wsup=wsup, study=name, n_trials=len(study.trials),
               params={k: round(float(best.params[k]), 5) for k in SPACE},
               obj20=float(best.value))
    m30 = fresh_verify(wsup, best.params)
    out["m30_fresh"] = {k: m30.get(k) for k in EV.METRIC_KEYS}
    out["survive30"] = bool(m30.get("survive"))
    fn = f"round_v2_w{int(round(wsup*100)):02d}.json"
    json.dump(out, open(os.path.join(HERE, fn), "w"), encoding="utf-8"),
    print(f"== ROUND wsup={wsup:.2f} obj20={out['obj20']:.1f} "
          f"survive={out['survive30']} -> {fn}", flush=True)
    return out


def main():
    n_stage1 = 500
    for a in sys.argv[1:]:
        if a.startswith("--n="):
            n_stage1 = int(a[4:])

    t0 = time.time()
    r70 = run_round(0.70, None, n_stage1)
    print(f"stage 1 done in {time.time()-t0:.0f} s", flush=True)

    if not r70["survive30"]:
        print("stage-1 winner failed fresh verify; trying top-5 ...")
        # get top 5 trials and fresh-verify each until one survives
        import optuna
        study = optuna.load_study(
            study_name="curr_w2l_persyn_v2", storage=DB)
        trials = sorted(study.trials, key=lambda t: t.value, reverse=True)
        for t in trials[:5]:
            m = fresh_verify(0.70, t.params)
            if m.get("survive"):
                r70 = dict(wsup=0.70, study="curr_w2l_persyn_v2",
                          params={k: round(float(v), 5)
                                  for k, v in t.params.items()},
                          obj20=t.value,
                          m30_fresh={k: m.get(k) for k in EV.METRIC_KEYS},
                          survive30=True)
                json.dump(r70, open(os.path.join(
                    HERE, "round_v2_w70.json"), "w"), indent=1)
                print(f"  top-#{t.number} SURVIVES fresh verify", flush=True)
                break
        else:
            print("no top-5 trial survives; honest stop", flush=True)
            return

    prev = r70["params"]
    for ws in LADDER:
        r = run_round(ws, prev, 100)
        if not r["survive30"]:
            print(f"WEAN DEATH at wsup={ws:.2f}", flush=True)
            return
        prev = r["params"]
    print("WEAN COMPLETE through wsup=0.10", flush=True)


if __name__ == "__main__":
    main()
