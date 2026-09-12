"""Bayesian optimization (Optuna TPE) of the ground-walking spinal network.

Motivation: SNS-Toolbox ships no optimizer; the sample-efficient standard
for expensive gait simulations is Bayesian optimization (Calandra 2014
biped BO; Antonova 2016; Ryu & Geyer 2021 CPG optimality). TPE handles our
~60 s/eval budget; sqlite storage makes runs resumable.

Each trial: mutate params.py module dicts in place (build_network binds the
SAME dict objects at import, so mutations propagate), run runner.main in
--eval mode (12 s ground walk, afferents ON, semi-supported rig), and score:

  reward  forward COM progress under the tether (dx)
          alternating bursts (burst_r)
          swing knee flexion depth |knee_min|
          hip range (hip_amp)
  penalty falling (knee... COM height < 0.70 soft, < 0.62 hard)
          trunk lean (tilt_max)
          network instability (NaN)

Usage: python optuna_walk.py [n_trials]      (resumable; sqlite keeps prior
trials - delete optuna_walk.db to restart fresh)
"""
from __future__ import annotations

import json
import sys

import optuna

import params
import runner as R

DB = "sqlite:///optuna_walk.db"
STUDY = "ground_walk_v2"


def set_params(p: dict):
    """Write one trial's parameters into the params module dicts.

    build_network imported G/TAU/W_* by reference (dict objects), so
    in-place mutation is visible to every later build()."""
    params.TAU["rg_adapt"] = p["rg_adapt"]
    params.G["descend_to_rg_e"] = p["desc_e"]
    params.G["rg_to_pf"] = p["rg_to_pf"]
    params.W_PF_MN["E2"]["ankle_pf"] = p["e2_pf"]
    params.W_PF_MN["F1"]["ankle_df"] = p["f1_df"]
    params.W_PF_MN["F1"]["knee_flex"] = p["f1_kf"]
    params.W_POSTURE["knee_ext"] = p["post_kneext"]
    params.W_POSTURE["hip_ext"] = p["post_hipext"]
    params.BAL["kx"] = p["kx"]


def objective(trial: optuna.Trial) -> float:
    p = dict(
        drive=trial.suggest_float("drive", 1.5, 3.5),
        rg_adapt=trial.suggest_float("rg_adapt", 0.8, 2.2),
        desc_e=trial.suggest_float("desc_e", 0.8, 1.6),
        rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        e2_pf=trial.suggest_float("e2_pf", 0.1, 0.4),
        f1_df=trial.suggest_float("f1_df", 0.4, 0.8),
        f1_kf=trial.suggest_float("f1_kf", 1.4, 2.2),
        post_kneext=trial.suggest_float("post_kneext", 0.1, 0.5),
        post_hipext=trial.suggest_float("post_hipext", 0.05, 0.5),
        kx=trial.suggest_float("kx", 100.0, 300.0),
    )
    set_params(p)
    m = R.main(["--eval", "--drive", f"{p['drive']:.4f}"])

    if m["nan"]:
        # unstable: partial credit for how long it survived
        return -50.0 + m["t_end"]

    # v2 objective: REQUIRES real swing flexion (capped reward saturates at
    # 35 deg, so a stiff-legged shuffle cannot win) and hip range
    knee_rew = 0.25 * min(35.0, abs(m["knee_min"]))
    score = (
        8.0 * m["dx"]                     # forward progress under the tether
        + 0.25 * m["burst_r"]             # alternating rhythm
        + knee_rew
        + 0.10 * min(50.0, m["hip_amp"])  # hip range (sagittal)
        - 0.4 * max(0.0, m["tilt_max"] - 20.0)   # trunk lean penalty
        - 6.0 * max(0.0, 0.70 - m["kz"])  # soft height penalty
    )
    if m["kz"] < 0.62:                    # hard fall
        score -= 10.0
    return score


def main(argv):
    n_trials = int(argv[0]) if argv and argv[0].isdigit() else 40
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        direction="maximize", storage=DB, study_name=STUDY,
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=5, n_startup_trials=8))
    # seed v2 with the v1 winner where the keys overlap
    try:
        with open("best_walk_params.json", encoding="utf-8") as f:
            v1 = json.load(f)["params"]
        seed = {k: v for k, v in v1.items() if k != "post_hipext"}
        seed["post_hipext"] = 0.22
        study.enqueue_trial(seed)
        print("enqueued v1-best as trial seed")
    except FileNotFoundError:
        pass
    study.optimize(objective, n_trials=n_trials, gc_after_trial=True)

    best = study.best_trial
    print(f"\n== best score {best.value:.3f} (trial {best.number})")
    for k, v in best.params.items():
        print(f"  {k:12s} {v:.4f}")
    with open("best_walk_params.json", "w", encoding="utf-8") as f:
        json.dump({"score": best.value, "params": best.params,
                   "trial": best.number, "study": STUDY}, f, indent=2)
    print("saved best_walk_params.json")


if __name__ == "__main__":
    main(sys.argv[1:])
