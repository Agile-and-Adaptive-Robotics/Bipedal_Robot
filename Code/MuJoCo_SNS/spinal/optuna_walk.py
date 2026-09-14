"""Bayesian optimization (Optuna TPE) of the ground-walking spinal network.

v4 (2026-09-12): KINEMATICS-MATCH objective (Ben: "keep fine-tuning until
the kinematics are similar to OpenSim"). The score is kine_ref.compare's
kine_score - cycle-normalized hip/knee/ankle shape RMSE + peak-knee +
range + stance-duty errors against subject01_walk1_ik.mot phased by the
measured GRF - with hard gates for falling/NaN/tilt (kinematics of a
fallen model are meaningless). v3's stability-shaped objective kept
winning with stiff shuffles; v4 makes the OpenSim match the target.

New knobs vs v3:
    desc_f     DRIVE -> RG-F conductance (swing-side drive; duty lever)
    e2_adapt   PF_SHAPE["E2"] adaptation multiplier (push-off window
               length -> stance duty)
Baseline = fitted_walk_params.json (fit_pf.py back-solve refit) exactly
like runner --fitted; weight knobs are multipliers on it, pf_gain scales
the whole table (back-solved weights ~10x smaller than hand-tuned).

Usage: python optuna_walk.py [n_trials]      (resumable; sqlite keeps
prior trials - delete optuna_walk.db to restart fresh)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import optuna

import params
import runner as R

HERE = Path(__file__).parent
DB = "sqlite:///optuna_walk.db"
STUDY = "ground_walk_v4b_kine"

FIT = None
BASE = dict(
    e2_pf=params.W_PF_MN["E2"]["ankle_pf"],
    f1_df=params.W_PF_MN["F1"]["ankle_df"],
    f1_kf=params.W_PF_MN["F1"]["knee_flex"],
    post_kneext=params.W_POSTURE["knee_ext"],
    post_hipext=params.W_POSTURE["hip_ext"],
    desc_f=params.G["descend_to_rg_f"],
    e2_adapt=params.PF_SHAPE["E2"][1],
)


def load_fitted_baseline():
    """Apply fitted W_PF_MN/W_POSTURE tables + adopt their knob values."""
    global FIT, BASE
    fit = json.loads((HERE / "fitted_walk_params.json").read_text("utf-8"))
    FIT = dict(W_PF_MN={p: dict(t) for p, t in fit["W_PF_MN"].items()},
               W_POSTURE=dict(fit["W_POSTURE"]))
    for ph, tbl in fit["W_PF_MN"].items():
        for g, w in tbl.items():
            params.W_PF_MN[ph][g] = float(w)
    for g, w in fit["W_POSTURE"].items():
        params.W_POSTURE[g] = float(w)
    BASE["e2_pf"] = params.W_PF_MN["E2"]["ankle_pf"]
    BASE["f1_df"] = params.W_PF_MN["F1"]["ankle_df"]
    BASE["f1_kf"] = params.W_PF_MN["F1"]["knee_flex"]
    BASE["post_kneext"] = params.W_POSTURE["knee_ext"]
    BASE["post_hipext"] = params.W_POSTURE["hip_ext"]
    print("fitted baseline loaded: " +
          " ".join(f"{k}={v:.3f}" for k, v in BASE.items()))


def set_params(p: dict):
    """Write one trial's parameters into the params module dicts.

    build_network binds the SAME dict objects at import, so in-place
    mutation propagates to every later build(). Weight knobs are
    multipliers on the fitted baseline; pf_gain scales the whole table;
    desc_f / e2_adapt are ABSOLUTE conductance / multiplier values."""
    gain = p["pf_gain"]
    for ph, tbl in FIT["W_PF_MN"].items():
        for g, w in tbl.items():
            params.W_PF_MN[ph][g] = w * gain
    for g, w in FIT["W_POSTURE"].items():
        params.W_POSTURE[g] = w * gain
    params.TAU["rg_adapt"] = p["rg_adapt"]
    params.G["descend_to_rg_e"] = p["desc_e"]
    params.G["descend_to_rg_f"] = p["desc_f"]
    params.G["rg_to_pf"] = p["rg_to_pf"]
    params.PF_SHAPE["E2"] = (params.PF_SHAPE["E2"][0], p["e2_adapt"])
    params.W_PF_MN["E2"]["ankle_pf"] = BASE["e2_pf"] * gain * p["e2_pf"]
    params.W_PF_MN["F1"]["ankle_df"] = BASE["f1_df"] * gain * p["f1_df"]
    params.W_PF_MN["F1"]["knee_flex"] = BASE["f1_kf"] * gain * p["f1_kf"]
    params.W_POSTURE["knee_ext"] = BASE["post_kneext"] * gain * p["post_kneext"]
    params.W_POSTURE["hip_ext"] = BASE["post_hipext"] * gain * p["post_hipext"]
    params.BAL["kx"] = p["kx"]


def objective(trial: optuna.Trial) -> float:
    p = dict(
        drive=trial.suggest_float("drive", 1.2, 3.2),
        rg_adapt=trial.suggest_float("rg_adapt", 0.8, 2.4),
        desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
        desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
        rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        pf_gain=trial.suggest_float("pf_gain", 0.4, 4.0, log=True),
        e2_pf=trial.suggest_float("e2_pf", 0.4, 2.0),
        f1_df=trial.suggest_float("f1_df", 0.4, 2.0),
        f1_kf=trial.suggest_float("f1_kf", 0.4, 2.2),
        e2_adapt=trial.suggest_float("e2_adapt", 0.8, 2.6),
        post_kneext=trial.suggest_float("post_kneext", 0.4, 1.8),
        post_hipext=trial.suggest_float("post_hipext", 0.4, 1.8),
        kx=trial.suggest_float("kx", 100.0, 300.0),
    )
    set_params(p)
    # NOTE: no --fitted here - load_fitted_baseline() already applied the
    # tables; a reload inside main() would wipe the trial's mutations
    m = R.main(["--eval", "--drive", f"{p['drive']:.4f}"])

    if m["nan"]:
        return -80.0 + m["t_end"]
    # v4b landscape: a frozen/no-rhythm run scores like the zero-motion
    # model it is (~-65, the kine penalty of not moving - the first v4
    # attempt used a -25 sentinel and TPE collapsed onto that plateau,
    # since every genuine walking attempt scored WORSE than -25)
    if m.get("kine") is None:
        return -65.0
    score = float(m["kine_score"])
    if m["kz"] < 0.62:                     # hard fall: kinematics invalid
        score -= 10.0
    if m["tilt_max"] > 40.0:               # limbo: shape no longer comparable
        score -= 5.0
    return score


def main(argv):
    n_trials = int(argv[0]) if argv and argv[0].isdigit() else 60
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    load_fitted_baseline()
    study = optuna.create_study(
        direction="maximize", storage=DB, study_name=STUDY,
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=11, n_startup_trials=10))
    # seed: the v3 winner expressed in v4 coordinates (its weight
    # multipliers = 1 at its own pf_gain; desc_f / e2_adapt untouched)
    seed = dict(drive=2.2, rg_adapt=1.9, desc_e=1.7, desc_f=1.4,
                rg_to_pf=2.4, pf_gain=1.0, e2_pf=1.0, f1_df=1.0, f1_kf=1.0,
                e2_adapt=1.0, post_kneext=1.0, post_hipext=1.0, kx=150.0)
    try:
        prev = json.loads((HERE / "best_walk_params.json").read_text("utf-8"))
        if prev.get("study") == "ground_walk_v3":
            mul = prev.get("multipliers", {})
            seed.update(drive=prev["params"]["drive"],
                        rg_adapt=prev["params"]["rg_adapt"],
                        desc_e=prev["params"]["desc_e"],
                        rg_to_pf=prev["params"]["rg_to_pf"],
                        kx=prev["params"]["kx"],
                        pf_gain=prev.get("pf_gain", 1.0))
            print("seeded v3 winner into v4")
    except (FileNotFoundError, KeyError):
        pass
    study.enqueue_trial(seed)
    print("enqueued seed trial")
    study.optimize(objective, n_trials=n_trials, gc_after_trial=True)

    best = study.best_trial
    print(f"\n== best kine_score {best.value:.3f} (trial {best.number})")
    for k, v in best.params.items():
        print(f"  {k:12s} {v:.4f}")
    set_params(best.params)
    g = best.params["pf_gain"]
    eff = dict(drive=best.params["drive"],
               rg_adapt=best.params["rg_adapt"],
               desc_e=best.params["desc_e"],
               desc_f=best.params["desc_f"],
               rg_to_pf=best.params["rg_to_pf"],
               e2_adapt=best.params["e2_adapt"],
               e2_pf=BASE["e2_pf"] * g * best.params["e2_pf"],
               f1_df=BASE["f1_df"] * g * best.params["f1_df"],
               f1_kf=BASE["f1_kf"] * g * best.params["f1_kf"],
               post_kneext=BASE["post_kneext"] * g * best.params["post_kneext"],
               post_hipext=BASE["post_hipext"] * g * best.params["post_hipext"],
               kx=best.params["kx"])
    with open(HERE / "best_walk_params.json", "w", encoding="utf-8") as f:
        json.dump({"score": best.value, "params": eff, "pf_gain": g,
                   "multipliers": best.params, "baseline": dict(BASE),
                   "trial": best.number, "study": STUDY}, f, indent=2)
    print("saved best_walk_params.json (runner --fitted --best reproduces)")


if __name__ == "__main__":
    main(sys.argv[1:])
