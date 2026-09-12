"""Bayesian optimization (Optuna TPE) of the ground-walking spinal network.

v3 (2026-09-11 evening): the balanced-pattern rerun. Ben's staged plan after
the v2 lesson ("stability penalties dominate motion rewards; the missing
piece is balanced activations, not more parameter search"):
  - baseline = fitted_walk_params.json (fit_pf.py: W_PF_MN / W_POSTURE
    per-group NNLS refit from the IK/NNLS back-solved human pattern)
  - the five weight knobs search RELATIVE multipliers [0.5, 1.8] around
    that baseline instead of absolute ranges
  - objective adds an E-duty term targeting the human 0.60 stance duty
    (runner --eval now reports `duty`)

Each trial: mutate params.py module dicts in place (build_network binds the
SAME dict objects at import, so mutations propagate), run runner.main in
--eval mode (12 s ground walk, afferents ON, semi-supported rig), score:

  reward  forward COM progress under the tether (dx)
          alternating bursts (burst_r)
          swing knee flexion depth |knee_min| (saturates at 35 deg)
          hip range (hip_amp), E-duty near 0.60
  penalty falling (COM height < 0.70 soft, < 0.62 hard)
          trunk lean (tilt_max), network instability (NaN)

Usage: python optuna_walk.py [n_trials]      (resumable; sqlite keeps prior
trials - delete optuna_walk.db to restart fresh)
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
STUDY = "ground_walk_v3"

# fitted baseline tables (loaded in main), plus knob values
FIT = None
BASE = dict(
    e2_pf=params.W_PF_MN["E2"]["ankle_pf"],
    f1_df=params.W_PF_MN["F1"]["ankle_df"],
    f1_kf=params.W_PF_MN["F1"]["knee_flex"],
    post_kneext=params.W_POSTURE["knee_ext"],
    post_hipext=params.W_POSTURE["hip_ext"],
)


def load_fitted_baseline():
    """Apply fitted W_PF_MN/W_POSTURE tables + adopt their 5 knob values."""
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

    build_network imported G/TAU/W_* by reference (dict objects), so
    in-place mutation is visible to every later build(). pf_gain scales
    the WHOLE fitted table (the back-solved weights carry honest human
    amplitudes ~0.1-0.3 while the network's gain structure was hand-tuned
    against a much larger table - the global gain finds the operating
    point); the 5 knob multipliers then refine individual entries."""
    gain = p["pf_gain"]
    for ph, tbl in FIT["W_PF_MN"].items():
        for g, w in tbl.items():
            params.W_PF_MN[ph][g] = w * gain
    for g, w in FIT["W_POSTURE"].items():
        params.W_POSTURE[g] = w * gain
    params.TAU["rg_adapt"] = p["rg_adapt"]
    params.G["descend_to_rg_e"] = p["desc_e"]
    params.G["rg_to_pf"] = p["rg_to_pf"]
    params.W_PF_MN["E2"]["ankle_pf"] = BASE["e2_pf"] * gain * p["e2_pf"]
    params.W_PF_MN["F1"]["ankle_df"] = BASE["f1_df"] * gain * p["f1_df"]
    params.W_PF_MN["F1"]["knee_flex"] = BASE["f1_kf"] * gain * p["f1_kf"]
    params.W_POSTURE["knee_ext"] = BASE["post_kneext"] * gain * p["post_kneext"]
    params.W_POSTURE["hip_ext"] = BASE["post_hipext"] * gain * p["post_hipext"]
    params.BAL["kx"] = p["kx"]


def objective(trial: optuna.Trial) -> float:
    p = dict(
        drive=trial.suggest_float("drive", 1.2, 3.2),
        rg_adapt=trial.suggest_float("rg_adapt", 0.8, 2.2),
        desc_e=trial.suggest_float("desc_e", 0.8, 1.6),
        rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        pf_gain=trial.suggest_float("pf_gain", 0.5, 8.0, log=True),
        e2_pf=trial.suggest_float("e2_pf", 0.5, 1.8),
        f1_df=trial.suggest_float("f1_df", 0.5, 1.8),
        f1_kf=trial.suggest_float("f1_kf", 0.5, 1.8),
        post_kneext=trial.suggest_float("post_kneext", 0.5, 1.8),
        post_hipext=trial.suggest_float("post_hipext", 0.5, 1.8),
        kx=trial.suggest_float("kx", 100.0, 300.0),
    )
    set_params(p)
    # NOTE: no --fitted here - load_fitted_baseline() at startup already
    # applied the tables, and a --fitted reload inside main() would wipe
    # this trial's mutations (verified: identical metrics with/without)
    m = R.main(["--eval", "--drive", f"{p['drive']:.4f}"])

    if m["nan"]:
        # unstable: partial credit for how long it survived
        return -50.0 + m["t_end"]

    # v2 shaping + v3 duty term (human stance duty ~0.60)
    knee_rew = 0.25 * min(35.0, abs(m["knee_min"]))
    duty_rew = 1.0 * max(0.0, 0.30 - abs(m.get("duty", 0.3) - 0.60))
    score = (
        8.0 * m["dx"]                     # forward progress under the tether
        + 0.25 * m["burst_r"]             # alternating rhythm
        + knee_rew
        + duty_rew
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
    load_fitted_baseline()
    study = optuna.create_study(
        direction="maximize", storage=DB, study_name=STUDY,
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=7, n_startup_trials=8))
    # seed with the fitted baseline itself (multipliers = 1); pf_gain ~3
    # bridges the hand-tuned-vs-human amplitude gap as a starting guess
    seed = dict(drive=2.2, rg_adapt=1.9, desc_e=1.7, rg_to_pf=2.4,
                pf_gain=3.0, e2_pf=1.0, f1_df=1.0, f1_kf=1.0,
                post_kneext=1.0, post_hipext=1.0, kx=150.0)
    try:
        best = json.loads((HERE / "best_walk_params.json").read_text("utf-8"))
        if best.get("study") == "ground_walk_v2":
            seed.update(drive=best["params"]["drive"],
                        rg_adapt=best["params"]["rg_adapt"],
                        desc_e=best["params"]["desc_e"],
                        rg_to_pf=best["params"]["rg_to_pf"],
                        kx=best["params"]["kx"])
            print("seeded v2-best drive/rg/kx on top of the fitted baseline")
    except (FileNotFoundError, KeyError):
        pass
    study.enqueue_trial(seed)
    print("enqueued fitted-baseline seed trial")
    study.optimize(objective, n_trials=n_trials, gc_after_trial=True)

    best = study.best_trial
    print(f"\n== best score {best.value:.3f} (trial {best.number})")
    for k, v in best.params.items():
        print(f"  {k:12s} {v:.4f}")
    set_params(best.params)
    g = best.params["pf_gain"]
    eff = dict(drive=best.params["drive"],
               rg_adapt=best.params["rg_adapt"],
               desc_e=best.params["desc_e"],
               rg_to_pf=best.params["rg_to_pf"],
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
    print(f"saved best_walk_params.json (effective values + pf_gain={g:.3f};"
          " run `runner --fitted --best`: --fitted loads the table, --best"
          " applies the gain + knobs + drive)")


if __name__ == "__main__":
    main(sys.argv[1:])
