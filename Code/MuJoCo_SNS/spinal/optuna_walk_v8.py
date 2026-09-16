"""v6 study (ground_walk_v8_pose): v5 + swing-knee quad suppression.

Same kine objective and search space as optuna_walk_v5.py PLUS:
    f1_kneext_inh  F1 -> KINH inhibitory interneuron -> knee_ext MN pools
                   (phase-gated by F1 = swing only), range 0-2, default 0
                   = topology absent (v5-identical behavior).
Seeded with the v5 winner (best_walk_params_v7.json multipliers + the new
gain at 0). Writes best_walk_params_v8.json (runner --best6 loads it) and
never touches the v5/v4b jsons. Same CSV/full22-capture machinery.

Usage: python optuna_walk_v6.py [n_trials]   (default 100)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import optuna

import params
import runner as R

HERE = Path(__file__).parent
DB = "sqlite:///optuna_walk.db"
STUDY = "ground_walk_v8_pose"
CSV = HERE / "v8_results.csv"
BEST6 = HERE / "best_walk_params_v8.json"

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
FULL_SCHED = dict(params.SCHEDULE)
TRIAL_STATE: dict[int, dict] = {}
BEST_SO_FAR = {"value": float("-inf"), "trial": -1}
CSV_FIELDS = ["ts", "trial", "phase_reset_e", "phase_reset_f",
              "f1_kneext_inh", "run", "kine_score", "duty", "knee_min",
              "rmse_hip", "rmse_knee", "rmse_ankle", "cadence", "tilt_max",
              "stayed_up", "npz"]


def load_fitted_baseline():
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
          " ".join(f"{k}={v:.3f}" for k, v in BASE.items()), flush=True)


def set_params(p: dict):
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
    params.G["phase_reset_e"] = float(p.get("phase_reset_e", 0.0))
    params.G["phase_reset_f"] = float(p.get("phase_reset_f", 0.0))
    params.G["f1_kneext_inh"] = float(p.get("f1_kneext_inh", 0.0))
    params.G["f1_anklepf_inh"] = float(p.get("f1_anklepf_inh", 0.0))
    params.G["renshaw"] = 0.5


def csv_row(row: dict):
    import csv as _csv
    new = not CSV.exists()
    with open(CSV, "a", newline="", encoding="utf-8") as f:
        w = _csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)
    print("csv: " + " ".join(f"{k}={row.get(k)}" for k in
                             ("trial", "run", "kine_score", "duty",
                              "knee_min", "cadence", "stayed_up")),
          flush=True)


def metrics_row(trial_no: int, p: dict, m: dict, run: str,
                npz: str = "") -> dict:
    import numpy as np
    k = m.get("kine") or {}
    walk_dur = (params.SCHEDULE["walk"][1] - params.SCHEDULE["walk"][0])
    return dict(ts=time.strftime("%Y-%m-%d %H:%M:%S"), trial=trial_no,
                phase_reset_e=f"{p.get('phase_reset_e', 0.0):.4f}",
                phase_reset_f=f"{p.get('phase_reset_f', 0.0):.4f}",
                f1_kneext_inh=f"{p.get('f1_kneext_inh', 0.0):.4f}",
                run=run,
                kine_score=f"{m.get('kine_score', float('nan')):.3f}",
                duty=f"{k.get('duty', float('nan')):.3f}",
                knee_min=f"{k.get('knee_min', float('nan')):.2f}",
                rmse_hip=f"{k.get('rmse_hip', float('nan')):.2f}",
                rmse_knee=f"{k.get('rmse_knee', float('nan')):.2f}",
                rmse_ankle=f"{k.get('rmse_ankle', float('nan')):.2f}",
                cadence=f"{k.get('n_cycles', 0) / max(walk_dur, 1e-9):.3f}",
                tilt_max=f"{m.get('tilt_max', float('nan')):.1f}",
                stayed_up=bool(not m.get("nan") and m.get("kz", 0) > 0.62),
                npz=npz)


def full22_capture(study: optuna.Study, trial_no: int, p: dict):
    import numpy as np
    import kine_ref
    params.SCHEDULE.update(FULL_SCHED)
    set_params(p)
    R.main(["--drive", repr(p["drive"])])
    npz = HERE / f"v8_best_trial{trial_no}.npz"
    if (HERE / "spinal_run.npz").exists():
        npz.unlink(missing_ok=True)
        (HERE / "spinal_run.npz").rename(npz)
    z = np.load(npz, allow_pickle=True)
    t, q_deg, neuro, com = z["t"], z["q"], z["neuro"], z["com"]
    walk_start = FULL_SCHED["walk"][0]
    n_done = int(np.sum(t > 0)) + 1
    k = kine_ref.compare(t[:n_done], q_deg[:n_done], neuro[:n_done],
                         walk_start)
    m = dict(kine_score=(k["kine_score"] if k else float("nan")),
             tilt_max=float(np.nanmax(q_deg[:n_done, 0])),
             kz=float(np.nanmin(com[:n_done, 2])),
             nan=bool(not np.all(np.isfinite(q_deg[:n_done]))),
             kine=k)
    csv_row(metrics_row(trial_no, p, m, "full22", npz.name))
    print(f"full22 capture: trial {trial_no} kine_score "
          f"{m['kine_score']:.3f} -> {npz.name}", flush=True)


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
        phase_reset_e=trial.suggest_float("phase_reset_e", 0.0, 2.0),
        phase_reset_f=trial.suggest_float("phase_reset_f", 0.0, 2.0),
        f1_kneext_inh=trial.suggest_float("f1_kneext_inh", 0.0, 2.0),
        f1_anklepf_inh=trial.suggest_float("f1_anklepf_inh", 0.0, 2.0),
    )
    set_params(p)
    m = R.main(["--eval", "--drive", repr(p["drive"])])
    TRIAL_STATE[trial.number] = dict(params=p, metrics=m)
    if m["nan"]:
        return -80.0 + m["t_end"]
    if m.get("kine") is None:
        return -65.0
    score = float(m["kine_score"])
    if m["kz"] < 0.62:
        score -= 10.0
    if m["tilt_max"] > 40.0:
        score -= 5.0
    return score


def callback(study: optuna.Study, trial: optuna.Trial):
    st = TRIAL_STATE.get(trial.number)
    if st is not None:
        csv_row(metrics_row(trial.number, st["params"], st["metrics"],
                            "eval"))
    best = study.best_trial
    if best.number == trial.number and \
            (best.value or float("-inf")) > BEST_SO_FAR["value"] + 1e-9:
        BEST_SO_FAR["value"] = float(best.value)
        BEST_SO_FAR["trial"] = trial.number
        st = TRIAL_STATE.get(trial.number) or dict(params=trial.params)
        print(f"NEW GLOBAL BEST {best.value:.3f} (trial {trial.number}) "
              "- running full 22 s capture", flush=True)
        try:
            full22_capture(study, trial.number, st["params"])
        except Exception as e:
            print(f"full22 capture FAILED: {e}", flush=True)
        finally:
            params.SCHEDULE.update(FULL_SCHED)


def main(argv):
    n_trials = int(argv[0]) if argv and argv[0].isdigit() else 100
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    load_fitted_baseline()
    study = optuna.create_study(
        direction="maximize", storage=DB, study_name=STUDY,
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=13, n_startup_trials=10))
    if len(study.trials) == 0:
        seed = dict(drive=2.2, rg_adapt=1.9, desc_e=1.7, desc_f=1.4,
                    rg_to_pf=2.4, pf_gain=1.0, e2_pf=1.0, f1_df=1.0,
                    f1_kf=1.0, e2_adapt=1.0, post_kneext=1.0, post_hipext=1.0,
                    kx=150.0, phase_reset_e=0.0, phase_reset_f=0.0,
                    f1_kneext_inh=0.0)
        try:
            prev = json.loads(
                (HERE / "best_walk_params_v7.json").read_text("utf-8"))
            mul = prev.get("multipliers", {})
            seed.update({k: float(mul[k]) for k in seed if k in mul})
            seed["f1_kneext_inh"] = 0.0
            print("seeded v7 winner into v8 (new start pose)", flush=True)
        except (FileNotFoundError, KeyError):
            pass
        study.enqueue_trial(seed)
        print("enqueued seed trial", flush=True)
    done = [t for t in study.trials
            if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    BEST_SO_FAR["value"] = float(max(t.value for t in done)) if done else \
        float("-inf")
    BEST_SO_FAR["trial"] = -1
    study.optimize(objective, n_trials=n_trials, gc_after_trial=True,
                   callbacks=[callback])

    best = study.best_trial
    print(f"\n== best kine_score {best.value:.3f} (trial {best.number})",
          flush=True)
    for k, v in best.params.items():
        print(f"  {k:15s} {v:.4f}")
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
               kx=best.params["kx"],
               phase_reset_e=best.params["phase_reset_e"],
               phase_reset_f=best.params["phase_reset_f"],
               f1_kneext_inh=best.params["f1_kneext_inh"],
               f1_anklepf_inh=best.params["f1_anklepf_inh"],
               renshaw=0.5)
    with open(BEST6, "w", encoding="utf-8") as f:
        json.dump({"score": best.value, "params": eff,
                   "pf_gain": g, "multipliers": best.params,
                   "baseline": dict(BASE), "trial": best.number,
                   "study": STUDY}, f, indent=2)
    print("saved best_walk_params_v8.json (runner --best8 reproduces)",
          flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
