"""v5 phase-reset study (ground_walk_v5_phase) — overnight workhorse.

Same kine objective as v4b (optuna_walk.py; no-rhythm = -65, NaN -80 +
t_end, fall -10, tilt>40 -5) with TWO new dimensions:
    phase_reset_e  HIP_EXT_SIG gain (stance-prolonging hip-extensor
                   length signal -> RG-E exc / RG-F inh), range 0-2
    phase_reset_f  HIP_FLEX_SIG gain (swing-triggering hip-flexor
                   shortening-velocity signal -> RG-F exc / RG-E inh), 0-2
Both default 0 = exact v4b behavior (regression-gated: the study-path
reproduction of trial 49 is bit-exact, see DESIGN.md).

Search space = v4b params + the two gains; seeded with the v4b winner
(multipliers + gains 0). Resumable (sqlite; re-runs only enqueue the
seed when the study is empty). Writes best_walk_params_v5.json (runner
--best5 loads it) and NEVER touches best_walk_params.json.

Per-trial rows -> v5_results.csv (trial, kine_score, duty, knee_min,
rmse_hip/knee/ankle, cadence, tilt_max, stayed_up + gains + run type).
On each NEW global best: the full 22 s schedule runs, spinal_run.npz is
kept as v5_best_trial<N>.npz, and a run=full22 row is appended.

Usage: python optuna_walk_v5.py [n_trials]   (default 100)
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import optuna

import kine_ref
import params
import runner as R

HERE = Path(__file__).parent
DB = "sqlite:///optuna_walk.db"
STUDY = "ground_walk_v5_phase"
CSV = HERE / "v5_results.csv"
BEST5 = HERE / "best_walk_params_v5.json"

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
FULL_SCHED = dict(params.SCHEDULE)      # runner --eval mutates SCHEDULE
TRIAL_STATE: dict[int, dict] = {}       # trial number -> metrics + params
BEST_SO_FAR = {"value": float("-inf"), "trial": -1}
CSV_FIELDS = ["ts", "trial", "phase_reset_e", "phase_reset_f", "run",
              "kine_score", "duty", "knee_min", "rmse_hip", "rmse_knee",
              "rmse_ankle", "cadence", "tilt_max", "stayed_up", "npz"]


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
    """One trial's parameters into the params module (in-place, same-dict
    binding as optuna_walk.set_params + the two v5 gains)."""
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


def csv_row(row: dict):
    new = not CSV.exists()
    with open(CSV, "a", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)
    print("csv: " + " ".join(f"{k}={row.get(k)}" for k in
                             ("trial", "run", "kine_score", "duty",
                              "knee_min", "cadence", "stayed_up")),
          flush=True)


def metrics_row(trial_no: int, p: dict, m: dict, run: str,
                npz: str = "") -> dict:
    k = m.get("kine") or {}
    walk_dur = (params.SCHEDULE["walk"][1] - params.SCHEDULE["walk"][0])
    return dict(ts=time.strftime("%Y-%m-%d %H:%M:%S"), trial=trial_no,
                phase_reset_e=f"{p.get('phase_reset_e', 0.0):.4f}",
                phase_reset_f=f"{p.get('phase_reset_f', 0.0):.4f}",
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
    """Full 22 s schedule at the new global best; keep the npz + CSV row."""
    params.SCHEDULE.update(FULL_SCHED)
    set_params(p)
    # walk_drive is runner.main's LOCAL (default 4.0) - the trial's drive
    # must be passed explicitly, same as the eval call
    R.main(["--drive", repr(p["drive"])])   # writes spinal_run.npz (+png)
    npz = HERE / f"v5_best_trial{trial_no}.npz"
    if (HERE / "spinal_run.npz").exists():
        npz.unlink(missing_ok=True)
        (HERE / "spinal_run.npz").rename(npz)
    # metrics from the saved run (R.main full mode returns None)
    z = np.load(npz, allow_pickle=True)
    t, q_deg, neuro = z["t"], z["q"], z["neuro"]
    com = z["com"]
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
    )
    set_params(p)
    # repr() = shortest round-trip float string: the trial is evaluated at
    # EXACTLY the drive that set_params stored, so best_walk_params_v5.json
    # + runner --best5 reproduces the winning eval BIT-EXACTLY. (v4 lesson:
    # the old f"{drive:.4f}" rounding made every --best reproduction differ
    # from its study score - v4b -61.174 study vs -61.818 --best - on a
    # chaotic sim a 3e-5 nA drive change moves kine_score by ~0.6.)
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
            params.SCHEDULE.update(FULL_SCHED)   # next --eval re-sets anyway


def main(argv):
    n_trials = int(argv[0]) if argv and argv[0].isdigit() else 100
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    load_fitted_baseline()
    study = optuna.create_study(
        direction="maximize", storage=DB, study_name=STUDY,
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(seed=11, n_startup_trials=10))
    if len(study.trials) == 0:
        # enqueue only when the study is truly empty (queued trials persist
        # in sqlite - re-enqueueing on a crashed relaunch would duplicate
        # the seed)
        # seed: the v4b winner in v5 coordinates (its multipliers + the
        # two phase-reset gains at 0 = the regression-gated baseline)
        seed = dict(drive=2.2, rg_adapt=1.9, desc_e=1.7, desc_f=1.4,
                    rg_to_pf=2.4, pf_gain=1.0, e2_pf=1.0, f1_df=1.0,
                    f1_kf=1.0, e2_adapt=1.0, post_kneext=1.0, post_hipext=1.0,
                    kx=150.0, phase_reset_e=0.0, phase_reset_f=0.0)
        try:
            prev = json.loads(
                (HERE / "best_walk_params.json").read_text("utf-8"))
            if prev.get("study") == "ground_walk_v4b_kine":
                mul = prev.get("multipliers", {})
                seed.update({k: float(mul[k]) for k in seed if k in mul})
                seed["phase_reset_e"] = 0.0
                seed["phase_reset_f"] = 0.0
                print("seeded v4b winner into v5 (gains 0)", flush=True)
        except (FileNotFoundError, KeyError):
            pass
        study.enqueue_trial(seed)
        print("enqueued seed trial", flush=True)
    # resume bookkeeping: never re-capture a stored global best (only
    # COMPLETED trials have values; a freshly enqueued seed would make
    # study.best_value raise "Record does not exist")
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
               phase_reset_f=best.params["phase_reset_f"])
    with open(BEST5, "w", encoding="utf-8") as f:
        json.dump({"score": best.value, "params": eff,
                   "pf_gain": g, "multipliers": best.params,
                   "baseline": dict(BASE), "trial": best.number,
                   "study": STUDY}, f, indent=2)
    print("saved best_walk_params_v5.json (runner --best5 reproduces; "
          "runner --fitted --best5 for the full chain)", flush=True)


if __name__ == "__main__":
    main(sys.argv[1:])
