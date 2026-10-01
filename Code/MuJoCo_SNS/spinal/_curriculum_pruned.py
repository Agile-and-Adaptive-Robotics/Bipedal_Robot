"""PRUNED-s3k retune campaign (Ben 2026-09-30): the 4-way pruned s3k
(combo config from the 09-28 prune matrix: interleg + contact + Ib +
rg-weak cut) RETUNED - short pass then long batches with Ben's stop
rule (2 consecutive batches < +0.5), three parallel studies.

The prune overlay (COMBO keys zeroed + --no-interleg on EVERY R.main)
mirrors prune_combo.py exactly; searched keys are only what REMAINS
meaningful in the pruned architecture.

Subcommands:
    python _curriculum_pruned.py air              # 30 trials, air rhythm
    python _curriculum_pruned.py walk <idx>       # batches of 40, stop
                                                   # rule, max 5 batches;
                                                   # idx in 1..3 (seed)
    python _curriculum_pruned.py finalize         # merge studies, write
                                                   # winner json + val
Artifacts: optuna_pruned.db, curriculum_pruned_air.json,
curriculum_pruned_walk.json, pruned_val_<mode>.jsonl.
"""
import io
import json
import os
import sys

for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

DB = "sqlite:///optuna_pruned.db"
NPZ = f"spinal_run_pruned_{os.getpid()}.npz"
os.environ["AARL_NPZ"] = NPZ
os.environ.pop("AARL_NET", None)

import time as _time

import numpy as np
import optuna

import optuna_walk_v10 as OW
import params as P
import runner as R
import _curriculum as CS3

# ---- the prune overlay (prune_combo.py's COMBO_KEYS, verbatim) ----
PRUNE_KEYS = {"heel_rge": 0.0, "toe_rge": 0.0, "contact_onset": 0.0,
              "contra_swing": 0.0, "pm_gain": 0.0, "pm_aff": 0.0,
              "ib_rge": 0.0, "ib_e_central": 0.0, "rg_weak_exc": 0.0}
PRUNE_ARGV = ["--no-interleg"]

T34 = json.loads(open(
    "reports_20260923/s3k_trial34_full_params.json",
    encoding="utf-8").read())["params"]

# v10 winner multipliers - OW.set_params READS keys from here (rg_adapt
# etc.); trial 34 does not carry them, so every config merge needs BASE.
BASE = dict(json.loads(open(
    "best_walk_params_v10.json", encoding="utf-8").read())["multipliers"])
BASE["renshaw"] = 0.5

# pinned (not searched): ky_scale/pelvis_ty/pm_T/pm_add/pm_ws and the
# Ia/II central gains + full_rules/no_cross stay at trial 34; the prune
# overlay kills what trial 34 had that the combo removes.
PIN = {k: T34[k] for k in ("ky_scale", "pelvis_ty", "pm_T", "pm_add",
                           "pm_ws", "ia_f_central", "ii_f_central",
                           "ii_e_central", "full_rules", "no_cross",
                           "v3_to_ibexc", "c1_gain")}
# c1_gain/v3_gain/ia_f_contra_f/contra_kinh are overwritten by the
# no_cross branch of set_stage anyway; carried for completeness.

AIR_KEYS = ("drive", "rg_nap_h", "desc_e", "desc_f", "rg_to_pf")
AIR_SEED = {k: float(T34[k]) for k in AIR_KEYS}

WALK_KEYS = ("drive", "rg_nap_h", "desc_e", "desc_f", "rg_to_pf",
             "ia_in", "ankle_post_walk_trim", "f1_anklepf_inh",
             "pf_gain")
WALK_SEED = {k: float(T34[k]) for k in WALK_KEYS}


def _config(p):
    """Full-dict -> set_stage(5) -> prune overlay. Returns the argv
    prefix every R.main call must carry. BASE (v10 multipliers incl
    rg_adapt) must be in the merge - OW.set_params reads it."""
    full = {**BASE, **{k: float(T34[k]) for k in T34},
            **PIN, **{k: float(v) for k, v in p.items()}}
    CS3.set_stage(5, dict(full))
    for k, v in PRUNE_KEYS.items():
        P.G[k] = float(v)
    return list(PRUNE_ARGV)


def _air_obj(trial):
    sug = dict(
        drive=trial.suggest_float("drive", 1.2, 3.2),
        rg_nap_h=trial.suggest_float("rg_nap_h", 0.15, 0.90),
        desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
        desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
        rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
    )
    argv = _config(sug)
    m = R.main(["--no-ground", "--time", "14",
                "--drive", repr(float(sug["drive"]))] + argv)
    z = np.load(NPZ, allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    names = [str(x) for x in z["neuro_names"]]
    joints = [str(x) for x in z["key_joints"]]
    i_rge, i_knee = names.index("RG_E_r"), joints.index("knee_angle_r")
    w = (t >= 5.0) & (t <= 17.0)
    if not np.all(np.isfinite(q[w])) or not np.all(np.isfinite(neuro[w])):
        return -200.0
    knee = q[w, i_knee]
    if not (-360.0 < float(knee.min()) < 360.0) or \
            not (-360.0 < float(knee.max()) < 360.0):
        return -200.0
    rge = neuro[w, i_rge]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    if rises > 30:
        return -200.0
    if rises < 3 or (float(rge.max()) - float(rge.min())) < 1.0:
        return -10.0 + 0.05 * (-float(knee.min()))
    score = 3.0 * rises + 0.5 * (-float(knee.min()))
    return float(score) if np.isfinite(score) else -200.0


def _walk_obj(trial):
    sug = dict(
        drive=trial.suggest_float("drive", 1.6, 2.6),
        rg_nap_h=trial.suggest_float("rg_nap_h", 0.15, 0.90),
        desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
        desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
        rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        ia_in=trial.suggest_float("ia_in", 0.0, 1.2),
        ankle_post_walk_trim=trial.suggest_float(
            "ankle_post_walk_trim", 0.05, 1.0),
        f1_anklepf_inh=trial.suggest_float("f1_anklepf_inh", 0.0, 1.2),
        pf_gain=trial.suggest_float("pf_gain", 0.15, 3.0, log=True),
    )
    argv = _config(sug)
    m = R.main(["--eval", "--drive", repr(float(sug["drive"]))] + argv)
    if m["nan"]:
        return -400.0
    if m.get("kine") is None:
        return -320.0
    score = max(float(m["kine_score"]), -315.0)
    if m["kz"] < 0.62:
        score -= 20.0
    if m["tilt_max"] > 40.0:
        score -= 10.0
    return score


def _study(name, seed, obj, seed_params):
    # PER-STUDY DB FILE: 4 concurrent processes creating studies in ONE
    # sqlite file hit lock contention ("Record does not exist",
    # sqlalche.me/e3q8) - separate files per study are the fix.
    db = f"sqlite:///optuna_pruned_{name}.db"
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    st = optuna.create_study(direction="maximize", storage=db,
                             study_name=name, load_if_exists=True,
                             sampler=optuna.samplers.TPESampler(
                                 seed=seed, n_startup_trials=8))
    if len(st.trials) == 0:
        st.enqueue_trial(dict(seed_params))   # FULL dict: every key
        print("seeded", json.dumps(seed_params), flush=True)
    return st


def run_air():
    st = _study("pruned_air", 11, _air_obj, AIR_SEED)
    st.optimize(_air_obj, n_trials=30, gc_after_trial=True)
    b = st.best_trial
    print(f"== pruned air best {b.value:.3f} (trial {b.number})")
    json.dump({"stage": "air", "score": b.value, "params": b.params,
               "trial": b.number, "study": "pruned_air", "db": DB},
              open("curriculum_pruned_air.json", "w",
                   encoding="utf-8"), indent=2)


def run_walk(idx):
    """Batches of 40; stop when 2 consecutive batches gain < +0.5 on
    the running best (Ben's generations stop rule); hard cap 5 batches
    = 200 trials per study."""
    name = f"pruned_walk_{idx}"
    seed = 40 + idx
    st = _study(name, seed, _walk_obj, WALK_SEED)
    stall = 0
    # prev_best must start None: after enqueue the study holds only a
    # WAITING trial - study.best_value on it raises "Record does not
    # exist" (cost the second relaunch).
    prev_best = None
    for b in range(5):
        st.optimize(_walk_obj, n_trials=40, gc_after_trial=True)
        best = st.best_value
        gain = (best - prev_best) if prev_best is not None else float("inf")
        print(f"[{name}] batch {b + 1}: best {best:.3f} "
              f"(gain {gain:+.3f})", flush=True)
        if gain < 0.5:
            stall += 1
        else:
            stall = 0
        prev_best = best
        if stall >= 2:
            print(f"[{name}] stop rule (2 batches < +0.5)", flush=True)
            break


def finalize():
    studies = []
    for idx in (1, 2, 3):
        try:
            studies.append(optuna.load_study(
                study_name=f"pruned_walk_{idx}",
                storage=f"sqlite:///optuna_pruned_pruned_walk_{idx}.db"))
        except KeyError:
            pass
    best_t, best_s = None, None
    for st in studies:
        for t in st.trials:
            if t.state.name == "COMPLETE" and \
                    (best_t is None or t.value > best_t.value):
                best_t, best_s = t, st.study_name
    if best_t is None:
        print("no completed walk trials yet", flush=True)
        return
    doc = {"score": best_t.value, "params": best_t.params,
           "trial": best_t.number, "study": best_s, "db": DB,
           "merged_from": [st.study_name for st in studies],
           "pin": PIN, "prune_keys": PRUNE_KEYS,
           "prune_argv": PRUNE_ARGV}
    json.dump(doc, open("curriculum_pruned_walk.json", "w",
                        encoding="utf-8"), indent=2)
    print(f"winner: {best_s} t{best_t.number} {best_t.value:.3f}",
          flush=True)
    # validation cells with the winner
    out = "pruned_val.jsonl"
    argv = _config(best_t.params)
    cells = [("WALK", ["--eval", "--drive",
                       repr(float(best_t.params["drive"]))]),
             ("STAND", ["--stand-eval", "8", "--rig-scale", "1.0"])]
    cells += [(f"PUSH {ax}", ["--stand-eval", "8", "--rig-scale", "1.0",
                               "--push", "40", "--push-time", "4.0",
                               "--push-axis", ax])
              for ax in ("+x", "-x", "+y", "-y")]
    for tag, a in cells:
        t0 = _time.time()
        m = R.main(a + argv) or {}
        met = {k: v for k, v in m.items()
               if isinstance(v, (int, float, bool))}
        rec = {"config": "pruned-winner", "mode": tag,
               "metrics": met, "wall_s": round(_time.time() - t0, 1)}
        with open(out, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"[val] {tag}: kine={met.get('kine_score')} "
              f"fell={met.get('bal_fell')} "
              f"sway={met.get('bal_sway')}", flush=True)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "air"
    if cmd == "air":
        run_air()
    elif cmd == "walk":
        run_walk(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    elif cmd == "finalize":
        finalize()
