"""SPIKING-MIRROR curriculum (2026-10-03 overnight campaign, Ben's ask:
"I want the whole curriculum for non-spiking and spiking variants after
the networks are verified").

Same staged ladder as _curriculum.py (1 air-deaff -> 2 air-aff ->
3 balance -> 4 ground walk -> 5 pf-variant), run on the AARL_NET=spiking
hybrid mirror (build_network_spiking.py, goal-1 verified 2026-10-02:
topology mirror exact, basin 12/12, untuned gait -236.35 BILATERAL).

Differences vs _curriculum.py (deliberate, all local to this file):
  * AARL_NET=spiking + a private AARL_NPZ pinned here (the s3k stage-1/2
    objective reads the hardcoded spinal_run.npz; this module reads ITS
    OWN npz and indexes columns BY NAME - the w2lvar convention).
  * air stages run --time 20 (window 5..19 s): the spiking mirror's
    untuned period is ~1.5-10x the analog's (goal-1 gate 2/4), so the
    rhythm gate needs the longer window to see >= 3 bursts.
  * DB = optuna_spiking.db, fresh dated study names curr_spk_s{S}_20261003.
  * stage-3 balance + stage-4/5 ground objectives reuse the exact
    _curriculum.py formulas (R.main metrics only - no npz coupling).

Usage: python _curriculum_spiking.py <stage> [n_trials]
"""
import json
import os
import sys

VARIANT = "spiking"
DB = "sqlite:///optuna_spiking.db"
NPZ = f"spinal_run_curr_spiking_{os.getpid()}.npz"
os.environ["AARL_NET"] = VARIANT
os.environ["AARL_NPZ"] = NPZ

import numpy as np  # noqa: E402
import optuna  # noqa: E402

import optuna_walk_v10 as OW  # noqa: E402
import params as P  # noqa: E402
import runner as R  # noqa: E402
# _curriculum wraps stdout at import (utf-8 console) - import it FIRST,
# then use its wrapper; do NOT add a second TextIOWrapper around the
# same buffer (GC closes the buffer - measured in the BPA spiking smoke)
import _curriculum as CS3  # noqa: E402  (set_stage reuse, read-only)

BASE_MUL = None

KEYS1 = ("drive", "rg_nap_h", "desc_e", "desc_f", "rg_to_pf")
KEYS2 = KEYS1 + ("heel_rge", "toe_rge", "ib_e_central", "ia_f_central",
                 "ii_f_central", "ii_e_central", "c1_gain", "v3_gain")
KEYS3 = ("vest_ext", "vest_flex_inh", "vest_prop", "rig_scale")
STAGE_KEYS = {1: KEYS1, 2: KEYS2, 3: KEYS3}
STUDY_NAMES = {1: "curr_spk_s1_air_deaff_20261003",
               2: "curr_spk_s2_air_aff_20261003",
               3: "curr_spk_s3_balance_20261003",
               4: "curr_spk_s4_ground_20261003"}
DEFAULT_N = {1: 25, 2: 25, 3: 25, 4: 25}

SEED_DEFAULTS = {
    1: dict(drive=2.5, rg_nap_h=0.35, desc_e=1.7, desc_f=1.4,
            rg_to_pf=2.4),
    2: dict(drive=2.5, rg_nap_h=0.35, desc_e=1.7, desc_f=1.4,
            rg_to_pf=2.4, heel_rge=0.0, toe_rge=0.0, ib_e_central=0.0,
            ia_f_central=0.0, ii_f_central=0.0, ii_e_central=0.0,
            c1_gain=0.1, v3_gain=0.0),
    3: dict(vest_ext=0.0, vest_flex_inh=0.0, vest_prop=0.0,
            rig_scale=1.0),
}


def _air_eval(p, deaff):
    """Shared stage-1/2 air eval. Returns the optuna score."""
    args = (["--no-ground", "--no-afferents", "--no-interleg"]
            if deaff else ["--no-ground"])
    args += ["--time", "20", "--drive", repr(float(p["drive"]))]
    R.main(args)
    z = np.load(NPZ, allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    names = [str(x) for x in z["neuro_names"]]
    joints = [str(x) for x in z["key_joints"]]
    i_rge = names.index("RG_E_r")
    i_knee = joints.index("knee_angle_r")
    w = (t >= 5.0) & (t <= 19.0)
    if not np.all(np.isfinite(q[w])) or not np.all(np.isfinite(neuro[w])):
        return -200.0
    knee = q[w, i_knee]
    if not (-360.0 < float(knee.min()) < 360.0) or \
            not (-360.0 < float(knee.max()) < 360.0):
        return -200.0
    rge = neuro[w, i_rge]
    on = rge > 0.5 * max(rge.max(), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    if rises > 40:
        return -200.0
    if rises < 3 or (float(rge.max()) - float(rge.min())) < 1.0:
        return -10.0 + 0.05 * (-float(knee.min()))
    score = 3.0 * rises + 0.5 * (-float(knee.min()))
    return float(score) if np.isfinite(score) else -200.0


def objective(stage):
    def obj(trial):
        sug = dict(
            drive=trial.suggest_float("drive", 1.2, 4.0),
            rg_nap_h=trial.suggest_float("rg_nap_h", 0.15, 0.90),
            desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
            desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
            rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        )
        if stage == 2:
            sug["heel_rge"] = trial.suggest_float("heel_rge", 0.0, 0.5)
            sug["toe_rge"] = trial.suggest_float("toe_rge", 0.0, 0.5)
            sug["ib_e_central"] = trial.suggest_float("ib_e_central",
                                                      0.0, 0.5)
            sug["ia_f_central"] = trial.suggest_float("ia_f_central",
                                                      0.0, 0.5)
            sug["ii_f_central"] = trial.suggest_float("ii_f_central",
                                                      0.0, 0.5)
            sug["ii_e_central"] = trial.suggest_float("ii_e_central",
                                                      0.0, 0.5)
            sug["c1_gain"] = trial.suggest_float("c1_gain", 0.1, 0.8)
            sug["v3_gain"] = trial.suggest_float("v3_gain", 0.0, 0.3)
        if stage == 3:
            sug["vest_ext"] = trial.suggest_float("vest_ext", 0.0, 1.0)
            sug["vest_flex_inh"] = trial.suggest_float("vest_flex_inh",
                                                       0.0, 1.0)
            sug["vest_prop"] = trial.suggest_float("vest_prop", 0.0, 1.0)
            sug["rig_scale"] = trial.suggest_float("rig_scale", 0.05,
                                                   1.0)
        p = {**BASE_MUL, **sug}
        CS3.set_stage(stage, p)
        if stage <= 2:
            return _air_eval(p, deaff=(stage == 1))
        if stage == 3:
            m = R.main(["--stand-eval", "8",
                        "--rig-scale", repr(float(p["rig_scale"]))])
            if m["nan"]:
                return -200.0
            if m["bal_fell"]:
                return -150.0
            return (100.0
                    - 400.0 * float(m["bal_sway"])
                    - 1.0 * float(m["bal_tilt_max"])
                    - 40.0 * abs(0.5 - float(m["bal_contact_sym"])))
        # stage 4: ground walk on the s3k-lineage objective (metrics
        # only; BASE pinned at the s3k winner like _curriculum stage 4)
        m = R.main(["--eval", "--drive", repr(p["drive"])])
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
    return obj


def main():
    global BASE_MUL
    stage = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_N[stage]
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5
    if stage >= 4:
        s3k = json.loads(open(
            "reports_20260923/s3k_trial34_full_params.json",
            encoding="utf-8").read())["params"]
        BASE_MUL.update({k: v for k, v in s3k.items()})
        if s3k.get("pf_gain") is None:
            BASE_MUL.pop("pf_gain", None)
    print(f"[spiking] pinned AARL_NET={os.environ['AARL_NET']} "
          f"AARL_NPZ={os.environ['AARL_NPZ']} DB={DB}", flush=True)
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", storage=DB,
                                study_name=STUDY_NAMES[stage],
                                load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(
                                    seed=21 + stage, n_startup_trials=8))
    if len(study.trials) == 0:
        sk = STAGE_KEYS[stage]
        seed = dict(SEED_DEFAULTS[stage])
        try:
            # warm start: this chain's own previous stage if present,
            # else the ANALOG chain's same-stage winner (good rhythm
            # warm start; the mirror maps the same knobs)
            prevw = json.loads(open(f"curriculum_spiking_stage{stage}"
                                    f".json", encoding="utf-8")
                               .read())["params"]
        except FileNotFoundError:
            prevw = {}
        if not prevw:
            for alt in (f"curriculum_spiking_stage{stage-1}.json",
                        f"curriculum_stage{min(stage, 2)}.json"):
                try:
                    prevw = json.loads(open(alt, encoding="utf-8")
                                       .read())["params"]
                    print(f"[spiking] warm start from {alt}", flush=True)
                    break
                except FileNotFoundError:
                    continue
        for k in sk:
            if k in prevw:
                seed[k] = float(prevw[k])
        study.enqueue_trial(seed)
        print("seeded", json.dumps(seed), flush=True)
    study.optimize(objective(stage), n_trials=n, gc_after_trial=True)
    best = study.best_trial
    print(f"== spiking stage {stage} best {best.value:.3f} "
          f"(trial {best.number})")
    print(json.dumps(best.params, indent=1))
    with open(f"curriculum_spiking_stage{stage}.json", "w",
              encoding="utf-8") as f:
        json.dump({"stage": stage, "score": best.value,
                   "params": best.params, "trial": best.number,
                   "study": STUDY_NAMES[stage], "db": DB, "npz": NPZ,
                   "aarl_net": VARIANT}, f, indent=2)
    print(f"saved curriculum_spiking_stage{stage}.json")


if __name__ == "__main__":
    main()
