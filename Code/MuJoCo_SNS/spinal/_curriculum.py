"""Staged tuning curriculum (Ben 2026-09-15): deafferented-air ->
afferented-air -> supported ground, each stage seeded by the previous,
only the NEW pathway gains searched per stage.

Stage 1 (air_deaff): deafferented air-stepping. Tunes: drive/rg_adapt/
    desc_e/rg_to_pf (rhythm core).
Stage 2 (air_aff): afferented air (heel/toe + load signals active, feet
    free). Adds: phase_reset_e/f, heel_rge, toe_rge.
Stage 3 (ground): supported ground walk. Adds: ib_rge (stance-Ib
    prolonger), ia_in (IaIN pathway), ankle_post_walk_trim.

Usage: python _curriculum.py <stage> [n_trials]
Each stage writes curriculum_stage<N>.json (its best) so the next stage
seeds from it. Resumable per-stage via distinct optuna study names.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

import optuna_walk_v10 as OW
import params as P
import runner as R

DB = "sqlite:///optuna_walk.db"
BASE_MUL = None  # v10 winner multipliers; loaded in main, merged per trial


def set_stage(stage, p):
    OW.load_fitted_baseline()
    OW.set_params(p)
    P.G["renshaw"] = 0.5
    P.G["ankle_post_walk_trim"] = float(p.get("ankle_post_walk_trim", 1.0))
    P.TAU["rg_nap_h"] = float(p.get("rg_nap_h", 0.35))
    if stage >= 2:
        P.G["heel_rge"] = float(p.get("heel_rge", 0.0))
        P.G["toe_rge"] = float(p.get("toe_rge", 0.0))
        P.G["ib_e_central"] = float(p.get("ib_e_central", 0.0))
        P.G["ia_f_central"] = float(p.get("ia_f_central", 0.0))
        P.G["ii_f_central"] = float(p.get("ii_f_central", 0.0))
        P.G["ii_e_central"] = float(p.get("ii_e_central", 0.0))
    if stage >= 3:
        P.G["ib_rge"] = float(p.get("ib_rge", 0.0))
        P.G["ia_in"] = float(p.get("ia_in", 0.0))
        P.G["ia_f_contra_f"] = float(p.get("ia_f_contra_f", 0.0))
        P.G["v3_to_ibexc"] = float(p.get("v3_to_ibexc", 0.0))


def objective(stage):
    def obj(trial):
        sug = dict(
            drive=trial.suggest_float("drive", 1.2, 3.2),
            rg_nap_h=trial.suggest_float("rg_nap_h", 0.15, 0.90),
            desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
            desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
            rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        )
        if stage >= 2:
            sug["heel_rge"] = trial.suggest_float("heel_rge", 0.0, 1.0)
            sug["toe_rge"] = trial.suggest_float("toe_rge", 0.0, 1.0)
            sug["ib_e_central"] = trial.suggest_float("ib_e_central",
                                                      0.0, 1.0)
            sug["ia_f_central"] = trial.suggest_float("ia_f_central",
                                                      0.0, 1.0)
            sug["ii_f_central"] = trial.suggest_float("ii_f_central",
                                                      0.0, 1.0)
            sug["ii_e_central"] = trial.suggest_float("ii_e_central",
                                                      0.0, 1.0)
        if stage >= 3:
            sug["ib_rge"] = trial.suggest_float("ib_rge", 0.0, 1.0)
            sug["ia_in"] = trial.suggest_float("ia_in", 0.0, 1.2)
            sug["ia_f_contra_f"] = trial.suggest_float("ia_f_contra_f",
                                                       0.0, 1.0)
            sug["v3_to_ibexc"] = trial.suggest_float("v3_to_ibexc",
                                                     0.0, 1.0)
            sug["ankle_post_walk_trim"] = trial.suggest_float(
                "ankle_post_walk_trim", 0.05, 1.0)
        # searched keys override; everything else pinned at v10 winner
        p = {**BASE_MUL, **sug}
        set_stage(stage, p)
        if stage == 1:
            m = R.main(["--no-ground", "--no-afferents", "--no-interleg",
                        "--time", "14", "--drive", repr(p["drive"])])
            import numpy as np
            z = __import__("numpy").load("spinal_run.npz", allow_pickle=True)
            t, q, neuro = z["t"], z["q"], z["neuro"]
            m = (t >= 5.0) & (t <= 17.0)
            if not np.all(np.isfinite(q[m])) or \
                    not np.all(np.isfinite(neuro[m])):
                return -200.0
            knee = q[m, 4]
            if not (-360.0 < float(knee.min()) < 360.0) or \
                    not (-360.0 < float(knee.max()) < 360.0):
                return -200.0  # unphysical RoM (finite but exploded)
            rge = neuro[m, 2]
            on = rge > 0.5 * max(rge.max(), 1e-9)
            rises = int(np.sum(np.diff(on.astype(int)) == 1))
            knee = q[m, 4]
            # air objective: rhythmic + deep knee swing flexion
            score = 3.0 * rises + 0.5 * (-float(knee.min()))
            if not np.isfinite(score) or rises > 30:
                return -200.0
            return float(score)
        m = R.main(["--eval", "--drive", repr(p["drive"])])
        if m["nan"]:
            return -110.0 + m["t_end"]
        if m.get("kine") is None:
            return -100.0
        score = float(m["kine_score"])
        if m["kz"] < 0.62:
            score -= 10.0
        if m["tilt_max"] > 40.0:
            score -= 5.0
        return score
    return obj


def main():
    global BASE_MUL
    stage = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    name = {1: "curr_s1_air_deaff", 2: "curr_s2_air_aff",
            3: "curr_s3_ground"}[stage]
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", storage=DB,
                                study_name=name, load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(
                                    seed=21 + stage, n_startup_trials=8))
    if len(study.trials) == 0:
        try:
            seed = json.loads(open(f"curriculum_stage{stage-1}.json",
                                   encoding="utf-8").read())["params"]
        except FileNotFoundError:
            prev = json.loads(open("best_walk_params_v10.json",
                                   encoding="utf-8").read())
            seed = dict(prev["multipliers"])
        seed = {k: v for k, v in seed.items()
                if k in ("drive", "rg_nap_h", "desc_e", "desc_f", "rg_to_pf",
                         "heel_rge", "toe_rge", "ib_e_central",
                         "ia_f_central", "ii_f_central", "ii_e_central",
                         "ib_rge", "ia_in", "ia_f_contra_f", "v3_to_ibexc",
                         "ankle_post_walk_trim")}
        study.enqueue_trial(seed)
        print("seeded", flush=True)
    study.optimize(objective(stage), n_trials=n, gc_after_trial=True)
    best = study.best_trial
    print(f"== stage {stage} best {best.value:.3f} (trial {best.number})")
    print(json.dumps(best.params, indent=1))
    with open(f"curriculum_stage{stage}.json", "w", encoding="utf-8") as f:
        json.dump({"stage": stage, "score": best.value,
                   "params": best.params, "trial": best.number,
                   "study": name}, f, indent=2)
    print(f"saved curriculum_stage{stage}.json")


if __name__ == "__main__":
    main()
