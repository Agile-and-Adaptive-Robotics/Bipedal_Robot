"""w2lvar VARIANT curriculum chain (goal4, 2026-09-25) - Ben's 5 stages:
    1 deafferented air, 2 afferented air, 3 standing balance,
    4 WALKING WITHOUT CONTACT, 5 WALKING WITH CONTACT.
Adapted from _curriculum.py (stock stack); see goal4_wiring.md.

Variant-specific decisions (each verified this session):
- The net is selected by env AARL_NET=w2lvar (build_network.py:1026);
  AARL_NPZ=spinal_run_w2lvar.npz redirects the run output. Both are
  pinned inside main() so the command line is self-contained (one
  detached command per stage) and the concurrent syn6 chain never
  clobbers this chain's npz or the protected spinal_run.npz.
- OWN database sqlite:///optuna_w2lvar.db; optuna_walk.db untouched.
- FRESH study names curr_w2lvar_s1_air_deaff .. curr_w2lvar_s5_walk_contact.
- ONLY knobs the variant consumes are searched. build_network_w2lvar.py
  applies its VARIANT_G overlay at build() time (heel_rge/toe_rge/
  ib_rge/joint_pf/ia_in/renshaw/full_rules/f1_kneext_inh/
  f1_anklepf_inh), so curriculum values for those keys would be
  clobbered -> not searched. vest_ext/vest_flex_inh excluded (the
  variant builds no VEST cells, net.vest=False). c1_gain/v3_gain and
  the *_central keys are stock-builder knobs the variant builder never
  reads. Searched set: stage 1 rhythm core (read by the builder) +
  per-muscle reflex motif gains ia_to_mn/ia_to_antagonist/ii_to_mn/
  ib_to_mn_inh (builder-read, NOT overlaid; ranges bracket their
  params.py defaults 0.6/0.4/0.4/0.35 - pre-existing keys, not new
  default-0 gains, so they seed at the defaults) + rhythm/laying knobs
  rg_mutual_inh/rg_weak_exc/pf_to_mn (stage >= 4) + contact knobs
  contact_onset/contra_swing/ib_group_exc/ib_exc_to_mn (stage 5).
- Stage 4 = --no-ground + the WALKING kine objective vs the reference
  cycle (pattern-match the gait in air before facing contact). The kz
  COM-floor penalty is stage-5 only (in suspended air the rig holds
  the pelvis up, so kz carries no information there).
- The air objective derives its columns from the npz's OWN neuro_names
  / key_joints arrays (RG_E_r / knee_angle_r) - s3k ordering is NOT
  assumed; the indices are printed once per process.
- renshaw: OW.set_params hard-sets 0.5, which matches the variant's
  VARIANT_G overlay (rc motif wired at 0.5).
Stage chaining via curriculum_w2lvar_stage{N}.json; seeds are FULL
dicts (every searched key explicit - enqueue_trial samples missing
keys). Usage:
    python _curriculum_w2lvar.py <stage> [n_trials]   (defaults 18/18/18/22/22)
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

import optuna_walk_v10 as OW
import params as P
import runner as R

VARIANT = "w2lvar"
DB = f"sqlite:///optuna_{VARIANT}.db"
NPZ = f"spinal_run_{VARIANT}.npz"
BASE_MUL = None  # v10 winner multipliers; loaded in main, merged per trial

# searched keys per stage (single source of truth for the objective, the
# seed builder, and the prev-winner filter). Strictly additive S2 < S4 <
# S5; stage 3 stands alone. NO mid-range-default elif chain: seed values
# come from SEED_DEFAULTS below (the stock (1,2) leakage trap).
KEYS1 = ("drive", "rg_nap_h", "desc_e", "desc_f", "rg_to_pf")
KEYS2 = KEYS1 + ("ia_to_mn", "ia_to_antagonist", "ii_to_mn",
                 "ib_to_mn_inh")
KEYS4 = KEYS2 + ("rg_mutual_inh", "rg_weak_exc", "pf_to_mn")
KEYS5 = KEYS4 + ("contact_onset", "contra_swing", "ib_group_exc",
                 "ib_exc_to_mn")
KEYS3 = ("vest_prop", "rig_scale")
STAGE_KEYS = {1: KEYS1, 2: KEYS2, 3: KEYS3, 4: KEYS4, 5: KEYS5}
DEFAULT_N = {1: 18, 2: 18, 3: 18, 4: 22, 5: 22}

# FULL-dict seed values per stage = the params.py defaults the gate-c
# rhythm smoke ran at (drive 2.5 / nap_h 0.35 / desc 1.7/1.4 / rg_to_pf
# 2.4), NOT v10 multipliers (those tuned the stock net). Keys absent
# from the previous stage json keep these values (trial 0 of each stage
# = the previous stage reproduced + new pathways at defaults).
_SEED1 = dict(drive=2.5, rg_nap_h=0.35, desc_e=1.7, desc_f=1.4,
              rg_to_pf=2.4)
SEED_DEFAULTS = {
    1: dict(_SEED1),
    2: {**_SEED1, "ia_to_mn": 0.6, "ia_to_antagonist": 0.4,
        "ii_to_mn": 0.4, "ib_to_mn_inh": 0.35},
    3: {"vest_prop": 0.0, "rig_scale": 1.0},
    4: {**_SEED1, "ia_to_mn": 0.6, "ia_to_antagonist": 0.4,
        "ii_to_mn": 0.4, "ib_to_mn_inh": 0.35,
        "rg_mutual_inh": 4.0, "rg_weak_exc": 0.0, "pf_to_mn": 2.0},
    5: {**_SEED1, "ia_to_mn": 0.6, "ia_to_antagonist": 0.4,
        "ii_to_mn": 0.4, "ib_to_mn_inh": 0.35,
        "rg_mutual_inh": 4.0, "rg_weak_exc": 0.0, "pf_to_mn": 2.0,
        "contact_onset": 0.0, "contra_swing": 0.0,
        "ib_group_exc": 0.5, "ib_exc_to_mn": 0.6},
}

STUDY_NAMES = {1: "curr_w2lvar_s1_air_deaff",
               2: "curr_w2lvar_s2_air_aff",
               3: "curr_w2lvar_s3_balance",
               4: "curr_w2lvar_s4_walk_nocontact",
               5: "curr_w2lvar_s5_walk_contact"}


def set_stage(stage, p):
    OW.load_fitted_baseline()
    OW.set_params(p)          # desc_e/desc_f -> descend_to_rg_e/f etc.;
                              # hard-sets renshaw 0.5 (matches VARIANT_G)
    P.TAU["rg_nap_h"] = float(p.get("rg_nap_h", 0.35))
    if stage >= 2:
        P.G["ia_to_mn"] = float(p.get("ia_to_mn", 0.6))
        P.G["ia_to_antagonist"] = float(p.get("ia_to_antagonist", 0.4))
        P.G["ii_to_mn"] = float(p.get("ii_to_mn", 0.4))
        P.G["ib_to_mn_inh"] = float(p.get("ib_to_mn_inh", 0.35))
    if stage >= 4:
        P.G["rg_mutual_inh"] = float(p.get("rg_mutual_inh", 4.0))
        P.G["rg_weak_exc"] = float(p.get("rg_weak_exc", 0.0))
        P.G["pf_to_mn"] = float(p.get("pf_to_mn", 2.0))
        # stock-stack rig/height env knobs are NOT searched here - make
        # sure a stale value from the parent shell cannot leak in
        os.environ.pop("AARL_KY", None)
        os.environ.pop("AARL_PELVIS_TY", None)
    if stage >= 5:
        P.G["contact_onset"] = float(p.get("contact_onset", 0.0))
        P.G["contra_swing"] = float(p.get("contra_swing", 0.0))
        P.G["ib_group_exc"] = float(p.get("ib_group_exc", 0.5))
        P.G["ib_exc_to_mn"] = float(p.get("ib_exc_to_mn", 0.6))
    if stage == 3:
        # vest_prop is runner-side (stance-gated II loop boost);
        # vest_ext/vest_flex_inh need VEST cells the variant lacks
        P.G["vest_prop"] = float(p.get("vest_prop", 0.0))


def objective(stage):
    def obj(trial):
        sug = dict(
            drive=trial.suggest_float("drive", 1.2, 4.0),
            rg_nap_h=trial.suggest_float("rg_nap_h", 0.15, 0.90),
            desc_e=trial.suggest_float("desc_e", 0.8, 1.8),
            desc_f=trial.suggest_float("desc_f", 0.7, 2.2),
            rg_to_pf=trial.suggest_float("rg_to_pf", 1.8, 3.0),
        )
        if stage >= 2:
            # per-muscle reflex motifs (builder-read G keys; ranges
            # bracket the params defaults they seed at)
            sug["ia_to_mn"] = trial.suggest_float("ia_to_mn", 0.0, 0.9)
            sug["ia_to_antagonist"] = trial.suggest_float(
                "ia_to_antagonist", 0.0, 0.6)
            sug["ii_to_mn"] = trial.suggest_float("ii_to_mn", 0.0, 0.6)
            sug["ib_to_mn_inh"] = trial.suggest_float("ib_to_mn_inh",
                                                      0.0, 0.5)
        if stage >= 4:
            sug["rg_mutual_inh"] = trial.suggest_float("rg_mutual_inh",
                                                       2.0, 6.0)
            # conditional topology: 0 = RG<->RG weak-excitation edges
            # absent (lab law - the seed keeps the current build)
            sug["rg_weak_exc"] = trial.suggest_float("rg_weak_exc",
                                                     0.0, 0.5)
            sug["pf_to_mn"] = trial.suggest_float("pf_to_mn", 0.5, 3.0)
        if stage >= 5:
            sug["contact_onset"] = trial.suggest_float("contact_onset",
                                                       0.0, 1.0)
            sug["contra_swing"] = trial.suggest_float("contra_swing",
                                                      0.0, 1.5)
            sug["ib_group_exc"] = trial.suggest_float("ib_group_exc",
                                                      0.0, 0.8)
            sug["ib_exc_to_mn"] = trial.suggest_float("ib_exc_to_mn",
                                                      0.0, 0.9)
        if stage == 3:
            sug["vest_prop"] = trial.suggest_float("vest_prop", 0.0, 1.0)
            sug["rig_scale"] = trial.suggest_float("rig_scale", 0.05,
                                                   1.0)
        # searched keys override; everything else pinned at v10 winner
        p = {**BASE_MUL, **sug}
        set_stage(stage, p)
        if stage in (1, 2):
            # AIR stepping (stage 1 deafferented + interleg off; stage 2
            # AFFERENTED + interleg ON)
            args = ["--no-ground", "--no-afferents", "--no-interleg",
                    "--time", "14", "--drive", repr(p["drive"])]
            if stage == 2:
                args = ["--no-ground", "--time", "14",
                        "--drive", repr(p["drive"])]
            m = R.main(args)
            import numpy as np
            z = np.load(NPZ, allow_pickle=True)
            t, q, neuro = z["t"], z["q"], z["neuro"]
            # column identity FROM THE NPZ (no assumed s3k ordering)
            names = [str(x) for x in z["neuro_names"]]
            joints = [str(x) for x in z["key_joints"]]
            i_rge, i_knee = names.index("RG_E_r"), \
                joints.index("knee_angle_r")
            if not hasattr(obj, "_cols_seen"):
                obj._cols_seen = True
                print(f"[w2lvar] npz columns verified: neuro_names="
                      f"{names} -> RG_E_r idx {i_rge}; key_joints -> "
                      f"knee_angle_r idx {i_knee}", flush=True)
            m = (t >= 5.0) & (t <= 17.0)
            if not np.all(np.isfinite(q[m])) or \
                    not np.all(np.isfinite(neuro[m])):
                return -200.0
            knee = q[m, i_knee]
            if not (-360.0 < float(knee.min()) < 360.0) or \
                    not (-360.0 < float(knee.max()) < 360.0):
                return -200.0  # unphysical RoM (finite but exploded)
            rge = neuro[m, i_rge]
            on = rge > 0.5 * max(rge.max(), 1e-9)
            rises = int(np.sum(np.diff(on.astype(int)) == 1))
            knee = q[m, i_knee]
            if rises > 30:
                return -200.0  # runaway flutter, not stepping
            # RHYTHM GATE (stock 2026-09-20): a static deep-flexion pose
            # must order BELOW every genuine rhythm (>=3 bursts)
            if rises < 3 or (float(rge.max()) - float(rge.min())) < 1.0:
                return -10.0 + 0.05 * (-float(knee.min()))
            # air objective: rhythmic + deep knee swing flexion
            score = 3.0 * rises + 0.5 * (-float(knee.min()))
            if not np.isfinite(score):
                return -200.0
            return float(score)
        if stage == 3:
            # STANDING-BALANCE eval (SCONE Tutorial-3a analog): 8 s
            # standing (--stand-eval 8, DRIVE 0) at this trial's rig
            # wean. Sentinels: NaN -200 < fall -150 < any stander.
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
        # stages 4-5: the WALKING kine objective. Stage 4 runs the model
        # with contact DISABLED (--no-ground) but scores the same gait
        # pattern (documented interpretation: pattern-match the gait in
        # air before facing contact); stage 5 is the normal ground eval.
        args = ["--eval", "--drive", repr(p["drive"])]
        if stage == 4:
            args = ["--no-ground", "--eval", "--drive", repr(p["drive"])]
        m = R.main(args)
        if m["nan"]:
            return -400.0
        if m.get("kine") is None:
            # no cycles on either leg - MUST sit below the worst genuine
            # walker (stock -320 discipline)
            return -320.0
        # clip so pathological walkers stay above the sentinels
        score = max(float(m["kine_score"]), -315.0)
        if stage == 5 and m["kz"] < 0.62:
            # COM floor: meaningful on the ground only (in air the rig
            # holds the pelvis up)
            score -= 20.0
        if m["tilt_max"] > 40.0:
            score -= 10.0
        return score
    return obj


def main():
    global BASE_MUL
    # SELF-CONTAINED pin: variant net + chain-private npz. The launch
    # step runs this script as one detached command; no env setup needed.
    os.environ["AARL_NET"] = VARIANT
    os.environ["AARL_NPZ"] = NPZ
    print(f"[w2lvar] pinned AARL_NET={os.environ['AARL_NET']} "
          f"AARL_NPZ={os.environ['AARL_NPZ']} DB={DB}", flush=True)
    stage = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    n = int(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_N[stage]
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5
    # NO s3k merge (the stock stage-4/5 move): s3k tuned the STOCK
    # topology; this variant chain is self-contained from its own
    # stage winners + the v10 multiplier shell that OW.set_params needs.
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", storage=DB,
                                study_name=STUDY_NAMES[stage],
                                load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(
                                    seed=21 + stage, n_startup_trials=8))
    if len(study.trials) == 0:
        # FULL-dict seed: every searched key explicit (missing keys get
        # SAMPLED by optuna - the JSON-rule lesson). Chain from this
        # variant's own stage jsons when they exist.
        sk = STAGE_KEYS[stage]
        seed = dict(SEED_DEFAULTS[stage])
        try:
            prevw = json.loads(open(f"curriculum_{VARIANT}_"
                                    f"stage{stage}.json",
                                    encoding="utf-8").read())["params"]
        except FileNotFoundError:
            try:
                prevw = json.loads(open(
                    f"curriculum_{VARIANT}_stage{stage-1}.json",
                    encoding="utf-8").read())["params"]
            except FileNotFoundError:
                prevw = {}
        for k in sk:
            if k in prevw:
                seed[k] = float(prevw[k])
        study.enqueue_trial(seed)
        print("seeded", json.dumps(seed), flush=True)
    study.optimize(objective(stage), n_trials=n, gc_after_trial=True)
    best = study.best_trial
    print(f"== stage {stage} best {best.value:.3f} (trial {best.number})")
    print(json.dumps(best.params, indent=1))
    with open(f"curriculum_{VARIANT}_stage{stage}.json", "w",
              encoding="utf-8") as f:
        json.dump({"stage": stage, "score": best.value,
                   "params": best.params, "trial": best.number,
                   "study": STUDY_NAMES[stage], "db": DB, "npz": NPZ,
                   "aarl_net": VARIANT}, f, indent=2)
    print(f"saved curriculum_{VARIANT}_stage{stage}.json")


if __name__ == "__main__":
    main()
