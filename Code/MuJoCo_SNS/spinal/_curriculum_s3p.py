"""Stage-3 PERTURBATION + METABOLIC curriculum (2026-10-20, Ben's
wave-2 expansion: 'the standing test must not be quiet standing').

Second stage-3 study per family (fresh names; the quiet curr_s3_balance
studies are NOT discarded): each trial runs a 7-eval standing battery -
quiet stand, 4 pelvis push pulses (forward/backward/lateral R/lateral L,
60 N = ~8.5% BW, 150 ms, at t=4 s), platform translate (6 cm +x) and
platform tilt (6 deg sagittal) - all 8 s, DRIVE 0, at the trial's rig
wean. Score = fall count + sway response to the PERTURBATION SET (not
just quiet sway) + a metabolic/co-contraction term (Ben: 'balance
EFFICIENTLY instead of just co-contracting everything'):

    nan anywhere                        -> -200   (sentinel)
    fell in ALL 7 battery evals         -> -150   (sentinel)
    else: 100 - 25*falls - 300*quiet_sway[m] - 1*quiet_tilt[deg]
             - 40*|0.5 - quiet_contact_sym|
             - 250*mean(post-push COM sway, 4 evals)
             - 250*mean(post-platform COM sway, 2 evals)
             - 30*bal_effort (mean squared muscle activation, quiet run)
             - 20*bal_cocon (antagonist co-contraction index, quiet run)
    clipped at -149 so partial standers order above the all-fall sentinel.

Searched = the family's quiet stage-3 keys PLUS the balance-mechanism
knobs from Ben's prior (wide stance, slight knee flexion, antagonist
co-activation): stance_width [0,10] deg hip abduction both hips, knee_bias
[0,15] deg extra knee flexion, coact [0,1] antagonist posture drive
(build_network), pelvis_ty [0.84,0.93] (a flexed-knee stand needs a lower
pelvis). New knobs default 0 = incumbent pose exactly (trial 0 = the
family's QUIET stage-3 winner reproduced under the perturbation battery).

Usage: python _curriculum_s3p.py <stock|w2lvar|syn6> [n_trials]
(defaults: stock 20, variants 18 - each trial = 7 evals of ~70 s).
"""
import io
import importlib
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna
import numpy as np

FAMILY = sys.argv[1] if len(sys.argv) > 1 else "stock"
N_DEF = {"stock": 20, "w2lvar": 18, "syn6": 18}
N = int(sys.argv[2]) if len(sys.argv) > 2 else N_DEF[FAMILY]

CFG = {
    "stock": dict(db="sqlite:///optuna_walk.db",
                  study="curr_s3p_perturb", mod="_curriculum",
                  quiet_json="curriculum_stage3.json",
                  npz="spinal_run_curr_s3p_20261020.npz", net=None),
    "w2lvar": dict(db="sqlite:///optuna_w2lvar.db",
                   study="curr_w2lvar_s3p_perturb",
                   mod="_curriculum_w2lvar",
                   quiet_json="curriculum_w2lvar_stage3.json",
                   npz="spinal_run_curr_w2lvar_s3p_20261020.npz",
                   net="w2lvar"),
    "syn6": dict(db="sqlite:///optuna_syn6.db",
                 study="curr_syn6_s3p_perturb", mod="_curriculum_syn6",
                 quiet_json="curriculum_syn6_stage3.json",
                 npz="spinal_run_curr_syn6_s3p_20261020.npz",
                 net="syn6"),
}[FAMILY]

os.environ["AARL_NPZ"] = CFG["npz"]
if CFG["net"]:
    os.environ["AARL_NET"] = CFG["net"]
else:
    os.environ.pop("AARL_NET", None)
print(f"[s3p:{FAMILY}] pinned AARL_NPZ={CFG['npz']} "
      f"AARL_NET={os.environ.get('AARL_NET', '<unset>')} "
      f"db={CFG['db']} study={CFG['study']}", flush=True)

C = importlib.import_module(CFG["mod"])
import params as P

# perturbation battery (all at event time 4 s, 8 s stands, DRIVE 0):
# 60 N ~ 8.5% of the ~72 kg model's body weight; 6 cm / 6 deg platform
# moves sit inside the literature ranges for support-surface protocols.
PUSH_N, PLAT_DX, PLAT_TILT, T_EVENT = 60.0, 0.06, 6.0, 4.0
BATTERY = [
    ("quiet", []),
    ("push_f", ["--push", repr(PUSH_N), "--push-axis", "+x"]),
    ("push_b", ["--push", repr(PUSH_N), "--push-axis", "-x"]),
    ("push_r", ["--push", repr(PUSH_N), "--push-axis", "+y"]),
    ("push_l", ["--push", repr(PUSH_N), "--push-axis", "-y"]),
    ("plat_dx", ["--plat-dx", repr(PLAT_DX)]),
    ("plat_tilt", ["--plat-tilt", repr(PLAT_TILT)]),
]

BASE_MUL = None


def load_params(p):
    """family set_stage + the s3p mechanism-knob loaders (JSON rule:
    every knob this study uses gets an explicit loader branch)."""
    C.set_stage(3, p)
    P.G["coact"] = float(p.get("coact", 0.0))
    if float(p.get("stance_width", 0.0)) > 0.0:
        os.environ["AARL_STANCE_W"] = repr(float(p["stance_width"]))
    else:
        os.environ.pop("AARL_STANCE_W", None)
    if float(p.get("knee_bias", 0.0)) > 0.0:
        os.environ["AARL_KNEE_BIAS"] = repr(float(p["knee_bias"]))
    else:
        os.environ.pop("AARL_KNEE_BIAS", None)
    if "pelvis_ty" in p:
        os.environ["AARL_PELVIS_TY"] = repr(float(p["pelvis_ty"]))
    else:
        os.environ.pop("AARL_PELVIS_TY", None)


def run_battery(rig_scale):
    """7 standing evals; returns list of metric dicts (battery order)."""
    out = []
    for _name, extra in BATTERY:
        args = (["--stand-eval", "8",
                 "--push-time", repr(T_EVENT),
                 "--rig-scale", repr(float(rig_scale))] + list(extra))
        m = C.R.main(args)
        out.append(m)
    return out


def score_battery(ms):
    falls = sum(1 for m in ms if m["bal_fell"])
    if any(m["nan"] for m in ms):
        return -200.0, falls
    if falls == len(ms):
        return -150.0, falls
    q = ms[0]

    def _ps(m):
        v = float(m["bal_push_sway"])
        return 0.30 if v < 0.0 else v   # -1 = degenerate: treat as bad

    push_sw = [_ps(m) for m, (nm, _x) in zip(ms, BATTERY)
               if nm.startswith("push")]
    plat_sw = [_ps(m) for m, (nm, _x) in zip(ms, BATTERY)
               if nm.startswith("plat")]
    s = (100.0
         - 25.0 * falls
         - 300.0 * float(q["bal_sway"])
         - 1.0 * float(q["bal_tilt_max"])
         - 40.0 * abs(0.5 - float(q["bal_contact_sym"]))
         - 250.0 * (sum(push_sw) / max(len(push_sw), 1))
         - 250.0 * (sum(plat_sw) / max(len(plat_sw), 1))
         - 30.0 * float(q["bal_effort"])
         - 20.0 * float(q["bal_cocon"]))
    return max(s, -149.0), falls


def main():
    global BASE_MUL
    prev = json.loads(open("best_walk_params_v10.json",
                           encoding="utf-8").read())
    BASE_MUL = dict(prev["multipliers"])
    BASE_MUL["renshaw"] = 0.5

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(direction="maximize", storage=CFG["db"],
                                study_name=CFG["study"],
                                load_if_exists=True,
                                sampler=optuna.samplers.TPESampler(
                                    seed=42, n_startup_trials=8))
    quiet_keys = ("vest_ext", "vest_flex_inh") if FAMILY == "stock" \
        else ("vest_prop",)
    sk = list(quiet_keys) + ["rig_scale", "stance_width", "knee_bias",
                             "coact", "pelvis_ty"]
    if len(study.trials) == 0:
        # FULL-dict seed: the family's QUIET stage-3 winner under the
        # perturbation battery (mechanism knobs 0 = incumbent pose).
        seed = {k: 0.0 for k in sk}
        seed["rig_scale"] = 1.0
        seed["pelvis_ty"] = 0.92
        try:
            qw = json.loads(open(CFG["quiet_json"],
                                 encoding="utf-8").read())["params"]
            for k in sk:
                if k in qw:
                    seed[k] = float(qw[k])
        except FileNotFoundError:
            print(f"[s3p:{FAMILY}] WARNING: {CFG['quiet_json']} missing; "
                  "seeding defaults", flush=True)
        study.enqueue_trial(seed)
        print(f"[s3p:{FAMILY}] seeded {json.dumps(seed)}", flush=True)

    def obj(trial):
        sug = {}
        if FAMILY == "stock":
            sug["vest_ext"] = trial.suggest_float("vest_ext", 0.0, 0.5)
            sug["vest_flex_inh"] = trial.suggest_float("vest_flex_inh",
                                                       0.0, 0.3)
        sug["vest_prop"] = trial.suggest_float("vest_prop", 0.0, 1.0)
        sug["rig_scale"] = trial.suggest_float("rig_scale", 0.05, 1.0)
        sug["stance_width"] = trial.suggest_float("stance_width",
                                                  0.0, 10.0)
        sug["knee_bias"] = trial.suggest_float("knee_bias", 0.0, 15.0)
        sug["coact"] = trial.suggest_float("coact", 0.0, 1.0)
        sug["pelvis_ty"] = trial.suggest_float("pelvis_ty", 0.84, 0.93)
        p = {**BASE_MUL, **sug}
        load_params(p)
        ms = run_battery(p["rig_scale"])
        val, falls = score_battery(ms)
        q = ms[0]
        trial.set_user_attr("falls", int(falls))
        trial.set_user_attr("quiet_sway", float(q["bal_sway"]))
        trial.set_user_attr("quiet_tilt", float(q["bal_tilt_max"]))
        trial.set_user_attr("quiet_effort", float(q["bal_effort"]))
        trial.set_user_attr("quiet_cocon", float(q["bal_cocon"]))
        trial.set_user_attr("push_sway_mean",
                            float(np.mean([m["bal_push_sway"]
                                           for m in ms[1:5]]))
                            if all(m["bal_push_sway"] >= 0
                                   for m in ms[1:5]) else -1.0)
        print(f"[s3p:{FAMILY}] t{trial.number} score {val:.2f} "
              f"falls {falls}/7 sway {float(q['bal_sway']):.3f} "
              f"tilt {float(q['bal_tilt_max']):.1f} "
              f"effort {float(q['bal_effort']):.3f} "
              f"cocon {float(q['bal_cocon']):.3f} "
              f"rig {p['rig_scale']:.2f}", flush=True)
        return val

    study.optimize(obj, n_trials=N, gc_after_trial=True)
    best = study.best_trial
    print(f"== [s3p:{FAMILY}] best {best.value:.3f} "
          f"(trial {best.number})")
    print(json.dumps(best.params, indent=1))
    with open(f"curriculum_s3p_{FAMILY}.json", "w", encoding="utf-8") as f:
        json.dump({"family": FAMILY, "stage": "3p", "score": best.value,
                   "params": best.params, "trial": best.number,
                   "study": CFG["study"], "db": CFG["db"],
                   "npz": CFG["npz"],
                   "aarl_net": CFG["net"],
                   "user_attrs": best.user_attrs}, f, indent=2)
    print(f"saved curriculum_s3p_{FAMILY}.json")


if __name__ == "__main__":
    main()
