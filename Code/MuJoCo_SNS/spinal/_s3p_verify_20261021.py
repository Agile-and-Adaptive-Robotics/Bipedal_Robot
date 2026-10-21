"""Fresh-process verification of the three stage-3p perturbation
winners (2026-10-21 resume session). Loads each s3p study's best trial
from its DB, rebuilds the EXACT s3p battery code path (set_stage(3)
+ mechanism-knob loaders + 7 runner evals + score_battery) and reports
the fresh score next to the DB value, plus per-case fall/sway. One
family per process (freshness; BLAS/attractor caveat).

Usage: python _s3p_verify_20261021.py <stock|w2lvar|syn6>
"""
import io
import importlib
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

FAMILY = sys.argv[1]

CFG = {
    "stock": dict(db="sqlite:///optuna_walk.db",
                  study="curr_s3p_perturb", mod="_curriculum", net=None),
    "w2lvar": dict(db="sqlite:///optuna_w2lvar.db",
                   study="curr_w2lvar_s3p_perturb",
                   mod="_curriculum_w2lvar", net="w2lvar"),
    "syn6": dict(db="sqlite:///optuna_syn6.db",
                 study="curr_syn6_s3p_perturb",
                 mod="_curriculum_syn6", net="syn6"),
}[FAMILY]

os.environ["AARL_NPZ"] = f"s3pverify_{FAMILY}_20261021.npz"
if CFG["net"]:
    os.environ["AARL_NET"] = CFG["net"]
else:
    os.environ.pop("AARL_NET", None)

import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
st = optuna.load_study(study_name=CFG["study"], storage=CFG["db"])
best = st.best_trial
print(f"[{FAMILY}] study {CFG['study']} best trial {best.number} "
      f"DB value {best.value:.4f}")
print(f"[{FAMILY}] params: {json.dumps(best.params)}", flush=True)

C = importlib.import_module(CFG["mod"])
import params as P

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


def load_params(p):
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


prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5

p = {**BASE_MUL, **best.params}
load_params(p)

import time
t0 = time.time()
ms = []
for name, extra in BATTERY:
    args = (["--stand-eval", "8",
             "--push-time", repr(T_EVENT),
             "--rig-scale", repr(float(p["rig_scale"]))] + list(extra))
    m = C.R.main(args)
    ms.append(m)
    print(f"[{FAMILY}] {name:9s} fell={m['bal_fell']} "
          f"push_sway={m['bal_push_sway']} nan={m['nan']}", flush=True)
wall = time.time() - t0


def score_battery(ms_):
    falls = sum(1 for m in ms_ if m["bal_fell"])
    if any(m["nan"] for m in ms_):
        return -200.0, falls
    if falls == len(ms_):
        return -150.0, falls
    q = ms_[0]

    def _ps(m):
        v = float(m["bal_push_sway"])
        return 0.30 if v < 0.0 else v

    push_sw = [_ps(m) for m, (nm, _x) in zip(ms_, BATTERY)
               if nm.startswith("push")]
    plat_sw = [_ps(m) for m, (nm, _x) in zip(ms_, BATTERY)
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


val, falls = score_battery(ms)
q = ms[0]
print(f"\n[{FAMILY}] FRESH s3p score {val:.4f} (DB {best.value:.4f}, "
      f"delta {val - best.value:+.4f}) falls {falls}/7 wall {wall:.0f}s")
print(f"[{FAMILY}] quiet: sway {float(q['bal_sway']):.4f} "
      f"tilt {float(q['bal_tilt_max']):.2f} "
      f"effort {float(q['bal_effort']):.3f} "
      f"cocon {float(q['bal_cocon']):.3f} "
      f"contact_sym {float(q['bal_contact_sym']):.3f}")

with open(f"_s3p_verify_{FAMILY}_20261021.json", "w",
          encoding="utf-8") as f:
    json.dump(dict(family=FAMILY, study=CFG["study"], db=CFG["db"],
                   trial=best.number, db_value=best.value,
                   params=best.params, fresh_score=val, falls=falls,
                   quiet={k: q[k] for k in
                          ("bal_sway", "bal_tilt_max", "bal_effort",
                           "bal_cocon", "bal_contact_sym", "bal_fell")},
                   per_case=[{"case": nm, "fell": bool(m["bal_fell"]),
                              "nan": bool(m["nan"]),
                              "push_sway": m["bal_push_sway"]}
                             for (nm, _x), m in zip(BATTERY, ms)],
                   wall_s=wall), f, indent=1, default=str)
print(f"[{FAMILY}] wrote _s3p_verify_{FAMILY}_20261021.json")
