"""Fresh-process standing-eval verification of the existing stage-3
balance winners (2026-10-20 audit). Loads each study's best trial
params straight from its DB, rebuilds the exact stage-3 objective code
path (set_stage(3, merged_params) -> runner --stand-eval 8 --rig-scale),
and prints the raw bal_* metrics so we can judge whether the winners
actually stand or just score. One case per process (freshness).

Usage: python _s3_verify_20261020.py <w2l_bal|w2l_fix|syn6_bal|stock>
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

CASE = sys.argv[1]

CASES = {
    # study db, study name, variant module or None (stock)
    "w2l_bal": ("sqlite:///optuna_w2lvar.db", "curr_w2lvar_s3_balance",
                "_curriculum_w2lvar", "w2lvar"),
    "w2l_fix": ("sqlite:///optuna_w2lvar.db",
                "curr_w2lvar_s3_fix_20261003", "_curriculum_w2lvar",
                "w2lvar"),
    "syn6_bal": ("sqlite:///optuna_syn6.db", "curr_syn6_s3_balance",
                 "_curriculum_syn6", "syn6"),
    "stock": ("sqlite:///optuna_walk.db", "curr_s3_balance",
              "_curriculum", None),
}
db, study_name, mod_name, variant = CASES[CASE]

# env BEFORE any spinal import (variant net selection + npz redirect)
os.environ["AARL_NPZ"] = f"s3verify_{CASE}_20261020.npz"
if variant:
    os.environ["AARL_NET"] = variant
else:
    os.environ.pop("AARL_NET", None)

import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
st = optuna.load_study(study_name=study_name, storage=db)
best = st.best_trial
print(f"[{CASE}] study {study_name} best trial {best.number} "
      f"value {best.value:.4f}")
print(f"[{CASE}] params: {json.dumps(best.params)}")

import importlib
C = importlib.import_module(mod_name)

prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5

p = {**BASE_MUL, **best.params}
C.set_stage(3, p)

import time
t0 = time.time()
m = C.R.main(["--stand-eval", "8",
              "--rig-scale", repr(float(p["rig_scale"]))])
wall = time.time() - t0

print(f"\n[{CASE}] RAW STANDING METRICS (8 s stand, DRIVE 0, "
      f"rig_scale={p['rig_scale']:.4f}, wall {wall:.1f}s):")
for k in ("nan", "bal_fell", "bal_sway", "bal_sway_rms",
          "bal_tilt_max", "bal_contact_sym", "bal_com_z_min",
          "bal_push_sway", "t_end", "kz", "tilt_max"):
    print(f"    {k:16s} = {m.get(k)}")

score = None
if not m["nan"] and not m["bal_fell"]:
    score = (100.0 - 400.0 * float(m["bal_sway"])
             - 1.0 * float(m["bal_tilt_max"])
             - 40.0 * abs(0.5 - float(m["bal_contact_sym"])))
print(f"[{CASE}] recomputed stage-3 score = "
      f"{score if score is not None else 'SENTINEL (nan/fell)'} "
      f"[DB value {best.value:.4f}]")

with open(f"_s3_verify_{CASE}_20261020.json", "w",
          encoding="utf-8") as f:
    json.dump(dict(case=CASE, study=study_name, db=db,
                   trial=best.number, db_value=best.value,
                   params=best.params,
                   metrics={k: (bool(v) if isinstance(v, bool) else v)
                            for k, v in m.items()
                            if not isinstance(v, dict)},
                   recomputed_score=score, wall_s=wall), f, indent=1,
              default=str)
print(f"[{CASE}] wrote _s3_verify_{CASE}_20261020.json")
