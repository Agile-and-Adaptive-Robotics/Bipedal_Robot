"""Ib ABLATION on the standing winners (2026-10-20, Ben: 'Ib load
feedback prevents collapse' must be numbers, not an assertion). Re-runs
each winner's quiet standing eval with EVERY extensor-Ib pathway gain
the family's builder reads ZEROED, then prints the metrics next to the
intact numbers recorded in _s3_verify_*_20261020.json.

Usage: python _s3_ib_ablate_20261020.py <w2l_bal|w2l_fix|syn6_bal|stock>
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

CASE = sys.argv[1]
CASES = {
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

os.environ["AARL_NPZ"] = f"s3ibablate_{CASE}_20261020.npz"
if variant:
    os.environ["AARL_NET"] = variant
else:
    os.environ.pop("AARL_NET", None)

import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
st = optuna.load_study(study_name=study_name, storage=db)
best = st.best_trial
print(f"[{CASE}] Ib-ABLATION of {study_name} best trial "
      f"{best.number} (DB value {best.value:.4f})")
print(f"[{CASE}] params: {json.dumps(best.params)}")

import importlib
C = importlib.import_module(mod_name)
import params as P

prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5
p = {**BASE_MUL, **best.params}
C.set_stage(3, p)

# ---- zero every Ib-pathway gain the family builder reads
ABLATE_KEYS = ("ib_rge", "ib_e_central", "ib_to_mn_inh", "ib_exc_to_mn",
               "ib_group_exc")
for k in ABLATE_KEYS:
    if k in P.G:
        print(f"[{CASE}] {k}: {P.G[k]} -> 0.0")
        P.G[k] = 0.0
try:   # w2lvar's own motif table
    import build_network_w2lvar as BW
    if BW.GAINS.get("w2l_ib_central", 0.0) > 0.0:
        print(f"[{CASE}] GAINS.w2l_ib_central: "
              f"{BW.GAINS['w2l_ib_central']} -> 0.0")
        BW.GAINS["w2l_ib_central"] = 0.0
except ImportError:
    pass

m = C.R.main(["--stand-eval", "8",
              "--rig-scale", repr(float(p["rig_scale"]))])
print(f"\n[{CASE}] Ib-ABLATED STANDING METRICS (8 s, rig_scale="
      f"{p['rig_scale']:.4f}):")
for k in ("nan", "bal_fell", "bal_sway", "bal_sway_rms",
          "bal_tilt_max", "bal_contact_sym", "bal_com_z_min"):
    print(f"    {k:16s} = {m.get(k)}")
try:
    ref = json.load(open(f"_s3_verify_{CASE}_20261020.json",
                         encoding="utf-8"))["metrics"]
    print(f"[{CASE}] INTACT (for comparison): fell={ref['bal_fell']} "
          f"sway={ref['bal_sway']:.4f} tilt={ref['bal_tilt_max']:.2f} "
          f"com_z_min={ref['bal_com_z_min']:.3f} "
          f"sym={ref['bal_contact_sym']:.3f}")
except FileNotFoundError:
    print(f"[{CASE}] (no intact reference file yet)")
with open(f"_s3_ibablate_{CASE}_20261020.json", "w",
          encoding="utf-8") as f:
    json.dump(dict(case=CASE, study=study_name, trial=best.number,
                   params=best.params, ablated=ABLATE_KEYS,
                   metrics={k: (bool(v) if isinstance(v, bool) else v)
                            for k, v in m.items()
                            if not isinstance(v, dict)}), f, indent=1,
              default=str)
print(f"[{CASE}] wrote _s3_ibablate_{CASE}_20261020.json")
