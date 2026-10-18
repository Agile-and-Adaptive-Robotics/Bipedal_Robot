"""One combo of the 2026-10-16 HUMAN L/R alpha+gamma drive sweep.

Usage: D:\\Anaconda\\envs\\myo\\python.exe drive_sweep_one_20261016.py <json-spec>
json-spec = {"id":..., "gamma_l":g, "gamma_r":g, "alpha_l":a, "alpha_r":a}

Base config = campaign anchor A (s3k trial-34 winner merged over the v10
multipliers, full_rules=0, set_stage(4, ...)) - the config gated bit-exact
at -161.56754173676563 (reports_20261016/DRIVE_VERIFY_20261016.md).
Knob semantics (runner.py step loop): u[GAMMA_c_s] = drive*G[gamma_s],
u[ALPHA_c_s] = drive*G[alpha_s] with drive = the s3k winner drive
2.0773198435110425, so knob values translate to port currents of
knob*2.077 nA. One combo per process (params.G is module-global).
Appends one JSON line to reports_20261016/drive_sweep_human/results.jsonl
and prints "SWEEPRESULT <json>" on stdout.
"""
import json
import os
import sys

spec = json.loads(sys.argv[1])

import _curriculum as C
import params as P
import runner as R

base = dict(json.loads(
    open("best_walk_params_v10.json", encoding="utf-8").read())
    ["multipliers"])
base["renshaw"] = 0.5
s3k = json.loads(open(
    "reports_20260923/s3k_trial34_full_params.json",
    encoding="utf-8").read())["params"]
base.update(s3k)
base["full_rules"] = 0.0          # anchor A recipe (the campaign baseline)
C.set_stage(4, base)
for k in ("gamma_l", "gamma_r", "alpha_l", "alpha_r"):
    P.G[k] = float(spec[k])       # DRIVE_VERIFY 2026-10-16 section-6 recipe

m = R.main(["--eval", "--drive", repr(base["drive"])])
k = m.get("kine") or {}
res = dict(
    spec=spec,
    drive=base["drive"],
    kine_score=m.get("kine_score"),
    nan=bool(m.get("nan")),
    cycles_r=int(k.get("n_cycles_r", 0)),
    cycles_l=int(k.get("n_cycles_l", 0)),
    duty_r=k.get("duty_r"),
    duty_l=k.get("duty_l"),
    contact_frac_r=k.get("contact_frac_r"),
    contact_frac_l=k.get("contact_frac_l"),
    ds=k.get("ds"),
    T_r=k.get("T_r"),
    T_l=k.get("T_l"),
    frozen_r=bool(k.get("frozen_r", False)),
    frozen_l=bool(k.get("frozen_l", False)),
    duty_neural=m.get("duty"),
    burst_r=m.get("burst_r"),
    bal_fell=bool(m.get("bal_fell")),
    bal_com_z_min=m.get("bal_com_z_min"),
    tilt_max=m.get("tilt_max"),
    knee_min=m.get("knee_min"),
    dx=m.get("dx"),
)
line = json.dumps(res)
outf = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "reports_20261016", "drive_sweep_human", "results.jsonl")
os.makedirs(os.path.dirname(outf), exist_ok=True)
with open(outf, "a", encoding="utf-8") as f:
    f.write(line + "\n")
print("SWEEPRESULT " + line, flush=True)
