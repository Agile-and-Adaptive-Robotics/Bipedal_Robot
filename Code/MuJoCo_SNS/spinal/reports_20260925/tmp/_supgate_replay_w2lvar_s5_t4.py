"""SUPERVISOR GATE: independent replay of the w2lvar s5 PRE-campaign winner
(trial 4, params read from reports_20260925/easteregg2/pre_w2lvar_stage5.json)
under the FIXED kine_ref. Mirrors tmp/replay_w2lvar_s5_rebase.py exactly
except the winner json is the PRE snapshot (the live curriculum json now
records the newer trial-68 winner) and the scratch npz is gate-specific.
Nothing chain-owned is written.
"""
import json
import os
import sys

# NOTE: no sys.stdout wrap here - _curriculum_w2lvar.py:47 wraps stdout at
# import; an extra wrapper here orphans the previous one and closes the
# shared buffer (ValueError: I/O operation on closed file).
os.chdir(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())

os.environ["AARL_NET"] = "w2lvar"
os.environ["AARL_NPZ"] = "scratch_supgate_w2lvar_rebase.npz"

import _curriculum_w2lvar as CW
import runner as R

SP = os.getcwd()
win = json.loads(open(os.path.join(SP, "reports_20260925", "easteregg2",
                                   "pre_w2lvar_stage5.json"),
                      encoding="utf-8").read())
prev = json.loads(open(os.path.join(SP, "best_walk_params_v10.json"),
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5
p = {**BASE_MUL, **win["params"]}
CW.set_stage(5, p)
args = ["--eval", "--drive", repr(p["drive"])]
print("supgate replay args:", args, flush=True)
print("trial-4 searched params:", json.dumps(win["params"]), flush=True)
m = R.main(args)
k = m.get("kine") or {}
print("nan =", m["nan"], "kine present =", m.get("kine") is not None)
print("kine_score =", repr(float(m["kine_score"])))
print("kz =", m.get("kz"), "tilt_max =", m.get("tilt_max"),
      "duty =", m.get("duty"))
sel = {}
for kk in ["contact_frac_r", "contact_frac_l", "n_cycles_r", "n_cycles_l",
           "bilateral", "knee_min"]:
    if kk in k:
        sel["kine_" + kk] = k[kk]
print("metrics:", sel)
TARGET = -165.33596254795964
got = float(m["kine_score"])
print(f"TARGET (tuner ee2_rebase_w2lvar_s5.log) = {TARGET!r}")
print(f"DELTA = {got - TARGET:+.3e}")
print("SUPGATE REPLAY",
      "MATCH" if abs(got - TARGET) < 1e-6 else "MISMATCH", flush=True)
