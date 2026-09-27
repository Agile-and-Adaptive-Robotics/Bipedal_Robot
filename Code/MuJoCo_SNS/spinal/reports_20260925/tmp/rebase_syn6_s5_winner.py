"""REBASELINE the syn6 stage-5 winner under the FIXED kine_ref (2026-09-26).

Copy of reports_20260925/tmp/replay_syn6_s5_winner.py (the 09-25
bit-exact exploit check) adapted for the kine_ref left-cycle fix:
the objective changed BY DESIGN, so the old template's
abs(score - win["score"]) < 5e-3 assertion no longer applies - this
script reports old -> new instead and keeps the guard rails as info
(m['nan'] False, m['kine'] not None, contact_frac > 0 on one leg,
RG_E_r rises >= 3 in the npz).
Scratch npz per ask: scratch_syn6_rebase.npz (never the chain npz).
"""
import os

os.environ["AARL_NET"] = "syn6"
os.environ["AARL_NPZ"] = "scratch_syn6_rebase.npz"

import json
import sys

sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

import numpy as np

import _curriculum_syn6 as CS   # also wraps stdout utf-8 at import
import runner as R

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
win = json.load(open(os.path.join(SP, "curriculum_syn6_stage5.json"),
                     encoding="utf-8"))
prev = json.loads(open(os.path.join(SP, "best_walk_params_v10.json"),
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.0
BASE_MUL["syn6"] = 1.0
BASE_MUL["syn6_brainstem"] = 0.0
p = {**BASE_MUL, **win["params"]}
CS.set_stage(5, p)
args = ["--eval", "--drive", repr(p["drive"])]
print("rebase args:", args, flush=True)
print("rebase winner params:", json.dumps(win["params"]), flush=True)
m = R.main(args)
k = m.get("kine") or {}
sel = {kk: m.get(kk) for kk in ["nan", "kine_score", "kz", "tilt_max",
                                "duty"]}
sel["contact_frac_r"] = k.get("contact_frac_r")
sel["contact_frac_l"] = k.get("contact_frac_l")
sel["n_cycles_r"] = k.get("n_cycles_r")
sel["n_cycles_l"] = k.get("n_cycles_l")
sel["bilateral"] = k.get("bilateral")
sel["T_r"] = k.get("T_r")
sel["lag_rl"] = k.get("lag_rl")
sel["duty_kine"] = k.get("duty")
sel["knee_min"] = k.get("knee_min")
sel["ds"] = k.get("ds")
print("metrics:", sel)
score = max(float(m["kine_score"]), -315.0)
if float(m["kz"]) < 0.62:
    score -= 20.0
if float(m["tilt_max"]) > 40.0:
    score -= 10.0
print(f"REBASE score={score:.4f}  OLD json score={win['score']:.4f}  "
      f"delta={score - win['score']:+.4f}  (kine_ref left-cycle fix "
      f"applied 2026-09-26)")
z = np.load(os.environ["AARL_NPZ"], allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(x) for x in z["neuro_names"]]
joints = [str(x) for x in z["key_joints"]]
i_rge, i_knee = names.index("RG_E_r"), joints.index("knee_angle_r")
msk = t >= 5.0
rge = neuro[msk, i_rge]
knee = q[msk, i_knee]
on = rge > 0.5 * max(float(rge.max()), 1e-9)
rises = int(np.sum(np.diff(on.astype(int)) == 1))
print(f"RG_E_r rises(after t=5)={rises} span="
      f"{float(rge.max()) - float(rge.min()):.3f}; knee range "
      f"{float(knee.min()):.1f}..{float(knee.max()):.1f} deg; cfg="
      f"{z['cfg']}")
rails = {
    "nan_false": not m["nan"],
    "kine_not_none": m.get("kine") is not None,
    "contact_some": max(k.get("contact_frac_r", 0.0),
                        k.get("contact_frac_l", 0.0)) > 0.0,
    "rises_ge_3": rises >= 3,
}
print("guard rails:", rails)
print("VERDICT:", "GENUINE GROUND WALK (rebased)" if all(rails.values())
      else "CHECK FAILED")
