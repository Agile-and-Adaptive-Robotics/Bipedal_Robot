"""Exploit/validity check for the syn6 stage-5 winner (trial 0, -197.215).

Ground-walking stage: replays the EXACT winner config
  set_stage(5, {BASE_MUL + winner params}) then
  R.main(["--eval","--drive", repr(drive)])
and recomputes the score:
  max(kine_score, -315) - 20 if kz < 0.62 - 10 if tilt_max > 40.
Guard rails: m['nan'] False, m['kine'] not None, contact_frac > 0 on at
least one leg (real ground contact - the stage-4 air walk had 0.0),
RG_E_r rises >= 3 in the npz. Prints the pattern metrics for the record.
"""
import os

os.environ["AARL_NET"] = "syn6"
os.environ["AARL_NPZ"] = "syn6_s5_replay_tmp.npz"

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
print("replay args:", args, flush=True)
print("replay searched params:", json.dumps(win["params"]), flush=True)
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
print(f"recomputed score={score:.4f}  json score={win['score']:.4f}  "
      f"delta={score - win['score']:+.4f}")
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
ok = (not m["nan"]) and (m.get("kine") is not None) and \
    (max(k.get("contact_frac_r", 0.0), k.get("contact_frac_l", 0.0))
     > 0.0) and \
    abs(score - win["score"]) < 5e-3 and rises >= 3
print("VERDICT:", "GENUINE GROUND WALK" if ok else "CHECK FAILED")
