"""REBASELINE the syn6 stage-4 winner under the FIXED kine_ref (2026-09-26).

Same recipe as rebase_syn6_s5_winner.py but for stage 4 (--no-ground +
--eval; NO kz penalty at stage 4 - _curriculum_syn6.py:238 applies the
COM-floor penalty to stage 5 only). Needed so the stage-4 adaptive
loop has a new-reference baseline: the study best (-237.253365, trial
14) was scored under the OLD left-ref and is not comparable to new
trials. Scratch npz: scratch_syn6_rebase_s4.npz (never the chain npz).
"""
import os

os.environ["AARL_NET"] = "syn6"
os.environ["AARL_NPZ"] = "scratch_syn6_rebase_s4.npz"

import json
import sys

sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

import numpy as np

import _curriculum_syn6 as CS   # also wraps stdout utf-8 at import
import runner as R

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
win = json.load(open(os.path.join(SP, "curriculum_syn6_stage4.json"),
                     encoding="utf-8"))
prev = json.loads(open(os.path.join(SP, "best_walk_params_v10.json"),
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.0
BASE_MUL["syn6"] = 1.0
BASE_MUL["syn6_brainstem"] = 0.0
p = {**BASE_MUL, **win["params"]}
CS.set_stage(4, p)
args = ["--no-ground", "--eval", "--drive", repr(p["drive"])]
print("rebase args:", args, flush=True)
print("rebase winner params:", json.dumps(win["params"]), flush=True)
m = R.main(args)
k = m.get("kine") or {}
sel = {kk: m.get(kk) for kk in ["nan", "kine_score", "kz", "tilt_max",
                                "duty"]}
sel["n_cycles_r"] = k.get("n_cycles_r")
sel["n_cycles_l"] = k.get("n_cycles_l")
sel["T_r"] = k.get("T_r")
sel["duty_kine"] = k.get("duty")
sel["knee_min"] = k.get("knee_min")
print("metrics:", sel)
# stage-4 objective EXACTLY as _curriculum_syn6.py:237-243 (no kz term)
score = max(float(m["kine_score"]), -315.0)
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
}
print("guard rails (stage 4 = air, no contact expected):", rails)
print("VERDICT:", "OK (rebased)" if all(rails.values()) else "CHECK FAILED")
