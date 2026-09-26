"""Exploit check for the syn6 stage-1 winner (trial 17, score 95.639).

Replays the EXACT winner configuration the study evaluated:
  set_stage(1, {BASE_MUL + winner params}) then
  R.main(["--no-ground","--no-afferents","--no-interleg","--time","14",
          "--drive", repr(drive)])
writes to a scratch npz (AARL_NPZ; bare name - runner.py:1630 joins it
with HERE), and recomputes the air objective's ingredients from the npz:
RG_E_r burst rises (0.5*max threshold), knee range in the t=5..17 window,
score = 3*rises + 0.5*(-knee_min). Exploit test: rises >= 3 and a moving
trace (not a static pose), score NOT on a sentinel floor.
"""
import os

os.environ["AARL_NET"] = "syn6"
os.environ["AARL_NPZ"] = "syn6_s1_replay_tmp.npz"

import json
import sys

sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

import numpy as np

import _curriculum_syn6 as CS   # also wraps stdout utf-8 at import
import runner as R

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
win = json.load(open(os.path.join(SP, "curriculum_syn6_stage1.json"),
                     encoding="utf-8"))
prev = json.loads(open(os.path.join(SP, "best_walk_params_v10.json"),
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.0
BASE_MUL["syn6"] = 1.0
BASE_MUL["syn6_brainstem"] = 0.0
p = {**BASE_MUL, **win["params"]}
CS.set_stage(1, p)
args = ["--no-ground", "--no-afferents", "--no-interleg", "--time", "14",
        "--drive", repr(p["drive"])]
print("replay args:", args, flush=True)
print("replay searched params:", json.dumps(win["params"]), flush=True)
R.main(args)
z = np.load(os.environ["AARL_NPZ"], allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(x) for x in z["neuro_names"]]
joints = [str(x) for x in z["key_joints"]]
i_rge, i_knee = names.index("RG_E_r"), joints.index("knee_angle_r")
print("replay npz cfg:", z["cfg"], "| RG_E_r idx", i_rge,
      "knee_angle_r idx", i_knee)
msk = (t >= 5.0) & (t <= 17.0)
rge = neuro[msk, i_rge]
knee = q[msk, i_knee]
on = rge > 0.5 * max(float(rge.max()), 1e-9)
rises = int(np.sum(np.diff(on.astype(int)) == 1))
score = 3.0 * rises + 0.5 * (-float(knee.min()))
print(f"RG_E_r span: {float(rge.min()):.4f}..{float(rge.max()):.4f} "
      f"(delta {float(rge.max()) - float(rge.min()):.4f})")
print(f"knee window: min {float(knee.min()):.2f} max {float(knee.max()):.2f}")
print(f"RISES={rises}")
print(f"recomputed score={score:.4f}  json score={win['score']:.4f}  "
      f"delta={score - win['score']:+.4f}")
print("VERDICT:", "REAL RHYTHM" if (rises >= 3 and
      (float(rge.max()) - float(rge.min())) >= 1.0) else "EXPLOIT/GATE")
