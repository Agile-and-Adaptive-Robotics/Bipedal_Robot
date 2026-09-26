"""Exploit/validity check for the syn6 stage-4 winner (trial 14, -237.253).

Contact-free walking stage: replays the EXACT winner config
  set_stage(4, {BASE_MUL + winner params}) then
  R.main(["--no-ground","--eval","--drive", repr(drive)])
and recomputes the score: max(kine_score, -315) minus the 10-point tilt
penalty if tilt_max > 40. Guard rails: m['nan'] False, m['kine'] not None
(cycles detected -> not a frozen/static pose). Extra oscillation evidence
from the replay npz: RG_E_r rises and knee range.
"""
import os

os.environ["AARL_NET"] = "syn6"
os.environ["AARL_NPZ"] = "syn6_s4_replay_tmp.npz"

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
print("replay args:", args, flush=True)
print("replay searched params:", json.dumps(win["params"]), flush=True)
m = R.main(args)
keys = ["nan", "kine", "kine_score", "kz", "tilt_max", "duty", "cycles",
        "cyc_r", "cyc_l", "knee_min", "hip_amp", "ankle_amp"]
print("metrics:", {k: m.get(k) for k in keys})
score = max(float(m["kine_score"]), -315.0)
if float(m["tilt_max"]) > 40.0:
    score -= 10.0
print(f"recomputed score={score:.4f}  json score={win['score']:.4f}  "
      f"delta={score - win['score']:+.4f}")
# extra oscillation evidence from the replay npz
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
      f"{float(knee.min()):.1f}..{float(knee.max()):.1f} deg")
ok = (not m["nan"]) and (m.get("kine") is not None) and \
    abs(score - win["score"]) < 5e-3 and rises >= 3
print("VERDICT:", "GENUINE AIR-GAIT PATTERN" if ok else "CHECK FAILED")
