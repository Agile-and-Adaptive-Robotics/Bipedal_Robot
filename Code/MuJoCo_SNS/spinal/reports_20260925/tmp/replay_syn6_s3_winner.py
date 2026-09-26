"""Exploit/validity check for the syn6 stage-3 winner (trial 11, 37.463).

Standing-balance stage: replays the EXACT winner config
  set_stage(3, {BASE_MUL + winner params}) then
  R.main(["--stand-eval","8","--rig-scale", repr(rig_scale)])
and recomputes the score from the returned metrics:
  100 - 400*bal_sway - 1.0*bal_tilt_max - 40*|0.5 - bal_contact_sym|
Guard rails: m['nan'] must be False, m['bal_fell'] must be False, and the
recomputed score must match the json. Reports rig_scale (the weaning level
the balance was proven at) explicitly.
"""
import os

os.environ["AARL_NET"] = "syn6"
os.environ["AARL_NPZ"] = "syn6_s3_replay_tmp.npz"

import json
import sys

sys.path.insert(0, r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

import _curriculum_syn6 as CS   # also wraps stdout utf-8 at import
import runner as R

SP = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
win = json.load(open(os.path.join(SP, "curriculum_syn6_stage3.json"),
                     encoding="utf-8"))
prev = json.loads(open(os.path.join(SP, "best_walk_params_v10.json"),
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.0
BASE_MUL["syn6"] = 1.0
BASE_MUL["syn6_brainstem"] = 0.0
p = {**BASE_MUL, **win["params"]}
CS.set_stage(3, p)
args = ["--stand-eval", "8", "--rig-scale",
        repr(float(p["rig_scale"]))]
print("replay args:", args, flush=True)
print("replay searched params:", json.dumps(win["params"]), flush=True)
m = R.main(args)
keys = ["nan", "bal_fell", "bal_sway", "bal_sway_rms", "bal_tilt_max",
        "bal_contact_sym", "kz", "tilt_max"]
print("metrics:", {k: m.get(k) for k in keys})
score = (100.0
         - 400.0 * float(m["bal_sway"])
         - 1.0 * float(m["bal_tilt_max"])
         - 40.0 * abs(0.5 - float(m["bal_contact_sym"])))
print(f"recomputed score={score:.4f}  json score={win['score']:.4f}  "
      f"delta={score - win['score']:+.4f}")
ok = (not m["nan"]) and (not m["bal_fell"]) and \
    abs(score - win["score"]) < 5e-3
print("VERDICT:", "GENUINE STAND" if ok else "CHECK FAILED")
