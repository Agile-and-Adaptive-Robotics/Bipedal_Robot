"""STEP 2 rebaseline: w2lvar stage-5 winner under the FIXED kine_ref.

The 09-25 stage-5 winner (curriculum_w2lvar_stage5.json, score -165.041)
was scored under the buggy left reference (left cycle cut at t=0.005 s,
np.interp edge-filling froze ref['l'] for the first 39.6% of the cycle;
fixed 2026-09-26 in kine_ref.load_reference via _first_pair_in_ik).
This script re-evaluates the EXACT winner config under the fixed
reference with the scratch npz (never the chain npz):
  set_stage(5, {BASE_MUL + winner params}) then
  R.main(["--eval", "--drive", repr(drive)])
and scores with the unchanged stage-5 objective:
  nan -> -400 | no cycles -> -320 |
  max(kine_score, -315) - 20 if kz < 0.62 - 10 if tilt_max > 40.
Loading mirrors harvest_w2lvar_s5.py (BASE = best_walk_params_v10
multipliers, renshaw 0.5) and the objective path is _curriculum_w2lvar
.set_stage itself. Report old -> new; NOTHING is written to the chain
jsons / npz.
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
os.chdir(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())

os.environ["AARL_NET"] = "w2lvar"
os.environ["AARL_NPZ"] = "scratch_w2lvar_rebase.npz"

import numpy as np

import _curriculum_w2lvar as CW
import runner as R

SP = os.getcwd()
win = json.loads(open(os.path.join(SP, "curriculum_w2lvar_stage5.json"),
                      encoding="utf-8").read())
prev = json.loads(open(os.path.join(SP, "best_walk_params_v10.json"),
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5
p = {**BASE_MUL, **win["params"]}
CW.set_stage(5, p)
args = ["--eval", "--drive", repr(p["drive"])]
print("rebase args:", args, flush=True)
print("winner searched params:", json.dumps(win["params"]), flush=True)
m = R.main(args)
k = m.get("kine") or {}
sel = {kk: m.get(kk) for kk in ["nan", "kine_score", "kz", "tilt_max",
                                "duty"]}
for kk in ["contact_frac_r", "contact_frac_l", "n_cycles_r", "n_cycles_l",
           "bilateral", "T_r", "lag_rl", "duty", "knee_min", "ds"]:
    if kk in k:
        sel["kine_" + kk] = k[kk]
print("metrics:", sel)
if m["nan"]:
    score, extra = -400.0, "NAN"
elif m.get("kine") is None:
    score, extra = -320.0, "NO CYCLES (frozen sentinel)"
else:
    score = max(float(m["kine_score"]), -315.0)
    kz_pen = float(m["kz"]) < 0.62
    tilt_pen = float(m["tilt_max"]) > 40.0
    if kz_pen:
        score -= 20.0
    if tilt_pen:
        score -= 10.0
    extra = f"kz_pen={kz_pen} tilt_pen={tilt_pen}"
print(f"REBASE: {extra}")
print(f"old score (buggy left ref) = {win['score']:.6f}")
print(f"new score (fixed left ref) = {score:.6f}")
print(f"delta = {score - win['score']:+.6f}")

# rhythm / pattern readout from the scratch npz (exploit-check context)
try:
    z = np.load(os.environ["AARL_NPZ"], allow_pickle=True)
    t, q, neuro = z["t"], z["q"], z["neuro"]
    names = [str(x) for x in z["neuro_names"]]
    joints = [str(x) for x in z["key_joints"]]
    i_rge = names.index("RG_E_r")
    i_knee_r = joints.index("knee_angle_r")
    i_knee_l = joints.index("knee_angle_l")
    msk = t >= 5.0
    rge = neuro[msk, i_rge]
    on = rge > 0.5 * max(float(rge.max()), 1e-9)
    rises = int(np.sum(np.diff(on.astype(int)) == 1))
    kr = q[msk, i_knee_r]
    kl = q[msk, i_knee_l]
    print(f"RG_E_r rises(t>=5)={rises} span="
          f"{float(rge.max()) - float(rge.min()):.3f}; knee_r "
          f"{float(kr.min()):.1f}..{float(kr.max()):.1f} deg; knee_l "
          f"{float(kl.min()):.1f}..{float(kl.max()):.1f} deg; cfg="
          f"{z['cfg']}")
except Exception as e:  # noqa: BLE001 - readout only, never fail rebase
    print(f"npz readout skipped: {e!r}")
