"""Quiet-eval rig-share for the syn6 stage-3p winner (rig_scale 0.054,
wide stance, knee bias, coact) - the near-rig-free headline case.
Runs ONE quiet 8 s stand into a dedicated npz, then applies the same
rig-share math as _s3_rigshare_20261021.py."""
import io
import importlib
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

os.environ["AARL_NPZ"] = "s3p_syn6_quietshare_20261021.npz"
os.environ["AARL_NET"] = "syn6"

import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
best = optuna.load_study(
    study_name="curr_syn6_s3p_perturb",
    storage="sqlite:///optuna_syn6.db").best_trial
print(f"syn6 s3p winner trial {best.number} params "
      f"{json.dumps(best.params)}", flush=True)

C = importlib.import_module("_curriculum_syn6")
import params as P

prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
p = {**dict(prev["multipliers"]), **best.params}
p["renshaw"] = 0.5
C.set_stage(3, p)
P.G["coact"] = float(p.get("coact", 0.0))
if float(p.get("stance_width", 0.0)) > 0.0:
    os.environ["AARL_STANCE_W"] = repr(float(p["stance_width"]))
if float(p.get("knee_bias", 0.0)) > 0.0:
    os.environ["AARL_KNEE_BIAS"] = repr(float(p["knee_bias"]))
os.environ["AARL_PELVIS_TY"] = repr(float(p["pelvis_ty"]))

m = C.R.main(["--stand-eval", "8", "--push-time", "4.0",
              "--rig-scale", repr(float(p["rig_scale"]))])
print(f"quiet metrics: fell={m['bal_fell']} sway={m['bal_sway']:.4f} "
      f"tilt={m['bal_tilt_max']:.2f}")

import numpy as np
W = 72.0 * 9.81
z = np.load("s3p_syn6_quietshare_20261021.npz", allow_pickle=True)
t = z["t"]
n = min(int(np.sum(t > 0)) + 1, len(t))
i0 = int(np.searchsorted(t[:n], 4.0))
ct = z["contact"][i0:n]
q = z["q"][i0:n]
com = z["com"][i0:n]
foot = 0.5 * (ct[:, 0] + ct[:, 1])
rec = dict(
    foot_load_mean_N=float(np.mean(foot)),
    foot_load_pct_weight=float(100.0 * np.mean(foot) / W),
    rig_share_pct=float(100.0 * (W - np.mean(foot)) / W),
    contact_duty_R_L=(float((ct[:, 0] > 20).mean()),
                      float((ct[:, 1] > 20).mean())),
    pelvis_tilt_max_deg=float(np.max(np.abs(q[:, 0]))),
    knee_flex_max_deg_R_L=(float(-np.min(q[:, 4])), float(-np.min(q[:, 9]))),
    com_z_mean=float(np.mean(com[:, 2])),
    com_z_min=float(np.min(com[:, 2])),
)
print(json.dumps(rec, indent=1))
with open("_s3p_syn6_rigshare_20261021.json", "w",
          encoding="utf-8") as f:
    json.dump(rec, f, indent=1)
print("wrote _s3p_syn6_rigshare_20261021.json")
