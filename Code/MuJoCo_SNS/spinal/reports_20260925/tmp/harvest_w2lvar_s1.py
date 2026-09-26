"""Harvest check for w2lvar stage 1: trial distribution + winner exploit
re-check (re-run winner params through the exact stage-1 objective path)."""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
sys.path.insert(0, os.getcwd())  # run-by-path puts tmp\ on sys.path, not cwd

import numpy as np
import optuna

st = optuna.load_study(study_name="curr_w2lvar_s1_air_deaff",
                       storage="sqlite:///optuna_w2lvar.db")
vals = [t.value for t in st.trials]
print(f"n_trials={len(st.trials)} "
      f"states={[t.state.name for t in st.trials].count('COMPLETE')} COMPLETE")
print("values:", [None if v is None else round(v, 3) for v in vals])
sent = sum(1 for v in vals if v is not None and v <= -199.5)
gate = sum(1 for v in vals if v is not None and -10.5 < v < -3.0)
real = sum(1 for v in vals if v is not None and v > 9.0)
print(f"sentinel(<=-199.5)={sent}  rhythm-gate(-10.5..-3)={gate}  "
      f"real-rhythm(>9)={real}")

win = json.loads(open("curriculum_w2lvar_stage1.json",
                      encoding="utf-8").read())
p = win["params"]
print("winner params:", p)

# EXACT stage-1 objective reproduction for the winner (code copied from
# _curriculum_w2lvar.py objective(), stage-1 branch, with the same
# set_stage prelude).
import optuna_walk_v10 as OW
import params as P
import runner as R

os.environ["AARL_NET"] = "w2lvar"
os.environ["AARL_NPZ"] = "spinal_run_w2lvar.npz"
BASE = dict(json.loads(open("best_walk_params_v10.json",
                            encoding="utf-8").read())["multipliers"])
BASE["renshaw"] = 0.5
full = {**BASE, **p}
OW.load_fitted_baseline()
OW.set_params(full)
P.TAU["rg_nap_h"] = float(p["rg_nap_h"])
os.environ.pop("AARL_KY", None)
os.environ.pop("AARL_PELVIS_TY", None)

args = ["--no-ground", "--no-afferents", "--no-interleg", "--time", "14",
        "--drive", repr(p["drive"])]
R.main(args)
z = np.load("spinal_run_w2lvar.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
names = [str(x) for x in z["neuro_names"]]
joints = [str(x) for x in z["key_joints"]]
i_rge = names.index("RG_E_r")
i_knee = joints.index("knee_angle_r")
m = (t >= 5.0) & (t <= 17.0)
knee = q[m, i_knee]
rge = neuro[m, i_rge]
on = rge > 0.5 * max(rge.max(), 1e-9)
rises = int(np.sum(np.diff(on.astype(int)) == 1))
kmin = float(knee.min())
score = 3.0 * rises + 0.5 * (-kmin)
print(f"REPRO: rises={rises} knee_min={kmin:.4f} rge_range="
      f"{float(rge.max()) - float(rge.min()):.4f} "
      f"recomputed_score={score:.6f} json_score={win['score']:.6f} "
      f"delta={abs(score - win['score']):.6f}")
print("exploit gate: rises>=3 ->", "PASS (real rhythm)" if rises >= 3
      else "FAIL (static-pose exploit)")
