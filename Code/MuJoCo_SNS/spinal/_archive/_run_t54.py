"""Rerun s3c trial 54 (corrected-objective best) and render the ISB
figures from its npz."""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
import _curriculum as C
import kine_ref as KR
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
s = optuna.load_study(study_name="curr_s3c_ground",
                      storage="sqlite:///optuna_walk.db")
params = s.trials[54].params
C.set_stage(3, {**mul, **params})
R.main(["--eval", "--drive", repr(params["drive"])])
z = np_holder = None
import numpy as np

z = np.load("spinal_run.npz", allow_pickle=True)
k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
               ref=KR.ref_cached(), contact=z["contact"])
z.close()
print("t54 corrected score:", k["kine_score"] if k else None)
if k:
    for key in ("duty_r", "duty_l", "n_cycles_r", "n_cycles_l",
                "frozen_l", "ds", "knee_min_r", "mean_hip_r",
                "mean_ankle_r", "mean_ankle_l", "T_r"):
        print(f"  {key}: {k.get(key)}")
