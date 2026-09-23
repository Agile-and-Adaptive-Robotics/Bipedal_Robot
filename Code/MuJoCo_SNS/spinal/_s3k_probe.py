"""Probe the s3k winner (severed-L/R bilateral retune)."""
import io
import json
import os
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import _curriculum as C
import kine_ref as KR
import optuna
import runner as R

optuna.logging.set_verbosity(optuna.logging.WARNING)
mul = json.loads((Path(__file__).parent / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((Path(__file__).parent / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
s = optuna.load_study(study_name="curr_s3k_nocross",
                      storage="sqlite:///optuna_walk.db")
best = max((t for t in s.trials if t.value is not None),
           key=lambda t: t.value)
p = {**mul, **st3, **best.params, "no_cross": 1.0}
C.set_stage(3, p)
os.environ["AARL_NPZ"] = "spinal_run_s3k.npz"
R.main(["--eval", "--drive", repr(best.params["drive"])])
z = np.load("spinal_run_s3k.npz", allow_pickle=True)
k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
               ref=KR.ref_cached(), contact=z["contact"])
z.close()
out = {kk: (round(k[kk], 3) if isinstance(k.get(kk), float)
            else k.get(kk))
       for kk in ("kine_score", "duty_r", "duty_l", "n_cycles_r",
                  "n_cycles_l", "frozen_r", "frozen_l", "ds",
                  "contact_frac_r", "contact_frac_l", "knee_min_r",
                  "knee_min_l", "mean_hip_r", "mean_hip_l", "T_r",
                  "T_l", "period_cv_r", "period_cv_l", "lag_rl")}
open("_s3k_probe_out.txt", "w", encoding="utf-8").write(repr(out))
print("WROTE _s3k_probe_out.txt")
