"""Probe the current s3g best: both-leg readout."""
import json
import os
import sys
from pathlib import Path

import numpy as np

import _curriculum as C
import kine_ref as KR
import optuna
import runner as R

optuna.logging.set_verbosity(optuna.logging.WARNING)

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5

s = optuna.load_study(study_name="curr_s3i_ground",
                      storage="sqlite:///optuna_walk.db")
best = max((t for t in s.trials if t.value is not None),
           key=lambda t: t.value)
C.set_stage(3, {**mul, **best.params})
os.environ["AARL_NPZ"] = "spinal_run_probe.npz"
R.main(["--eval", "--drive", repr(best.params["drive"])])
z = np.load("spinal_run_probe.npz", allow_pickle=True)
k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
               ref=KR.ref_cached(), contact=z["contact"])
z.close()
out = {kk: k.get(kk) for kk in
       ("kine_score", "duty_r", "duty_l", "n_cycles_r", "n_cycles_l",
        "frozen_l", "knee_min_l", "ds", "contact_frac_l")}
with open("_s3g_probe_out.txt", "w", encoding="utf-8") as f:
    f.write(repr(out) + "\n")
print("WROTE _s3g_probe_out.txt")
