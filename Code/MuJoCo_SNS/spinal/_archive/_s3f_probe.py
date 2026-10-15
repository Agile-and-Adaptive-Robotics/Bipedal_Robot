"""Probe the current s3f best: both-leg readout (left cycles = the
goal). Side npz; safe while the chain runs."""
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

s = optuna.load_study(study_name="curr_s3f_ground",
                      storage="sqlite:///optuna_walk.db")
best = max((t for t in s.trials if t.value is not None),
           key=lambda t: t.value)
params = best.params
print(f"probing s3f t{best.number} value {best.value:.1f}: "
      f"pm {params.get('pm_gain'):.2f} T {params.get('pm_T'):.2f}")
C.set_stage(3, {**mul, **params})
os.environ["AARL_NPZ"] = "spinal_run_probe.npz"
R.main(["--eval", "--drive", repr(params["drive"])])
z = np.load("spinal_run_probe.npz", allow_pickle=True)
k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
               ref=KR.ref_cached(), contact=z["contact"])
z.close()
if k:
    print(f"score {k['kine_score']:.1f} | duty r/l "
          f"{k.get('duty_r'):.2f}/{k.get('duty_l')} | cyc r/l "
          f"{k.get('n_cycles_r')}/{k.get('n_cycles_l')} | frozen "
          f"r/l {k.get('frozen_r')}/{k.get('frozen_l')} | knee_min_l "
          f"{k.get('knee_min_l', float('nan'))} | ds "
          f"{k.get('ds', float('nan')):.2f} | lag "
          f"{k.get('lag_rl', float('nan'))}")
else:
    print("no kine")
