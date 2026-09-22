"""Probe s3d: leaders' pf_gain/contra_swing + a rerun of the current
best with the both-leg readout (left cycles = the architectural goal).
Side npz to avoid touching the running chain's spinal_run.npz."""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
import numpy as np

import _curriculum as C
import kine_ref as KR
import runner as R
import os

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]

s = optuna.load_study(study_name="curr_s3d_ground",
                      storage="sqlite:///optuna_walk.db")
done = [(t.value, t.number, t.params) for t in s.trials
        if t.value is not None]
done.sort(key=lambda x: -x[0])
print("leaders (value, trial, pf_gain, contra_swing, drive, pelvis):")
for v, n, p in done[:8]:
    print(f"  {v:8.1f} t{n:3d}  pf {p.get('pf_gain', -1):.2f}  "
          f"cs {p.get('contra_swing', -1):.2f}  d {p['drive']:.2f}  "
          f"ty {p.get('pelvis_ty', -1):.3f}")

v, n, params = done[0]
print(f"\nrerunning best t{n} for the left-leg verdict...")
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(3, {**mul, **params})
os.environ["AARL_NPZ"] = "spinal_run_probe.npz"
R.main(["--eval", "--drive", repr(params["drive"])])
z = np.load("spinal_run_probe.npz", allow_pickle=True)
k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
               ref=KR.ref_cached(), contact=z["contact"])
z.close()
if k:
    print(f"score {k['kine_score']:.1f} | duty r/l "
          f"{k.get('duty_r'):.2f}/{k.get('duty_l')} | "
          f"cyc r/l {k.get('n_cycles_r')}/{k.get('n_cycles_l')} | "
          f"frozen_l {k.get('frozen_l')} | knee_min_l "
          f"{k.get('knee_min_l', float('nan'))}")
else:
    print("no kine on rerun")
