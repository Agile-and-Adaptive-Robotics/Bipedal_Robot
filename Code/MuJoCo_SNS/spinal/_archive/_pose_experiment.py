"""POSE-LATCH TEST (2026-09-21): run the overall corrected-best config
(s3c t54) with the SYMMETRIC start pose. If the planted/frozen leg
follows the POSE (or both legs cycle from a symmetric start), the
frozen-left failure is a pose latch, not an architecture limit.
Also runs the s3d-best (t33) both ways for comparison."""
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

db = optuna.load_study(study_name="curr_s3d_ground",
                       storage="sqlite:///optuna_walk.db")
t33 = dict(db.trials[33].params)
s3c = optuna.load_study(study_name="curr_s3c_ground",
                        storage="sqlite:///optuna_walk.db")
t54 = dict(s3c.trials[54].params)

os.environ["AARL_NPZ"] = "spinal_run_pose.npz"

for label, params, pose in (("t54-s3c-best SYMMETRIC", t54,
                             "symmetric"),
                            ("t33-s3d-best SYMMETRIC", t33,
                             "symmetric")):
    os.environ["AARL_POSE"] = pose
    C.set_stage(3, {**mul, **params})
    m = R.main(["--eval", "--drive", repr(params["drive"])])
    z = np.load("spinal_run_pose.npz", allow_pickle=True)
    k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                   ref=KR.ref_cached(), contact=z["contact"])
    z.close()
    if k:
        print(f"[{label}] score {k['kine_score']:.1f} | duty r/l "
              f"{k.get('duty_r'):.2f}/{k.get('duty_l')} | cyc r/l "
              f"{k.get('n_cycles_r')}/{k.get('n_cycles_l')} | frozen "
              f"r/l {k.get('frozen_r')}/{k.get('frozen_l')} | "
              f"contactfrac r/l {k.get('contact_frac_r'):.2f}/"
              f"{k.get('contact_frac_l'):.2f} | tilt "
              f"{m['tilt_max']:.1f} kz {m['kz']:.2f}", flush=True)
    else:
        print(f"[{label}] NO KINE (nan {m['nan']}, tilt "
              f"{m['tilt_max']:.1f}, kz {m['kz']:.2f})", flush=True)
