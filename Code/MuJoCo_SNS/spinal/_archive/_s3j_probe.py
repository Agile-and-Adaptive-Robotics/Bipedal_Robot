"""Probe the s3j winner (full-Deng connectome retune) both-leg state."""
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
s = optuna.load_study(study_name="curr_s3j_fullrules",
                      storage="sqlite:///optuna_walk.db")
best = max((t for t in s.trials if t.value is not None),
           key=lambda t: t.value)
p = {**mul, **st3, **best.params}
C.set_stage(3, p)
os.environ["AARL_NPZ"] = "spinal_run_s3j.npz"
R.main(["--eval", "--drive", repr(best.params["drive"])])
z = np.load("spinal_run_s3j.npz", allow_pickle=True)
k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
               ref=KR.ref_cached(), contact=z["contact"])
z.close()
out = {kk: k.get(kk) for kk in
       ("kine_score", "duty_r", "duty_l", "n_cycles_r", "n_cycles_l",
        "frozen_r", "frozen_l", "ds", "contact_frac_r",
        "contact_frac_l")}
open("_s3j_probe_out.txt", "w", encoding="utf-8").write(repr(out))
print("WROTE _s3j_probe_out.txt")
