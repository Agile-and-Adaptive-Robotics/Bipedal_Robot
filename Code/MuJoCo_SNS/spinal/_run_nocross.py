"""Run the SEVERED-L/R condition (all crossed communication zeroed)
of the s3i winner - the config where the left leg stepped for the
first time - and save its run npz for figures + animation.
No stdout wrapper: _curriculum re-wraps on import (closes ours)."""
import json
import os
import sys
from pathlib import Path

import _curriculum as C
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(3, {**mul, **st3,
                "contra_swing": 0.0, "contra_kinh": 0.0,
                "ia_f_contra_f": 0.0, "c1_gain": 0.1, "v3_gain": 0.0})
os.environ["AARL_NPZ"] = "spinal_run_nocross.npz"
R.main(["--eval", "--drive", repr(st3["drive"])])
print("saved spinal_run_nocross.npz (severed-L/R, s3i t20 config)")
