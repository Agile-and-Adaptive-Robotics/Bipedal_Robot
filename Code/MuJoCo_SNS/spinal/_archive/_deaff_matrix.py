"""UNILATERAL DEAFFERENTATION MATRIX (Ben 2026-09-21: "if one leg is
frozen and the other stepping and you don't deafferent the
contralateral communication, I will have questions on whether you
have bad architecture on one side").

Five conditions on the s3i winner (trial 20, -148.7):
  base      : as tuned
  no-cross  : all CROSSED communication zeroed (contra_swing,
              contra_kinh, ia_f_contra_f, c1, v3) - commissural isolation
  deaff-L   : LEFT leg fully deafferented (all its Ia/II/Ib + heel/toe/
              load/AFF inputs = 0) - tests "the left is jammed by its
              own afferents"
  deaff-R   : RIGHT leg deafferented - tests "the right steps only
              because the left's afferents drive it"
  both      : crossed OFF + left deafferented (isolated right side)

Readout per condition: score, per-leg duty/cycles/frozen, contact.
Side npz. Each eval ~2.5 min.
"""
import io
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]

import _curriculum as C
import kine_ref as KR
import runner as R

CASES = [
    ("base", {}, ""),
    ("no-cross", {"contra_swing": 0.0, "contra_kinh": 0.0,
                  "ia_f_contra_f": 0.0, "c1_gain": 0.1, "v3_gain": 0.0},
     ""),
    ("deaff-L", {}, "l"),
    ("deaff-R", {}, "r"),
    ("no-cross+deaff-L", {"contra_swing": 0.0, "contra_kinh": 0.0,
                          "ia_f_contra_f": 0.0, "c1_gain": 0.1,
                          "v3_gain": 0.0}, "l"),
]

C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
rows = []
for name, over, deaff in CASES:
    C.set_stage(3, {**mul, **st3, **over})
    os.environ["AARL_DEAFF"] = deaff
    os.environ["AARL_NPZ"] = "spinal_run_deaff.npz"
    m = R.main(["--eval", "--drive", repr(st3["drive"])])
    z = np.load("spinal_run_deaff.npz", allow_pickle=True)
    k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                   ref=KR.ref_cached(), contact=z["contact"])
    z.close()
    row = {"case": name, "score": None if k is None
           else round(k["kine_score"], 1),
           "duty_r": None if k is None else
           round(k.get("duty_r") or -1, 2),
           "duty_l": None if k is None else
           round(k["duty_l"], 2) if k.get("duty_l") is not None else None,
           "cyc_r": None if k is None else k.get("n_cycles_r"),
           "cyc_l": None if k is None else k.get("n_cycles_l"),
           "frozen_l": None if k is None else k.get("frozen_l"),
           "frozen_r": None if k is None else k.get("frozen_r"),
           "nan": m["nan"], "kz": round(m["kz"], 2),
           "tilt": round(m["tilt_max"], 1)}
    rows.append(row)
    print("RESULT", row, flush=True)

with open("_deaff_matrix_out.json", "w", encoding="utf-8") as f:
    json.dump(rows, f, indent=1)
print("matrix complete -> _deaff_matrix_out.json")
