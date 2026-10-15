"""Smoke gates for the contact-onset change (2026-09-20).

Gate 1 (bit-identity): stage-3 config with contact_onset=0 must
reproduce the pre-change diag metrics BIT-EXACTLY
(dx 0.01837978106322702, kz 0.8940922282853601, duty 1.0).
Gate 2 (finite + effect): contact_onset=0.6 must stay finite and be
seen in the run (report metrics; qualitative readout only).
"""
import io
import json
import sys
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import _curriculum as C
import params as P
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]


def run(onset):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    p = {**mul, **st3, "contact_onset": onset}
    C.set_stage(3, p)
    print(f"contact_onset = {P.G['contact_onset']}")
    m = R.main(["--eval", "--drive", repr(st3["drive"])])
    print(f"  -> kz {m['kz']:.17g}  duty {m['duty']}  burst_r "
          f"{m['burst_r']}  dx {m['dx']:.17g}  nan {m['nan']}")
    return m


print("=== gate 1: onset 0 must be bit-identical ===")
m0 = run(0.0)
ok = (abs(m0["dx"] - 0.01837978106322702) < 1e-15
      and abs(m0["kz"] - 0.8940922282853601) < 1e-15
      and m0["duty"] == 1.0)
print("GATE 1 (bit-identity):", "PASS" if ok else "FAIL")

print("=== gate 2: onset 0.6 finite ===")
m1 = run(0.6)
z = np.load("spinal_run.npz", allow_pickle=True)
fin = bool(np.all(np.isfinite(z["q"])) and np.all(np.isfinite(z["neuro"])))
print("GATE 2 (finite):", "PASS" if (fin and not m1["nan"]) else "FAIL")
print(f"  duty {m1['duty']}  burst_r {m1['burst_r']}  "
      f"kine_score {m1.get('kine_score')}")
