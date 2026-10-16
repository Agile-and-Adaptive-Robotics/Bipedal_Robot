"""Smoke gates for contra_kinh (2026-09-21).
Gate 1 (bit-identity): s3d-best (t33) config with contra_kinh=0 must
reproduce the s3d study value EXACTLY (-191.9877...).
Gate 2: contra_kinh=0.5 finite; print the LEFT-leg verdict."""
import json
import os
import sys
from pathlib import Path

import numpy as np

import _curriculum as C
import kine_ref as KR
import params as P
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
os.environ["AARL_NPZ"] = "spinal_run_smoke3.npz"


def run(ck):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    p = {**mul, **st3, "contra_kinh": ck}
    C.set_stage(3, p)
    print(f"contra_kinh = {P.G['contra_kinh']}")
    m = R.main(["--eval", "--drive", repr(st3["drive"])])
    z = np.load("spinal_run_smoke3.npz", allow_pickle=True)
    k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                   ref=KR.ref_cached(), contact=z["contact"])
    z.close()
    ks = float(m["kine_score"])
    if k:
        print(f"  score {ks:.6f} duty_r {k.get('duty_r')} cyc r/l "
              f"{k.get('n_cycles_r')}/{k.get('n_cycles_l')} "
              f"frozen_l {k.get('frozen_l')} nan {m['nan']}")
    else:
        print(f"  score {ks:.6f} kine None nan {m['nan']}")
    return ks, (k or {})


print("=== gate 1: contra_kinh 0 must be bit-identical to s3d t33 ===")
k0, d0 = run(0.0)
ok = abs(k0 - (-191.98767443280454)) < 1e-9
print("GATE 1 (bit-identity):", "PASS" if ok else f"FAIL ({k0})")
print("=== gate 2: contra_kinh 0.5 ===")
k1, d1 = run(0.5)
print("GATE 2 (finite):", "PASS" if np.isfinite(k1) else "FAIL")
