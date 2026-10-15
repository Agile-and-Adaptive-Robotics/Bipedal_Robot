"""Smoke gates for the phase machine v2 weight-shift (2026-09-21).
Gate 1: pm_ws=0 on the s3f t32 config must reproduce -164.4258... exactly
(pm_ws=0 keeps v1 behavior: no gate-hold, no abductor prep).
Gate 2: pm_ws=0.5 finite + left-leg readout."""
import json
import os
import sys
from pathlib import Path

import numpy as np

import _curriculum as C
import kine_ref as KR
import params as P
import runner as R
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
s = optuna.load_study(study_name="curr_s3f_ground",
                      storage="sqlite:///optuna_walk.db")
REF = s.trials[32].value
print(f"s3f t32 exact: {REF!r}")
os.environ["AARL_NPZ"] = "spinal_run_smoke5.npz"


def run(ws):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    p = {**mul, **st3, "pm_ws": ws}
    C.set_stage(3, p)
    print(f"pm_ws = {P.G['pm_ws']}")
    m = R.main(["--eval", "--drive", repr(st3["drive"])])
    z = np.load("spinal_run_smoke5.npz", allow_pickle=True)
    k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                   ref=KR.ref_cached(), contact=z["contact"])
    z.close()
    ks = float(m["kine_score"])
    if k:
        print(f"  score {ks:.10f} | duty r/l {k.get('duty_r')} / "
              f"{k.get('duty_l')} | cyc r/l "
              f"{k.get('n_cycles_r')}/{k.get('n_cycles_l')} | "
              f"frozen_l {k.get('frozen_l')} | nan {m['nan']}")
    else:
        print(f"  score {ks} kine None nan {m['nan']}")
    return ks, (k or {})


print("=== gate 1: pm_ws 0 bit-identity ===")
k0, _ = run(0.0)
ok = abs(k0 - REF) < 1e-9
print("GATE 1:", "PASS" if ok else f"FAIL ({k0})")
print("=== gate 2: pm_ws 0.5 ===")
k1, d1 = run(0.5)
print("GATE 2 (finite):", "PASS" if np.isfinite(k1) else "FAIL")
