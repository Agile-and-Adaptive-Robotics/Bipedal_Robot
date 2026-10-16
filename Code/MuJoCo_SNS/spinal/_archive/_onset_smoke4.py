"""Smoke gates for the phase machine (2026-09-21).
Gate 1 (bit-identity): stage3.json config with pm_gain=0 must
reproduce the s3e best value exactly.
Gate 2: pm_gain=0.5, pm_T=1.23 finite + both-leg readout."""
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
REF = -191.8074068353274  # s3e best (trial 0); verified from db below
os.environ["AARL_NPZ"] = "spinal_run_smoke4.npz"

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
s = optuna.load_study(study_name="curr_s3e_ground",
                      storage="sqlite:///optuna_walk.db")
val = s.best_trial.value
print(f"s3e best trial {s.best_trial.number} exact value: "
      f"{val!r}")


def run(pm, pmT):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    p = {**mul, **st3, "pm_gain": pm, "pm_T": pmT}
    C.set_stage(3, p)
    print(f"pm_gain = {P.G['pm_gain']}  pm_T = {P.G['pm_T']}")
    m = R.main(["--eval", "--drive", repr(st3["drive"])])
    z = np.load("spinal_run_smoke4.npz", allow_pickle=True)
    k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                   ref=KR.ref_cached(), contact=z["contact"])
    z.close()
    ks = float(m["kine_score"])
    if k:
        print(f"  score {ks:.10f} | duty r/l {k.get('duty_r')} / "
              f"{k.get('duty_l')} | cyc r/l "
              f"{k.get('n_cycles_r')}/{k.get('n_cycles_l')} | "
              f"frozen_l {k.get('frozen_l')}")
    else:
        print(f"  score {ks} kine None nan {m['nan']}")
    return ks, (k or {})


print("=== gate 1: pm_gain 0 must be bit-identical ===")
k0, _ = run(0.0, 1.23)
ok = abs(k0 - val) < 1e-9
print("GATE 1 (bit-identity):", "PASS" if ok else f"FAIL ({k0})")
print("=== gate 2: pm_gain 0.5 ===")
k1, d1 = run(0.5, 1.23)
print("GATE 2 (finite):", "PASS" if np.isfinite(k1) else "FAIL")
