"""Full-literature-connectome test (2026-09-21, Ben: "wire it the way
the working Deng model does"). s3i winner config with full_rules=1:
gate 1 = full_rules 0 bit-identity vs -148.6878643; gate 2 = ON:
air rhythm smoke + ground both-leg readout."""
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
import params as P
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
os.environ["AARL_NPZ"] = "spinal_run_fr.npz"


def run(fr, air=False):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    C.set_stage(3, {**mul, **st3, "full_rules": fr})
    P.G["full_rules"] = fr
    args = (["--no-ground", "--no-afferents", "--time", "10"]
            if air else ["--eval", "--drive", repr(st3["drive"])])
    m = R.main(args)
    z = np.load("spinal_run_fr.npz", allow_pickle=True)
    t, q, neuro, contact = z["t"], z["q"], z["neuro"], z["contact"]
    names = [str(s) for s in z["neuro_names"]]
    if air:
        w = (t >= 3.0) & (t <= 10.0)
        rge = neuro[w, names.index("RG_E_r")]
        rgl = neuro[w, names.index("RG_E_l")]
        on = (rge > 0.5 * max(rge.max(), 1e-9)).astype(int)
        rises = int(np.sum(np.diff(on) == 1))
        x = float(np.corrcoef(rge - rge.mean(),
                              rgl - rgl.mean())[0, 1])
        print(f"  AIR: finite {bool(np.all(np.isfinite(neuro[w])))} "
              f"RG_E_r swings {rge.max()-rge.min():.2f} mV, bursts "
              f"{rises}, R/L corr {x:+.2f}")
    else:
        k = KR.compare(t, q, neuro, 2.0, ref=KR.ref_cached(),
                       contact=contact)
        ks = None if k is None else round(k["kine_score"], 1)
        print(f"  GROUND: score {ks} | duty r/l {k.get('duty_r'):.2f}/"
              f"{k.get('duty_l')} | cyc r/l "
              f"{k.get('n_cycles_r')}/{k.get('n_cycles_l')} | frozen_l "
              f"{k.get('frozen_l')} | nan {m['nan']}")


print("=== gate 1: full_rules 0 bit-identity ===")
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(3, {**mul, **st3, "full_rules": 0.0})
P.G["full_rules"] = 0.0
m = R.main(["--eval", "--drive", repr(st3["drive"])])
ok = abs(m["kine_score"] - (-148.6878643)) < 1e-6
print("GATE 1 (bit-identity):", "PASS" if ok else f"FAIL ({m['kine_score']})")
print("=== gate 2: full_rules 1 - air rhythm ===")
run(1.0, air=True)
print("=== gate 3: full_rules 1 - ground ===")
run(1.0, air=False)
