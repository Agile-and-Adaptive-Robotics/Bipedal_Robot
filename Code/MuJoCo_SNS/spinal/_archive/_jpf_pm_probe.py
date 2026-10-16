"""Plan-B probe: s3f t32 config + joint-layer PF (joint_pf=1) on/off,
with the phase machine active. Does the HIP/KNEE/ANK half-center layer
change the frozen-left picture? Side npz; safe alongside s3g."""
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
os.environ["AARL_NPZ"] = "spinal_run_jpf.npz"

for jpf in (0.0, 1.0):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    C.set_stage(3, {**mul, **st3})
    P.G["joint_pf"] = jpf
    m = R.main(["--eval", "--joint-pf", "1" if jpf else "0",
                "--drive", repr(st3["drive"])])
    z = np.load("spinal_run_jpf.npz", allow_pickle=True)
    k = KR.compare(z["t"], z["q"], z["neuro"], 2.0,
                   ref=KR.ref_cached(), contact=z["contact"])
    z.close()
    if k:
        print(f"[joint_pf {jpf:.0f}] score {k['kine_score']:.1f} | "
              f"duty r/l {k.get('duty_r'):.2f}/{k.get('duty_l')} | "
              f"cyc r/l {k.get('n_cycles_r')}/{k.get('n_cycles_l')} | "
              f"frozen_l {k.get('frozen_l')} | nan {m['nan']}",
              flush=True)
    else:
        print(f"[joint_pf {jpf:.0f}] NO KINE (nan {m['nan']})",
              flush=True)
