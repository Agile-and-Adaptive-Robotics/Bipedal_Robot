"""Smoke gates for contra_swing (2026-09-21).
Gate 1 (bit-identity): t54 config with contra_swing=0 must reproduce
the t54 eval score EXACTLY (-181.58612797682238).
Gate 2 (finite + effect): contra_swing=0.5 finite, metrics differ.
NOTE: no stdout wrapper here - _curriculum wraps on import (closing a
pre-existing wrapper; the documented double-wrap trap)."""
import json
import sys
from pathlib import Path

import _curriculum as C
import params as P
import runner as R

HERE = Path(__file__).parent
mul = json.loads((HERE / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((HERE / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]


def run(cs):
    C.BASE_MUL = dict(mul)
    C.BASE_MUL["renshaw"] = 0.5
    p = {**mul, **st3, "contra_swing": cs}
    C.set_stage(3, p)
    print(f"contra_swing = {P.G['contra_swing']}")
    m = R.main(["--eval", "--drive", repr(st3["drive"])])
    k = m.get("kine")
    print(f"  kine_score {m['kine_score']:.17g}  nan {m['nan']} "
          f"kz {m['kz']:.4f}")
    if k:
        print(f"  duty_r {k.get('duty_r')} cyc_r {k.get('n_cycles_r')} "
              f"cyc_l {k.get('n_cycles_l')} "
              f"frozen_l {k.get('frozen_l')}")
    return m


print("=== gate 1: contra_swing 0 must be bit-identical to t54 ===")
m0 = run(0.0)
ok1 = abs(m0["kine_score"] - (-181.58612797682238)) < 1e-9
print("GATE 1 (bit-identity):", "PASS" if ok1 else "FAIL")
print("=== gate 2: contra_swing 0.5 finite ===")
m1 = run(0.5)
print("GATE 2 (finite):", "PASS" if not m1["nan"] else "FAIL")
