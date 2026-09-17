"""Isolate the ground rhythm death: stage-1 winner config under
(a) air + afferents ON (ground off)
(b) ground full run (22 s, watch the joints + RG cycles)
(c) ground + --no-afferents (contacts but no afferent encoding)
"""
import json

import _curriculum as C
import runner as R

s1 = json.loads(open("curriculum_stage1.json", encoding="utf-8").read())
mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]
p = {**mul, **s1["params"]}

print("=== (a) AIR + AFFERENTS (ground off) ===")
C.set_stage(1, p)
R.main(["--no-ground", "--no-interleg", "--time", "14",
        "--drive", repr(p["drive"])])

print("\n=== (b) GROUND full run ===")
R.main(["--time", "22", "--drive", repr(p["drive"])])

print("\n=== (c) GROUND + no afferents ===")
R.main(["--no-afferents", "--time", "22", "--drive", repr(p["drive"])])
