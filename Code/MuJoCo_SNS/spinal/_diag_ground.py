"""Diagnose the stage-2/3 -100 flatline: run the ground eval manually for
(a) the stage-1 winner at zero new gains, (b) the stage-3 winner params;
print the full metrics dict + a rhythm summary."""
import json

import _curriculum as C
import runner as R

s1 = json.loads(open("curriculum_stage1.json", encoding="utf-8").read())
s3 = json.loads(open("curriculum_stage3.json", encoding="utf-8").read())
mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]

for tag, st, params in (("stage1 winner, gains 0 (ground eval)", 3, s1["params"]),
                        ("stage3 winner (ground eval)", 3, s3["params"])):
    p = {**mul, **params}
    C.set_stage(st, p)
    print(f"\n=== {tag}: drive {p['drive']:.3f}")
    m = R.main(["--eval", "--drive", repr(p["drive"])])
    print(f"metrics: nan={m['nan']} t_end={m.get('t_end')} "
          f"kz={m.get('kz')} tilt_max={m.get('tilt_max')} "
          f"kine={'None' if m.get('kine') is None else 'present'} "
          f"kine_score={m.get('kine_score')}")
    if m.get("kine"):
        k = m["kine"]
        print(f"  n_cycles={k.get('n_cycles')} duty={k.get('duty')}")
