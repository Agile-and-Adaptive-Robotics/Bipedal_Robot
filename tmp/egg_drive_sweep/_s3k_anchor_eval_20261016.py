"""Reproduce the s3k anchor evals bit-exactly (campaign 2026-10-16).

FLAVOR A (the recorded baseline anchor): the AGENTS.md 2026-09-26 recipe
  set_stage(4, {**v10 multipliers, **curriculum_stage3.json params,
                full_rules: 0}) + runner --eval at the winner drive
  expected kine_score = -161.56754173676563 (recorded post-left-ref-fix).
FLAVOR B (production topology, full_rules=1 as curr_s3k_nocross ran it):
  same merge WITHOUT the full_rules override; no recorded number - this
  run records it (pre-change observation 2026-10-16: -25.0).

Run with PYTHONIOENCODING=utf-8. Run in separate processes (params.G is
module-global; one flavor per process keeps the merge unambiguous).
"""
import json
import os
import sys

flavor = sys.argv[1] if len(sys.argv) > 1 else "A"

import _curriculum as C
import runner as R

base = dict(json.loads(
    open("best_walk_params_v10.json", encoding="utf-8").read())
    ["multipliers"])
base["renshaw"] = 0.5
s3k = json.loads(open(
    "reports_20260923/s3k_trial34_full_params.json",
    encoding="utf-8").read())["params"]
base.update(s3k)
if flavor == "A":
    base["full_rules"] = 0.0     # the AGENTS 2026-09-26 anchor recipe
C.set_stage(4, base)
m = R.main(["--eval", "--drive", repr(base["drive"])])
print(f"ANCHOR[{flavor}] kine_score:", repr(m.get("kine_score")))
if flavor == "A":
    print("ANCHOR[A] expected  : -161.56754173676563")
    print("ANCHOR[A] BIT-EXACT:",
          repr(m.get("kine_score")) == "-161.56754173676563")
