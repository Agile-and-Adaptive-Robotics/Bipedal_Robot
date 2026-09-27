"""Re-measure the s3k stock-chain winner ground eval under the FIXED
kine_ref (2026-09-26 left-cycle fix). Faithful re-run of
_pf_layer_variant_test.py gate 3 (the documented
-160.23425729850192 baseline) with TWO adaptations forced by repo
state:
  1) C.set_stage(4, ...) instead of the stale set_stage(3, ...):
     the 2026-09-24-night reorder (Ben: "standing should be stage 3
     and walking stage 4") moved the ground-walk key set from stage 3
     to stage 4; old stage-3 block == new stage-4 block, verified
     line-by-line vs commit 2b1d679a.
  2) AARL_NPZ=scratch_s3k_rescore.npz so the protected spinal_run.npz
     is never touched (runner.py:1629 AARL_NPZ redirect).
"""
import json
import os
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).parent          # reports_20260925/tmp
SPINAL = HERE.parent.parent           # .../spinal
sys.path.insert(0, str(SPINAL))
os.chdir(SPINAL)

assert "AARL_NET" not in os.environ, "run with AARL_NET unset"
os.environ["AARL_NPZ"] = "scratch_s3k_rescore.npz"
print(f"AARL_NPZ = {os.environ['AARL_NPZ']} (spinal_run.npz untouched)")

import muscle_map as MM  # noqa: F401  (same import set as the gate script)
import params as P
import build_network as BN  # noqa: F401
import _curriculum as C
import runner as R

mul = json.loads((SPINAL / "best_walk_params_v10.json")
                 .read_text(encoding="utf-8"))["multipliers"]
st3 = json.loads((SPINAL / "curriculum_stage3.json")
                 .read_text(encoding="utf-8"))["params"]
C.BASE_MUL = dict(mul)
C.BASE_MUL["renshaw"] = 0.5
C.set_stage(4, {**mul, **st3, "full_rules": 0.0})
P.G["full_rules"] = 0.0
print(f"drive={st3['drive']!r} full_rules={P.G['full_rules']!r} "
      f"renshaw={P.G['renshaw']!r} heel_rge={P.G['heel_rge']!r} "
      f"ib_rge={P.G['ib_rge']!r} no_cross_applied="
      f"{P.G['contra_swing'] == 0.0 and P.G['v3_gain'] == 0.0}")
m = R.main(["--eval", "--drive", repr(st3["drive"])])
ks = m["kine_score"]
print(f"NEW kine_score (fixed left ref) = {ks!r}")
print("old (buggy left ref)           = -160.23425729850192")
print(f"delta = {ks - (-160.23425729850192):+.6f}")
