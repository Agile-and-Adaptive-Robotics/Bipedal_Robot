"""Dump the set_params-path param state (study path) for diffing."""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna_walk as OW
import runner as R

prev = json.loads(open("best_walk_params.json", encoding="utf-8").read())
mul = prev["multipliers"]
OW.load_fitted_baseline()
OW.set_params(mul)

# reuse runner's dump block by invoking main with a tiny schedule? No:
# dump directly with the same structure.
js = dict(
    W_PF_MN={ph: dict(t) for ph, t in OW.FIT["W_PF_MN"].items()},
)
import params as P
state = dict(
    W_PF_MN=P.W_PF_MN, W_POSTURE=P.W_POSTURE, TAU=P.TAU, G=P.G,
    PF_SHAPE=P.PF_SHAPE, BAL=P.BAL, AFF=P.AFF, MOD=P.MOD,
    walk_drive=float(f"{mul['drive']:.4f}"),
    SCHEDULE={k: list(v) for k, v in P.SCHEDULE.items()},
)
with open(os.environ["RUNNER_DUMP_STATE"], "w", encoding="utf-8") as f:
    json.dump(state, f, indent=1, default=float)
print("dumped study-path state")
