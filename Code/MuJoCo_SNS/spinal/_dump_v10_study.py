"""Dump the study-path (set_params) state for the v10 winner."""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna_walk_v10 as OW
import params as P

prev = json.loads(open("best_walk_params_v10.json", encoding="utf-8").read())
mul = prev["multipliers"]
OW.load_fitted_baseline()
OW.set_params(mul)
state = dict(
    W_PF_MN=P.W_PF_MN, W_POSTURE=P.W_POSTURE, TAU=P.TAU, G=P.G,
    PF_SHAPE=P.PF_SHAPE, BAL=P.BAL, AFF=P.AFF, MOD=P.MOD,
    walk_drive=float(prev["params"]["drive"]),
    SCHEDULE={k: list(v) for k, v in P.SCHEDULE.items()},
)
with open(os.environ["RUNNER_DUMP_STATE"], "w", encoding="utf-8") as f:
    json.dump(state, f, indent=1, default=float)
print("dumped study-path state")
