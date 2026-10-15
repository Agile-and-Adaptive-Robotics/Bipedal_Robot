"""Point curriculum_stage3.json at the CORRECTED best (s3c trial 54)
so the s3d seed chains from it. Keeps the study-record fields."""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

s = optuna.load_study(study_name="curr_s3c_ground",
                      storage="sqlite:///optuna_walk.db")
t54 = s.trials[54]
out = {
    "stage": 3,
    "score": -181.58612797682238,   # corrected kine_ref v2 (frozen-leg
    "trial": 54,                    # hole closed); study record was
    "study": "curr_s3c_ground",     # t26 -119.3 under the holed v2
    "note": "corrected-objective best (frozen-leg hole fixed); "
            "study's recorded best under holed score was t26 -119.3",
    "params": t54.params,
}
with open("curriculum_stage3.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print("curriculum_stage3.json -> t54",
      json.dumps({k: round(v, 3) for k, v in t54.params.items()}))
