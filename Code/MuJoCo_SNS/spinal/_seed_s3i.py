"""Reseed curriculum_stage3.json from s3h t39 + pm_aff 0.8."""
import json
import sys

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
s = optuna.load_study(study_name="curr_s3h_ground",
                      storage="sqlite:///optuna_walk.db")
t = s.trials[39]
p = dict(t.params)
p["pm_aff"] = 0.8
out = {"stage": 3, "score": t.value, "trial": 39,
       "study": "curr_s3h_ground", "params": p}
json.dump(out, open("curriculum_stage3.json", "w"), indent=2)
print("stage3.json = s3h t39 + pm_aff 0.8;",
      "pm_add", round(p["pm_add"], 2), "pm", round(p["pm_gain"], 2),
      "ky", round(p["ky_scale"], 3))
