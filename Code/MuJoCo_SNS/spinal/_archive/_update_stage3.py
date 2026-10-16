"""Point curriculum_stage3.json at a study's best trial (arg1 = study
name, arg2 = trial number) so the next stage seeds from it."""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

study = sys.argv[1] if len(sys.argv) > 1 else "curr_s3d_ground"
trial = int(sys.argv[2]) if len(sys.argv) > 2 else 33

s = optuna.load_study(study_name=study, storage="sqlite:///optuna_walk.db")
t = s.trials[trial]
out = {"stage": 3, "score": t.value, "trial": trial, "study": study,
       "params": t.params}
with open("curriculum_stage3.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=2)
print(f"curriculum_stage3.json -> {study} t{trial} value {t.value:.3f}; "
      f"pf {t.params.get('pf_gain'):.2f} cs "
      f"{t.params.get('contra_swing'):.2f}")
