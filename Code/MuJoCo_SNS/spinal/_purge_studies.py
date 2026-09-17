"""Dump curriculum trial histories, then delete contaminated studies.

Pre-lamination curriculum trials: s2 trials 0-24, s3 trials 0-29
(direct-synapse architecture); later trials are laminated but under-tuned.
Deleting: the 3 curr_* studies (internally MIXED architectures = TPE-
contaminated; _curriculum.py would RESUME them via load_if_exists) +
v7/v8/v9/v10 (tuned against the pre-lamination wiring).
Kept: v1-v6 single-architecture history (not in Ben's delete list,
irreversible, not resumable by the curriculum).
"""
import json

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
DB = "sqlite:///optuna_walk.db"

DUMP = ["curr_s1_air_deaff", "curr_s2_air_aff", "curr_s3_ground"]
hist = {}
for name in DUMP:
    try:
        st = optuna.load_study(study_name=name, storage=DB)
    except KeyError:
        continue
    hist[name] = [
        {"number": t.number, "value": t.value, "params": t.params}
        for t in st.trials
    ]
with open("curriculum_prelam_history_20260916.json", "w",
          encoding="utf-8") as f:
    json.dump(hist, f, indent=1)
print("dumped", {k: len(v) for k, v in hist.items()})

DELETE = [
    "curr_s1_air_deaff",
    "curr_s2_air_aff",
    "curr_s3_ground",
    "ground_walk_v7_rom",
    "ground_walk_v8_pose",
    "ground_walk_v8b_normal",   # absent if not found
    "ground_walk_v9_flat",
    "ground_walk_v10_transient",
    "ground_walk_v11",          # absent if not found
]
for name in DELETE:
    try:
        optuna.delete_study(study_name=name, storage=DB)
        print("deleted:", name)
    except KeyError:
        print("not present (skip):", name)

print("--- remaining ---")
for s in optuna.get_all_study_summaries(DB):
    best = f"{s.best_trial.value:.3f}" if s.best_trial is not None else "-"
    print(f"{s.study_name:35s} trials={s.n_trials:4d} best={best}")
