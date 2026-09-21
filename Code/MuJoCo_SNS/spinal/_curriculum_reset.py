"""Purge the all-sentinel curriculum studies (Ben GO 2026-09-18).

Dumps every trial of curr_s2_air_aff / curr_s3_ground to a dated archive
json, then deletes the studies so the rerun starts fresh.  Stage 1
(curr_s1_air_deaff, score 137.254) is VALID and untouched.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

DB = "sqlite:///optuna_walk.db"
ARCHIVE = "curriculum_sentinel_archive_20260918.json"

out = {}
for name in ("curr_s2_air_aff", "curr_s3_ground"):
    try:
        study = optuna.load_study(study_name=name, storage=DB)
    except KeyError:
        print(f"{name}: not present, nothing to purge")
        continue
    trials = [{"number": t.number, "value": t.value, "params": t.params}
              for t in study.trials]
    out[name] = trials
    print(f"{name}: archived {len(trials)} trials "
          f"(best {max((t['value'] or -1e9) for t in trials)})")
    optuna.delete_study(study_name=name, storage=DB)
    print(f"{name}: DELETED from db")

if out:
    with open(ARCHIVE, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print(f"archive: {ARCHIVE}")
print("purge complete; stage-1 study untouched")
