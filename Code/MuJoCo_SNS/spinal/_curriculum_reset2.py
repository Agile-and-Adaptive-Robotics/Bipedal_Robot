"""Purge the EXPLOITED curriculum studies (2026-09-20).

curr_s2_air_aff (rerun vintage, 25 trials): won by the static
deep-flexion pose (best 36.906 = 0.5*(-knee_min) with rises=0 - the
rhythm metric had no gate). curr_s3_ground (50 trials, all -100) seeded
from that non-rhythmic winner. Both archived here then deleted.
Stage 1 (curr_s1_air_deaff, 137.254) is VALID and untouched.
Reruns use fresh names curr_s2b_air_aff / curr_s3b_ground.
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

DB = "sqlite:///optuna_walk.db"
ARCHIVE = "curriculum_exploit_archive_20260920.json"

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
