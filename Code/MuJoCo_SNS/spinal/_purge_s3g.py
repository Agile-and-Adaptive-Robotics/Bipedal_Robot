"""Purge the 13-trial curr_s3g_ground run (superseded: it ran with the
rigid lateral rig; ky_scale + 15% gate added). Archive + delete."""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
DB = "sqlite:///optuna_walk.db"
name = "curr_s3g_ground"
study = optuna.load_study(study_name=name, storage=DB)
trials = [{"number": t.number, "value": t.value, "params": t.params}
          for t in study.trials]
with open("curriculum_s3g_rigidrig_archive_20260921.json", "w",
          encoding="utf-8") as f:
    json.dump({name: trials}, f, indent=1)
optuna.delete_study(study_name=name, storage=DB)
print(f"archived + deleted {name} ({len(trials)} trials, rigid-rig "
      f"vintage)")
