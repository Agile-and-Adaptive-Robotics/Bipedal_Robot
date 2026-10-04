"""Read-only peek at the optuna DBs used by the tuning-seed round.

Lists study names + trial counts so fresh dated study names can be
checked for collisions before launch. Read-only: load_study only,
no optimize, no writes (per the python-skill sqlite hygiene note).
"""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
for db in ("optuna_walk.db", "optuna_syn6.db", "optuna_w2lvar.db"):
    print(f"== {db}")
    try:
        st = optuna.get_all_study_names(storage=f"sqlite:///{SPINAL}\\{db}")
    except Exception as e:  # noqa: BLE001 - report, do not force
        print(f"  ERROR: {type(e).__name__}: {e}")
        continue
    for name in st:
        try:
            s = optuna.load_study(study_name=name,
                                  storage=f"sqlite:///{SPINAL}\\{db}")
            n = len(s.trials)
            done = sum(1 for t in s.trials if t.state.name == "COMPLETE")
            best = None
            try:
                best = s.best_value
            except Exception:  # noqa: BLE001 - empty studies have no best
                pass
            print(f"  {name}: trials={n} complete={done} best={best}")
        except Exception as e:  # noqa: BLE001
            print(f"  {name}: ERROR {type(e).__name__}: {e}")
