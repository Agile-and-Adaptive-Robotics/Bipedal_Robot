# Recon: enumerate optuna studies + winners on easteregg2 (read-only).
# Run from repo root: D:\Anaconda\envs\myo\python.exe <this> in spinal\
import json, os, sys, glob

os.environ.setdefault("OPTUNA_WARN_MULTIPLE_STUDY_OBJECT", "1")
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

print("=== DB CENSUS ===")
import sqlite3
for db in sorted(glob.glob("*.db")):
    try:
        con = sqlite3.connect(db)
        names = [r[0] for r in con.execute(
            "SELECT study_name FROM studies ORDER BY study_id").fetchall()]
        con.close()
        print(f"\n--- {db}: {len(names)} studies ---")
        for name in names:
            try:
                study = optuna.load_study(study_name=name, storage=f"sqlite:///{db}")
                try:
                    best = study.best_trial
                except Exception as e:
                    print(f"  {name}: n={len(study.trials)} NO_BEST ({type(e).__name__})")
                    continue
                p = best.params
                pshort = {k: (round(v, 4) if isinstance(v, float) else v)
                          for k, v in list(p.items())[:6]}
                print(f"  {name}: n={len(study.trials)} best_trial={best.number} "
                      f"value={best.value!r} params({len(p)})={pshort}")
            except Exception as e:
                print(f"  {name}: LOAD_FAIL {type(e).__name__}: {e}")
    except Exception as e:
        print(f"\n--- {db}: DB_FAIL {type(e).__name__}: {e}")

print("\n=== WINNER JSONS ===")
for f in sorted(glob.glob("decomp_winner_*.json")) + ["GAIT_DECOMPOSITION_20261006.json"]:
    try:
        d = json.load(open(f))
        print(f"\n--- {f} ---")
        s = json.dumps(d, indent=1)
        print(s[:1800])
        if len(s) > 1800:
            print(f"...[truncated, {len(s)} chars total]")
    except Exception as e:
        print(f"\n--- {f}: READ_FAIL {type(e).__name__}: {e}")

print("\n=== BEST*_PARAMS SCORE FIELDS ===")
for f in sorted(glob.glob("best_walk_params*.json")):
    try:
        d = json.load(open(f))
        keys = {k: d[k] for k in d if not isinstance(d[k], (dict, list))}
        print(f"{f}: top-level scalars={keys}")
    except Exception as e:
        print(f"{f}: READ_FAIL {e}")

print("\n=== CAMPAIGNS DIR ===")
for p in sorted(glob.glob("campaigns/*")):
    print(" ", p, "(dir)" if os.path.isdir(p) else f"({os.path.getsize(p)} B)")
