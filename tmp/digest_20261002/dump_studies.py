"""Dump all study results from the optuna dbs used by the overnight curriculum."""
import optuna
import sqlite3
from pathlib import Path

SPINAL = Path(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
DBS = [
    ("optuna_walk.db", SPINAL),
    ("optuna_syn6.db", SPINAL),
    ("optuna_w2lvar.db", SPINAL),
    ("optuna_pruned_pruned_air_fix_20261003.db", SPINAL),
    ("optuna_pruned_pruned_walk_1_fix_20261003.db", SPINAL),
    ("optuna_pruned_pruned_walk_2_fix_20261003.db", SPINAL),
]
# spiking + bpa dbs might be elsewhere
for pat in ["optuna_spk*.db", "optuna_bpa*.db"]:
    for p in SPINAL.glob(pat):
        DBS.append((p.name, SPINAL))
    for p in (SPINAL / "w2l_mujoco").glob(pat):
        DBS.append((p.name, p / "w2l_mujoco" if False else SPINAL / "w2l_mujoco"))

seen = set()
for dbname, cwd in DBS:
    dbpath = cwd / dbname
    if not dbpath.exists() or dbname in seen:
        continue
    seen.add(dbname)
    try:
        study_names = optuna.get_all_study_names(storage=f"sqlite:///{dbpath}")
    except Exception:
        continue
    for sname in study_names:
        if "20261003" not in sname and "fix" not in sname.lower():
            continue
        try:
            study = optuna.load_study(study_name=sname, storage=f"sqlite:///{dbpath}")
            trials = study.trials
            complete = [t for t in trials if t.state.name == "COMPLETE"]
            if not complete:
                continue
            best = min(complete, key=lambda t: t.value)
            print(f"\n=== {dbname} :: {sname} ===")
            print(f"    trials: {len(trials)} total / {len(complete)} complete")
            print(f"    best:   {best.value:.4f} (trial {best.number})")
            # top 5
            top5 = sorted(complete, key=lambda t: t.value)[:5]
            for t in top5:
                print(f"    t{t.number:3d}: {t.value:.4f}")
            # key params of the best
            if best.params:
                ps = {k: round(v, 4) if isinstance(v, float) else v for k, v in best.params.items()}
                print(f"    params: {ps}")
        except Exception as e:
            print(f"    ERROR: {e}")
