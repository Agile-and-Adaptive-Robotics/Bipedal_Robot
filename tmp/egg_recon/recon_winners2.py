# Recon part 2: db census (fixed) + report-folder score mining. Run in spinal\.
import os, re, glob

here = os.path.dirname(os.path.abspath(__file__))
os.chdir(here)

import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)

print("=== DB CENSUS (best per study) ===")
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
                print(f"  {name}: n={len(study.trials)} best_trial={best.number} value={best.value!r}")
            except Exception as e:
                print(f"  {name}: LOAD_FAIL {type(e).__name__}: {e}")
    except Exception as e:
        print(f"\n--- {db}: DB_FAIL {type(e).__name__}: {e}")

print("\n=== REPORT DIRS ===")
for d in ["reports_20260923", "reports_20260924", "reports_20260925",
          "reports_spiking_20261002", "campaigns/20260930"]:
    for p in sorted(glob.glob(d + "/*")):
        if os.path.isfile(p):
            print(f"  {p}  ({os.path.getsize(p)} B)")

print("\n=== SCORE GREP (lines citing kine/lib scores) ===")
pat = re.compile(r"(-?18[5-9]\.\d|-?19\d\.\d|-?20[0-9]\.\d|of-record|gait_lib|s3k|syn6|pruned)", re.I)
for d in ["reports_20260923", "reports_20260924", "reports_20260925",
          "reports_spiking_20261002", "campaigns/20260930"]:
    for p in sorted(glob.glob(d + "/**/*.md", recursive=True)):
        try:
            txt = open(p, encoding="utf-8", errors="replace").read()
        except Exception:
            continue
        hits = [ln.strip() for ln in txt.splitlines() if pat.search(ln)]
        if hits:
            print(f"\n--- {p} ({len(hits)} hit lines) ---")
            for h in hits[:25]:
                print("   ", h[:220].encode("ascii", "replace").decode("ascii"))
