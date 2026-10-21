"""Stage-3 balance audit pull (2026-10-20): dump the state of the
balance studies in the three optuna DBs so the winners can be re-run
fresh through the standing eval. Read-only; no studies created.

Usage: python _s3_pull_20261020.py > _s3_pull_20261020_out.txt
"""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

TARGETS = [
    ("sqlite:///optuna_walk.db", ["curr_s3_balance", "curr_s4_balance"]),
    ("sqlite:///optuna_w2lvar.db", ["curr_w2lvar_s3_balance",
                                    "curr_w2lvar_s3_fix_20261003"]),
    ("sqlite:///optuna_syn6.db", ["curr_syn6_s3_balance"]),
]

out = {}
for db, names in TARGETS:
    print(f"\n===== {db} =====")
    try:
        all_names = list(optuna.get_all_study_names(storage=db))
    except Exception as e:
        print(f"  (db unreadable: {e})")
        continue
    print(f"  all studies: {all_names}")
    for name in names:
        if name not in all_names:
            print(f"  -- {name}: DOES NOT EXIST (never launched)")
            out[f"{db}::{name}"] = None
            continue
        st = optuna.load_study(study_name=name, storage=db)
        trs = st.trials
        finished = [t for t in trs
                    if t.state == optuna.trial.TrialState.COMPLETE]
        print(f"  -- {name}: {len(trs)} trials "
              f"({len(finished)} complete)")
        vals = [t.value for t in finished]
        if vals:
            import collections
            hist = collections.Counter(
                "fell" if v == -150.0 else ("nan" if v == -200.0
                                            else "stand")
                for v in vals)
            print(f"     value histogram: {dict(hist)}")
            best = st.best_trial
            print(f"     BEST trial {best.number} value "
                  f"{best.value:.4f}")
            print(f"     params: {json.dumps(best.params)}")
            # top-5 for context (are there ANY non-sentinel trials?)
            top5 = sorted(finished, key=lambda t: -t.value)[:5]
            for t in top5:
                print(f"       t{t.number}: {t.value:.4f}")
            out[f"{db}::{name}"] = dict(
                n=len(trs), best_value=best.value,
                best_number=best.number, best_params=best.params,
                hist={k: int(v) for k, v in hist.items()})

with open("_s3_pull_20261020.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)
print("\nwrote _s3_pull_20261020.json")
