"""query_best.py <db> <study_name>  -> prints study best_value + n_trials.

Tiny helper for the w2lvar tuner workflow (STEP 3/4 batch bookkeeping).
Read-only: opens the sqlite storage, prints, exits. Counted-trial census
(COMPLETE states) is included so sentinel plateaus are visible.
"""
import sys

import optuna

db, study_name = sys.argv[1], sys.argv[2]
st = optuna.load_study(study_name=study_name,
                       storage=f"sqlite:///{db}")
vals = [t.value for t in st.trials]
comp = sum(1 for t in st.trials if t.state.name == "COMPLETE")
print(f"study={study_name} n_trials={len(vals)} COMPLETE={comp}")
if st.best_trial is not None:
    print(f"best_value={st.best_value:.6f} (trial {st.best_trial.number})")
else:
    print("best_value=None (no completed trials)")
print("values:", [None if v is None else round(v, 3) for v in vals])
