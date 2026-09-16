"""Print the top-N completed trials of a study with selected params.

Usage: python _study_top.py <study_name> [n]
"""
import sys

import optuna

study = optuna.load_study(study_name=sys.argv[1],
                          storage="sqlite:///optuna_walk.db")
n = int(sys.argv[2]) if len(sys.argv) > 2 else 8
done = [t for t in study.trials
        if t.state == optuna.trial.TrialState.COMPLETE]
done.sort(key=lambda t: t.value, reverse=True)
print(f"{study.study_name}: {len(done)} complete, best {study.best_value:.3f}"
      f" (trial {study.best_trial.number})")
keys = ["phase_reset_e", "phase_reset_f", "f1_kneext_inh", "drive",
        "rg_adapt", "desc_e", "desc_f", "pf_gain", "rg_to_pf", "e2_adapt"]
head = "trial    score   " + "  ".join(f"{k[:9]:>9}" for k in keys)
print(head)
for t in done[:n]:
    row = f"{t.number:5d} {t.value:8.3f}  " + "  ".join(
        f"{t.params.get(k, float('nan')):9.3f}" for k in keys)
    print(row)
