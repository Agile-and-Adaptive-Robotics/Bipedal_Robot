"""List optuna studies in optuna_walk.db (name, n trials, best value)."""
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
summary = optuna.get_all_study_summaries("sqlite:///optuna_walk.db")
for s in summary:
    best = "-"
    if s.best_trial is not None:
        best = f"{s.best_trial.value:.3f}"
    print(f"{s.study_name:35s} trials={s.n_trials:4d} best={best}")
