"""Read-only peek at curriculum study progress (safe vs the running
optimizer: sqlite read of completed trials)."""
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
for name in ("curr_s1_air_deaff", "curr_s2_air_aff", "curr_s3_ground"):
    try:
        st = optuna.load_study(study_name=name,
                               storage="sqlite:///optuna_walk.db")
    except KeyError:
        print(f"{name}: not created yet")
        continue
    done = [t for t in st.trials
            if t.state == optuna.trial.TrialState.COMPLETE
            and t.value is not None]
    if not done:
        print(f"{name}: 0 complete trials")
        continue
    best = max(done, key=lambda t: t.value)
    print(f"{name}: {len(st.trials)} trials ({len(done)} complete), "
          f"best {best.value:.3f} @ trial {best.number}")
