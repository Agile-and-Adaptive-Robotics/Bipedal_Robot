"""Status probe: state of the resume-target studies (2026-10-03 overnight campaign)."""
import optuna

TARGETS = [
    ("optuna_walk.db", "curr_s4_nocross_mutfix_20261003"),
    ("optuna_syn6.db", "curr_syn6_s5_reseed_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s1_fix_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s2_fix_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s3_fix_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s1_air_deaff"),
    ("optuna_w2lvar.db", "curr_w2lvar_s5_walk_contact"),
    ("optuna_pruned_pruned_air_fix_20261003.db",
     "pruned_air_fix_20261003"),
    ("optuna_pruned_pruned_walk_1_fix_20261003.db",
     "pruned_walk_1_fix_20261003"),
    ("optuna_pruned_pruned_walk_2_fix_20261003.db",
     "pruned_walk_2_fix_20261003"),
]

for db, name in TARGETS:
    try:
        st = optuna.load_study(study_name=name, storage="sqlite:///" + db)
        trials = st.get_trials(deepcopy=False)
        done = [t for t in trials if t.state.name == "COMPLETE"]
        running = [t for t in trials if t.state.name == "RUNNING"]
        dirs = [d.name for d in st.directions]
        key = max if "MAXIMIZE" in dirs else min
        best_v = best_n = None
        if done:
            b = key(done, key=lambda t: t.value)
            best_v, best_n = b.value, b.number
        print(f"{db:26s} {name:36s} dir={dirs[0][:3]} n={len(trials):3d} "
              f"complete={len(done):3d} run={len(running)} best={best_v} (t{best_n})")
    except Exception as e:
        print(f"{db:26s} {name:36s} ERROR {e}")
