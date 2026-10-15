"""Fail stale RUNNING trials in the campaign studies (killed processes).
Optional argv[1]: substring filter on study name (matched against study
name OR db filename). NEVER run this while a live process may be
mid-trial on any listed db (2026-10-03 lesson: an externally failed
live trial crashes study.optimize on its return).
"""
import sys

import optuna

TARGETS = [
    ("optuna_walk.db", "curr_s4_nocross_mutfix_20261003"),
    ("optuna_syn6.db", "curr_syn6_s5_reseed_20261003"),
    ("optuna_spiking.db", "curr_spk_s1_air_deaff_20261003"),
    ("optuna_spiking.db", "curr_spk_s2_air_aff_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s1_fix_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s2_fix_20261003"),
    ("optuna_w2lvar.db", "curr_w2lvar_s3_fix_20261003"),
    ("optuna_pruned_pruned_air_fix_20261003.db",
     "pruned_air_fix_20261003"),
    ("optuna_pruned_pruned_walk_1_fix_20261003.db",
     "pruned_walk_1_fix_20261003"),
]
flt = sys.argv[1] if len(sys.argv) > 1 else ""
for db, name in TARGETS:
    if flt and flt not in name and flt not in db:
        continue
    try:
        st = optuna.load_study(study_name=name, storage="sqlite:///" + db)
    except KeyError:
        continue
    for t in st.trials:
        if t.state.name == "RUNNING":
            try:
                st.tell(t.number, state=optuna.trial.TrialState.FAIL)
                print(f"{name}: failed stale trial {t.number}")
            except Exception as e:
                print(f"{name}: trial {t.number}: {e}")
    done = [t for t in st.trials if t.state.name == "COMPLETE"]
    if done:
        best = max(done, key=lambda t: t.value)
        print(f"{name}: {len(done)} complete, best "
              f"{best.value:.4f} (t{best.number})")
    else:
        print(f"{name}: none complete")
