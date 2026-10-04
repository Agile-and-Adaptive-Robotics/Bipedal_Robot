"""Final read-only summary of the 2026-10-03 tuning-seed round.

Per study: trial counts by state, trial-0 (enqueued incumbent seed)
value, best trial value + number. Read-only.
"""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
RUNS = [("curr_s4_nocross_mutfix_20261003", "optuna_walk.db"),
        ("curr_syn6_s5_reseed_20261003", "optuna_syn6.db")]

for name, db in RUNS:
    st = optuna.load_study(study_name=name,
                           storage=f"sqlite:///{SPINAL}\\{db}")
    states = {}
    for t in st.trials:
        states[t.state.name] = states.get(t.state.name, 0) + 1
    comp = [t for t in st.trials if t.state.name == "COMPLETE"]
    best = max(comp, key=lambda t: t.value)
    t0 = next((t for t in st.trials if t.number == 0), None)
    print(f"== {name} ({db})")
    print(f"   trials={len(st.trials)} states={states}")
    if t0 is not None and t0.value is not None:
        print(f"   trial0 (seed incumbent) value={t0.value!r}")
    print(f"   best value={best.value!r} (trial {best.number})")
    print(f"   best params={best.params}")
