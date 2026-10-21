"""Stage-3 balance campaign db stats (2026-10-21 resume session).
Trial counts, best trials, fall histograms for every balance study."""
import io
import json
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

TARGETS = [
    ("sqlite:///optuna_walk.db", "curr_s3_balance"),
    ("sqlite:///optuna_walk.db", "curr_s3p_perturb"),
    ("sqlite:///optuna_w2lvar.db", "curr_w2lvar_s3_balance"),
    ("sqlite:///optuna_w2lvar.db", "curr_w2lvar_s3p_perturb"),
    ("sqlite:///optuna_syn6.db", "curr_syn6_s3_balance"),
    ("sqlite:///optuna_syn6.db", "curr_syn6_s3p_perturb"),
]

out = {}
for db, name in TARGETS:
    try:
        st = optuna.load_study(study_name=name, storage=db)
    except Exception as e:
        out[f"{db}::{name}"] = f"LOAD FAIL: {e}"
        print(f"{db}::{name} LOAD FAIL: {e}")
        continue
    trials = st.trials
    n_done = [t for t in trials if t.state.name == "COMPLETE"]
    n_run = [t for t in trials if t.state.name == "RUNNING"]
    best = max(n_done, key=lambda t: t.value) if n_done else None
    hist = {}
    for t in n_done:
        fell = t.user_attrs.get("falls")
        if fell is None:
            v = t.value
            k = "fell" if v <= -149.5 else "stand"
        else:
            k = f"falls={fell}/7" if fell else "0falls"
        hist[k] = hist.get(k, 0) + 1
    rec = dict(n_total=len(trials), n_complete=len(n_done),
               n_running=len(n_run),
               best=(dict(number=best.number, value=best.value,
                          params=best.params,
                          user_attrs=best.user_attrs) if best else None),
               hist=hist)
    out[f"{db}::{name}"] = rec
    print(f"\n== {db}::{name}")
    print(f"   total={rec['n_total']} complete={rec['n_complete']} "
          f"running={rec['n_running']}")
    print(f"   hist={hist}")
    if best:
        print(f"   best t{best.number} = {best.value:.3f}")
        print(f"        {json.dumps(best.params)}")
        print(f"        attrs={json.dumps(best.user_attrs)}")

with open("_s3_stats_20261021.json", "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1, default=str)
print("\nwrote _s3_stats_20261021.json")
