"""Read-only peek at the curriculum optuna studies: per-trial values +
params, to see whether any stage-1/2 trial showed real rhythm (the
objective is 3*rises + 0.5*(-knee_min); rises >= 4 needs the E side to
actually burst)."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
st = optuna.load_study(study_name="curr_s1_air_deaff",
                       storage="sqlite:///optuna_walk.db")
print(f"=== curr_s1_air_deaff: {len(st.trials)} trials ===")
for tr in st.trials:
    if tr.value is None:
        continue
    # infer rises: 3r + 0.5k = v with k in [-90, 0] -> r in [0..30]
    print(f"  t{tr.number:3d} v={tr.value:8.3f} "
          f"drive={tr.params.get('drive', 0):.2f}")
print("best:", st.best_value, "trial", st.best_trial.number,
      st.best_trial.params)

for name in ("curr_s2_air_aff", "curr_s3_ground"):
    try:
        s2 = optuna.load_study(study_name=name,
                               storage="sqlite:///optuna_walk.db")
    except KeyError:
        print(f"=== {name}: MISSING ===")
        continue
    vals = [t.value for t in s2.trials if t.value is not None]
    print(f"=== {name}: {len(s2.trials)} trials, "
          f"best {max(vals) if vals else float('nan'):.3f} ===")
    top = sorted((t for t in s2.trials if t.value is not None),
                 key=lambda t: -t.value)[:5]
    for tr in top:
        # knee_min implied if rises==0: k = 2*v
        print(f"  t{tr.number:3d} v={tr.value:8.3f} "
              f"drive={tr.params.get('drive', 0):.2f} "
              f"heel={tr.params.get('heel_rge', -1):.2f}")
