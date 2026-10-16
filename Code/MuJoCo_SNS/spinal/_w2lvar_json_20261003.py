"""Write the FRESH w2lvar stage jsons from the fix studies (used when a
chain died before main()'s json write)."""
import json
import sys

import optuna

MAP = {1: ("optuna_w2lvar.db", "curr_w2lvar_s1_fix_20261003"),
       2: ("optuna_w2lvar.db", "curr_w2lvar_s2_fix_20261003"),
       3: ("optuna_w2lvar.db", "curr_w2lvar_s3_fix_20261003")}

stage = int(sys.argv[1]) if len(sys.argv) > 1 else 1
db, name = MAP[stage]
st = optuna.load_study(study_name=name, storage="sqlite:///" + db)
done = [t for t in st.trials if t.state.name == "COMPLETE" and
        t.value is not None]
if not done:
    raise SystemExit(f"{name}: no complete trials")
for t in st.trials:
    if t.state.name == "RUNNING":
        try:
            st.tell(t.number, state=optuna.trial.TrialState.FAIL)
            print(f"failed stale {t.number}")
        except Exception:
            pass
b = max(done, key=lambda t: t.value)
doc = {"stage": stage, "score": b.value, "params": b.params,
       "trial": b.number, "study": name, "db": "sqlite:///" + db,
       "npz": "spinal_run_w2lvar.npz", "aarl_net": "w2lvar"}
json.dump(doc, open(f"curriculum_w2lvar_stage{stage}.json", "w",
                    encoding="utf-8"), indent=2)
print(f"stage {stage}: best {b.value:.4f} (t{b.number}) -> "
      f"curriculum_w2lvar_stage{stage}.json")
