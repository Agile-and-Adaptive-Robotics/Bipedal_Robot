"""Refresh the campaign winner jsons from the dbs (the killed resume
processes never reached main()'s json write)."""
import json

import optuna

JOBS = [
    ("optuna_walk.db", "curr_s4_nocross_mutfix_20261003",
     ["curriculum_stage4.json", "curriculum_stage4_mutfix_20261003.json"],
     4, None),
    ("optuna_syn6.db", "curr_syn6_s5_reseed_20261003",
     ["curriculum_syn6_stage5.json",
      "curriculum_syn6_stage5_reseed_20261003.json"], 5, "syn6"),
]
for db, name, outs, stage, variant in JOBS:
    st = optuna.load_study(study_name=name, storage="sqlite:///" + db)
    done = [t for t in st.trials if t.state.name == "COMPLETE" and
            t.value is not None]
    b = max(done, key=lambda t: t.value)
    doc = {"stage": stage, "score": b.value, "params": b.params,
           "trial": b.number, "study": name, "db": "sqlite:///" + db}
    if variant:
        doc.update({"npz": f"spinal_run_{variant}.npz",
                    "aarl_net": variant})
    for out in outs:
        json.dump(doc, open(out, "w", encoding="utf-8"), indent=2)
    print(f"{name}: best {b.value:.4f} (t{b.number}) -> "
          + ", ".join(outs))
