import json
p = r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\connectome_templates.json"
t = json.load(open(p, encoding="utf-8"))
print("template keys:", [k for k in t if "bilateralrg" in k])
b = t.get("bilateralrg2")
if b:
    print("bilateralrg2 census: neurons", len(b.get("neurons", [])),
          "synapses", len(b.get("synapses", [])),
          "inputs", len(b.get("inputs", [])),
          "outputs", len(b.get("outputs", [])))

import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)
db = "sqlite:///D:/GitHub/Bipedal_Robot/Code/MuJoCo_SNS/spinal/reports_20261020/w2l_contact2/optuna_w2l_contact2.db".replace("\\", "/")
studies = optuna.get_all_study_names(storage=db)
print("studies:", studies)
for s in studies:
    st = optuna.load_study(study_name=s, storage=db)
    bt = st.best_trial
    print(f"  {s}: trials={len(st.trials)} best_obj={bt.value:.3f}")
