"""Pre-launch check: state of curr_w2lvar_s1_air_deaff in optuna_w2lvar.db."""
import os
from collections import Counter

import optuna

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
st = optuna.load_study(study_name="curr_w2lvar_s1_air_deaff",
                       storage="sqlite:///optuna_w2lvar.db")
print("n_trials:", len(st.trials))
print(Counter(t.state.name for t in st.trials))
for t in st.trials:
    print(t.number, t.state.name, t.value, t.params)
print("other studies in db:",
      sorted(s.study_name for s in optuna.study.get_all_studies(
          storage="sqlite:///optuna_w2lvar.db")))
