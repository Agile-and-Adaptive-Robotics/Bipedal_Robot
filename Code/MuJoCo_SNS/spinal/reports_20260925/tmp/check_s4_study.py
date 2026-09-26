"""Pre-launch check: state of curr_w2lvar_s4_walk_nocontact."""
import os

import optuna

os.chdir(r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal")
st = optuna.load_study(study_name="curr_w2lvar_s4_walk_nocontact",
                       storage="sqlite:///optuna_w2lvar.db")
print("n_trials:", len(st.trials))
for t in st.trials:
    print(t.number, t.state.name, t.value, t.params)
