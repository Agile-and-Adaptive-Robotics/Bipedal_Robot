"""Print the number of trials a curriculum stage still needs (target 100)."""
import sys

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)
name = {1: "curr_s1_air_deaff", 2: "curr_s2_air_aff",
        3: "curr_s3_ground"}[int(sys.argv[1])]
try:
    st = optuna.load_study(study_name=name,
                           storage="sqlite:///optuna_walk.db")
    print(max(0, 100 - len(st.trials)))
except KeyError:
    print(100)
