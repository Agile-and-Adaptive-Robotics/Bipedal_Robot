"""Spiking study status."""
import sys

import optuna

for name in ("curr_spk_s1_air_deaff_20261003",
             "curr_spk_s2_air_aff_20261003"):
    try:
        st = optuna.load_study(study_name=name,
                               storage="sqlite:///optuna_spiking.db")
    except KeyError:
        print(f"== {name}: MISSING")
        continue
    print(f"== {name}")
    for t in st.trials:
        d = ((t.datetime_complete - t.datetime_start).total_seconds() / 60
             if t.datetime_complete else None)
        print(t.number, t.state.name, t.value,
              t.datetime_start, f"{d:.1f}min" if d else "running")
