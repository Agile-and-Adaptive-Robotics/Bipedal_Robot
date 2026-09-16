"""Generate optuna_walk_v10.py from optuna_walk_v9.py (transient PRESET
+ ankle trim study)."""
import io

t = io.open("optuna_walk_v9.py", encoding="utf-8").read()
t = t.replace("ground_walk_v9_flat", "ground_walk_v10_transient")
t = t.replace("best_walk_params_v9.json", "best_walk_params_v10.json")
t = t.replace("v9_best_trial", "v10_best_trial")
t = t.replace("v9_results.csv", "v10_results.csv")
t = t.replace("runner --best9 reproduces", "runner --best10 reproduces")
t = t.replace(
    'f1_anklepf_inh=trial.suggest_float("f1_anklepf_inh", 0.0, 2.0),',
    'f1_anklepf_inh=trial.suggest_float("f1_anklepf_inh", 0.0, 2.0),\n'
    '        ankle_post_walk_trim=trial.suggest_float('
    '"ankle_post_walk_trim", 0.05, 1.0),')
t = t.replace(
    'params.G["f1_anklepf_inh"] = float(p.get("f1_anklepf_inh", 0.0))',
    'params.G["f1_anklepf_inh"] = float(p.get("f1_anklepf_inh", 0.0))\n'
    '    params.G["ankle_post_walk_trim"] = '
    'float(p.get("ankle_post_walk_trim", 1.0))')
t = t.replace(
    'f1_anklepf_inh=best.params["f1_anklepf_inh"],\n'
    '               renshaw=0.5)',
    'f1_anklepf_inh=best.params["f1_anklepf_inh"],\n'
    '               ankle_post_walk_trim=best.params['
    '"ankle_post_walk_trim"],\n'
    '               renshaw=0.5)')
t = t.replace('(HERE / "best_walk_params_v8.json").read_text("utf-8"))',
              '(HERE / "best_walk_params_v9.json").read_text("utf-8"))')
t = t.replace(
    "seeded v8 hand-pose winner into v8b (normal.mot pose, 16 s eval)",
    "seeded v9 winner into v10 (transient PRESET + ankle trim)")
io.open("optuna_walk_v10.py", "w", encoding="utf-8").write(t)
n = t.count("ankle_post_walk_trim")
print(f"v10 written ({n} trim references)")
