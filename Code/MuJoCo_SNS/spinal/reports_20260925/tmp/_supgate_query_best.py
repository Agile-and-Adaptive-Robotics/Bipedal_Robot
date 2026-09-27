"""SUPERVISOR GATE: independent winner audit of the 10 variant studies
(curr_w2lvar_s1..s5 in optuna_w2lvar.db, curr_syn6_s1..s5 in optuna_syn6.db).
Read-only via optuna: db argmax over COMPLETE trials vs recorded winner in
curriculum_*.json (score + trial + params). Prints per-study trial counts,
best value, repeated-value census (plateau/sentinel visibility).
"""
import io
import json
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
os.chdir(r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal")

import optuna

optuna.logging.set_verbosity(optuna.logging.WARNING)

VARIANTS = {
    "w2lvar": "optuna_w2lvar.db",
    "syn6": "optuna_syn6.db",
}

overall_ok = True
for var, db in VARIANTS.items():
    storage = f"sqlite:///{db}"
    names = sorted(s.study_name for s in
                   optuna.get_all_study_summaries(storage=storage))
    print(f"=== {var}: {db}")
    print(f"    studies: {names}")
    for i in range(1, 6):
        jf = f"curriculum_{var}_stage{i}.json"
        pref = f"curr_{var}_s{i}"
        cand = [n for n in names if n.startswith(pref)]
        if len(cand) != 1:
            print(f"--- stage {i}: STUDY NAME RESOLUTION FAILED: {cand}")
            overall_ok = False
            continue
        st = optuna.load_study(study_name=cand[0], storage=storage)
        done = [t for t in st.trials
                if t.state == optuna.trial.TrialState.COMPLETE]
        if not done:
            print(f"--- stage {i}: NO COMPLETE TRIALS")
            overall_ok = False
            continue
        best = max(done, key=lambda t: t.value)
        j = json.loads(open(jf, encoding="utf-8").read())
        tops = sorted(done, key=lambda t: -t.value)[:3]
        floors = {}
        for t in done:
            key = round(t.value, 3)
            floors[key] = floors.get(key, 0) + 1
        hot = {k: c for k, c in floors.items() if c >= 3}
        t_eq = best.number == j.get("trial")
        v_delta = abs(best.value - float(j["score"]))
        p_eq = None
        try:
            p_eq = all(abs(float(best.params[kk]) - float(j["params"][kk]))
                       < 1e-12 for kk in j["params"]) and \
                set(j["params"]) <= set(best.params)
        except Exception as e:  # noqa: BLE001
            p_eq = f"ERR {e!r}"
        ok = t_eq and v_delta < 1e-6 and p_eq is True
        overall_ok = overall_ok and ok
        print(f"--- stage {i}: study={cand[0]} n_trials={len(st.trials)} "
              f"n_complete={len(done)}")
        print(f"    db best : trial {best.number} value {best.value!r} "
              f"(last completion {best.datetime_complete})")
        print(f"    json win: trial {j.get('trial')} score {j.get('score')!r}")
        print(f"    params match: {p_eq}")
        print(f"    MATCH trial={t_eq} value_delta={v_delta:.3e} -> "
              f"{'OK' if ok else 'FAIL'}")
        print(f"    top3: {[(t.number, round(t.value, 6)) for t in tops]}")
        print(f"    values shared by >=3 trials: "
              f"{dict(sorted(hot.items(), key=lambda kv: -kv[1])[:6])}")
print()
print("OVERALL:", "PASS" if overall_ok else "FAIL", flush=True)
