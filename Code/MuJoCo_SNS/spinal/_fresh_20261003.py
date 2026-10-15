"""FRESH-STUDY driver, 2026-10-03 post-mutual-inhibition-fix campaign.

The 2026-10-03 wiring audit (WIRING_RULINGS_20261003.md, change 1) added
the missing IaIN<->IaIN / IBIN<->IBIN mutual-inhibition edges to the
s3k/w2lvar builders, so every pre-fix study score for those builders is
stale. This driver launches FRESH DATED studies through the canonical
curriculum modules WITHOUT modifying them (study-name redirect via the
modules' own name tables / helper calls).

Usage:
    python _fresh_20261003.py w2lvar <stage 1..5> <n_trials>
    python _fresh_20261003.py pruned air   <n_trials>
    python _fresh_20261003.py pruned walk  <idx 1..3> <n_trials>

w2lvar: study curr_w2lvar_s{S}_fix_20261003 in optuna_w2lvar.db; the
canonical curriculum_w2lvar_stage{S}.json is backed up once to
curriculum_w2lvar_stage{S}_prefix_20261003.json (pre-fix record) and then
rewritten by the module's own main() so the fresh chain seeds itself.
pruned: study pruned_air_fix_20261003 / pruned_walk_{I}_fix_20261003 in
their own per-study db files (the module's contention fix).
"""
import json
import os
import shutil
import sys

TAG = "fix_20261003"

model = sys.argv[1]

if model == "w2lvar":
    stage = int(sys.argv[2])
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    import _curriculum_w2lvar as CW

    dated = f"curr_w2lvar_s{stage}_{TAG}"
    CW.STUDY_NAMES[stage] = dated
    CW.NPZ = f"spinal_run_{dated}.npz"
    canon = f"curriculum_w2lvar_stage{stage}.json"
    backup = f"curriculum_w2lvar_stage{stage}_prefix_20261003.json"
    if os.path.exists(canon) and not os.path.exists(backup):
        shutil.copyfile(canon, backup)
        print(f"[fresh] backed up {canon} -> {backup}")
    if stage >= 2:
        # the module's seed chain reads curriculum_w2lvar_stage{N}.json
        # FIRST - that is the OLD pre-fix winner. Move it aside so the
        # FileNotFoundError fallback picks stage N-1 = the FRESH winner
        # written by this chain's previous stage run.
        prev_fresh = f"curriculum_w2lvar_stage{stage-1}.json"
        if os.path.exists(canon) and os.path.exists(prev_fresh):
            os.replace(canon, backup)
            print(f"[fresh] moved stale {canon} aside "
                  f"(seed chains from the fresh stage "
                  f"{stage-1} winner)")
    sys.argv = ["_curriculum_w2lvar.py", str(stage), str(n)]
    CW.main()
    shutil.copyfile(canon, f"curriculum_w2lvar_stage{stage}_{TAG}.json")
    print(f"[fresh] copied winner to curriculum_w2lvar_stage{stage}_{TAG}.json")

elif model == "pruned":
    cmd = sys.argv[2]
    if cmd == "air":
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 30
        import _curriculum_pruned as CP

        name = f"pruned_air_{TAG}"
        st = CP._study(name, 11, CP._air_obj, CP.AIR_SEED)
        st.optimize(CP._air_obj, n_trials=n, gc_after_trial=True)
        b = st.best_trial
        print(f"== pruned air best {b.value:.3f} (trial {b.number})")
        json.dump({"stage": "air", "score": b.value, "params": b.params,
                   "trial": b.number, "study": name,
                   "db": f"optuna_pruned_{name}.db", "pin": CP.PIN,
                   "prune_keys": CP.PRUNE_KEYS, "prune_argv": CP.PRUNE_ARGV},
                  open(f"curriculum_pruned_air_{TAG}.json", "w",
                       encoding="utf-8"), indent=2)
    elif cmd == "walk":
        idx = int(sys.argv[3])
        n = int(sys.argv[4]) if len(sys.argv) > 4 else 40
        import _curriculum_pruned as CP

        # seed the walk study from the fresh post-fix AIR winner
        seed = dict(CP.WALK_SEED)
        try:
            airw = json.loads(open(f"curriculum_pruned_air_{TAG}.json",
                                   encoding="utf-8").read())["params"]
            for k in seed:
                if k in airw:
                    seed[k] = float(airw[k])
        except FileNotFoundError:
            pass
        name = f"pruned_walk_{idx}_{TAG}"
        st = CP._study(name, 40 + idx, CP._walk_obj, seed)
        st.optimize(CP._walk_obj, n_trials=n, gc_after_trial=True)
        b = st.best_trial
        print(f"== {name} best {b.value:.3f} (trial {b.number})")
        json.dump({"score": b.value, "params": b.params, "trial": b.number,
                   "study": name, "db": f"optuna_pruned_{name}.db",
                   "pin": CP.PIN, "prune_keys": CP.PRUNE_KEYS,
                   "prune_argv": CP.PRUNE_ARGV},
                  open(f"curriculum_pruned_walk_{idx}_{TAG}.json", "w",
                       encoding="utf-8"), indent=2)
    else:
        raise SystemExit(f"unknown pruned cmd {cmd}")

else:
    raise SystemExit(f"unknown model {model}")
