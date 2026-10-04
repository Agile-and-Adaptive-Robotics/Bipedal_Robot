"""Pre-flight for the 2026-10-03 tuning-seed round. Read-only.

Phase a: env imports, fresh-name collision check in BOTH dbs, s3k
seed coverage (imports _curriculum).
Phase b: syn6 seed coverage (imports _curriculum_syn6).

Two phases because each curriculum module re-wraps sys.stdout at import
and a second wrap closes the first's shared buffer (python-skill
double-wrap trap); each phase imports exactly one of them.
"""
import json
import sys

SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)

PHASE = sys.argv[1] if len(sys.argv) > 1 else "a"
ok = True

if PHASE == "a":
    # import the stdout-wrapping module FIRST: anything printed before
    # its wrap is lost when the piped (block-buffered) stdout object is
    # replaced (observed this session)
    import _curriculum as C  # noqa: E402  (wraps stdout on import)

    import optuna  # noqa: E402

    print("optuna", optuna.__version__)
    import mujoco  # noqa: E402

    print("mujoco", mujoco.__version__)

    NEW = {"curr_s4_nocross_mutfix_20261003": "optuna_walk.db",
           "curr_syn6_s5_reseed_20261003": "optuna_syn6.db"}
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    for name, db in NEW.items():
        names = optuna.get_all_study_names(
            storage=f"sqlite:///{SPINAL}\\{db}")
        free = name not in names
        print(f"[names] {db}: '{name}' free = {free}")
        ok &= free

    s3k = json.loads(open(f"{SPINAL}\\reports_20260923\\"
                          f"s3k_trial34_full_params.json",
                          encoding="utf-8").read())["params"]
    miss = [k for k in C.STAGE_KEYS[4] if k not in s3k]
    print(f"[seed] s3k t34 params cover stage-4 keys: "
          f"{'PASS' if not miss else 'MISSING ' + str(miss)}")
    ok &= not miss
elif PHASE == "b":
    import _curriculum_syn6 as S  # noqa: E402  (wraps stdout on import)

    doc = json.loads(open(f"{SPINAL}\\curriculum_syn6_stage5.json",
                          encoding="utf-8").read())
    miss = [k for k in S.STAGE_KEYS[5] if k not in doc["params"]]
    print(f"[seed] syn6 stage-5 winner (score {doc['score']}, study "
          f"{doc.get('study')}) covers stage-5 keys: "
          f"{'PASS' if not miss else 'MISSING ' + str(miss)}")
    ok &= not miss
else:
    print("unknown phase", PHASE)
    sys.exit(2)

print("PREFLIGHT", PHASE, "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
