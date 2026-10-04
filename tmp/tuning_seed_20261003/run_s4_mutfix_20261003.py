"""Tuning-seed round 2026-10-03: stock s3k-lineage stage 4 (ground walk,
KEYS_WALK) under a FRESH DATED study name.

Why this run: the 2026-10-03 wiring audit completed the one-directional
IaIN/IBIN mutual inhibition under full_rules (WIRING_RULINGS_20261003.md
change 1), so every recorded full_rules>0 score is a pre-fix value; this
study re-evals the s3k trial-34 incumbent (enqueued as trial 0) and
retunes from it under the corrected topology.

This launcher RUNS THE EXISTING STUDY INFRASTRUCTURE: objective,
set_stage, search ranges, sentinels, and seed construction are imported
verbatim from _curriculum (nothing in that module or in any network code
is modified). The ONLY deltas, each deliberate and documented:
  1. study name  curr_s4_nocross_mutfix_20261003  (dated; the canonical
     curr_s4_nocross study has never been created, but a dated name also
     carries the post-fix provenance the audit asks for),
  2. output json  curriculum_stage4_mutfix_20261003.json  (the canonical
     curriculum_stage4.json stays absent),
  3. AARL_NPZ redirect so this run cannot clobber spinal_run.npz,
  4. optimize(timeout=5100 s) so no new trial is LAUNCHED after 85 min
     (the ask's hard 90-min budget with margin for the in-flight trial).
No number in any searched range, seed, or gain is changed.
"""
import json
import os
import sys

SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
os.chdir(SPINAL)
sys.path.insert(0, SPINAL)

STUDY = "curr_s4_nocross_mutfix_20261003"
OUT_JSON = "curriculum_stage4_mutfix_20261003.json"
NPZ = "spinal_run_s4_mutfix_20261003.npz"
STAGE = 4
N_TRIALS = 18
TIMEOUT_S = 5100

os.environ["AARL_NPZ"] = NPZ

# NOTE: do NOT wrap sys.stdout here; _curriculum wraps it on import and a
# second wrap would close the first (python-skill double-wrap trap).
import _curriculum as C  # noqa: E402  (import after path setup)
import optuna  # noqa: E402

print(f"[s4-mutfix] study={STUDY} db={C.DB} npz={NPZ} "
      f"n_trials={N_TRIALS} timeout_s={TIMEOUT_S}", flush=True)

# ---- verbatim from _curriculum.main() (lines 401-413) ----
prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
C.BASE_MUL = dict(prev["multipliers"])
C.BASE_MUL["renshaw"] = 0.5
s3k = json.loads(open(
    "reports_20260923/s3k_trial34_full_params.json",
    encoding="utf-8").read())["params"]
C.BASE_MUL.update({k: v for k, v in s3k.items()})
if s3k.get("pf_gain") is None:
    C.BASE_MUL.pop("pf_gain", None)

optuna.logging.set_verbosity(optuna.logging.WARNING)
study = optuna.create_study(direction="maximize",
                            storage=C.DB, study_name=STUDY,
                            load_if_exists=True,
                            sampler=optuna.samplers.TPESampler(
                                seed=21 + STAGE, n_startup_trials=8))
if len(study.trials) == 0:
    # ---- verbatim seed construction from _curriculum.main() ----
    sk = C.STAGE_KEYS[STAGE]
    seed = {k: 0.0 for k in sk}
    # stages 4-5 chain from the s3k PRODUCTION winner directly
    prevw = s3k
    for k in sk:
        if k in prevw:
            seed[k] = float(prevw[k])
    if "c1_gain" in sk:
        seed["c1_gain"] = max(seed["c1_gain"], 0.1)
    if "full_rules" in sk:
        seed["full_rules"] = 1.0
    if "no_cross" in sk:
        seed["no_cross"] = 1.0
    study.enqueue_trial(seed)
    print("seeded", json.dumps(seed), flush=True)

study.optimize(C.objective(STAGE), n_trials=N_TRIALS,
               gc_after_trial=True, timeout=TIMEOUT_S)
best = study.best_trial
print(f"== stage {STAGE} [{STUDY}] best {best.value:.3f} "
      f"(trial {best.number})", flush=True)
print(json.dumps(best.params, indent=1), flush=True)
with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump({"stage": STAGE, "score": best.value,
               "params": best.params, "trial": best.number,
               "study": STUDY, "db": C.DB, "npz": NPZ,
               "n_trials_in_study": len(study.trials),
               "provenance": "post-mutual-inhibition-fix continuation "
                              "of the s3k lineage (WIRING_RULINGS_20261003 "
                              "change 1); seeded from s3k trial 34"}, f,
              indent=2)
print(f"saved {OUT_JSON}", flush=True)
