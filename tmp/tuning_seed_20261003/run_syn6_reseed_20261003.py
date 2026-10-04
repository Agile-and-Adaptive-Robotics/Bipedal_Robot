"""Tuning-seed round 2026-10-03: syn6 stage 5 (walk with contact) under
a FRESH DATED study name, seeded fresh from the recorded stage-5 winner.

Audit condition satisfied: WIRING_RULINGS_20261003.md verdict for syn6
is PASS on the synergy implementation and the model was NOT affected by
the mutual-inhibition fix (build unchanged at 794/382/3033), so the
ask's condition for reseeding syn6 ("only if the wiring audit did not
flag its synergy layer as broken") holds. The S1 heel-side conflict is
a documented ruling with a retune recommendation that requires a network
change; per this round's constraints no network code is modified, so the
rescan runs on the current wiring.

This launcher RUNS THE EXISTING STUDY INFRASTRUCTURE: objective,
set_stage, search ranges, sentinels, env pins, and seed construction are
imported verbatim from _curriculum_syn6 (nothing modified anywhere). The
ONLY deltas, each deliberate and documented:
  1. study name  curr_syn6_s5_reseed_20261003  (dated, fresh),
  2. output json  curriculum_syn6_stage5_reseed_20261003.json  so the
     recorded canonical curriculum_syn6_stage5.json (curr_syn6_s5 winner)
     is NOT overwritten,
  3. optimize(timeout=5100 s) so no new trial is LAUNCHED after 85 min
     (the ask's hard 90-min budget with margin for the in-flight trial).
AARL_NET / AARL_NPZ are pinned exactly as _curriculum_syn6.main() pins
them. No number in any searched range, seed, or gain is changed.
"""
import json
import os
import sys

SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
os.chdir(SPINAL)
sys.path.insert(0, SPINAL)

STUDY = "curr_syn6_s5_reseed_20261003"
OUT_JSON = "curriculum_syn6_stage5_reseed_20261003.json"
STAGE = 5
N_TRIALS = 16
TIMEOUT_S = 5100

# NOTE: do NOT wrap sys.stdout here; _curriculum_syn6 wraps it on import
# and a second wrap would close the first (python-skill double-wrap trap).
import _curriculum_syn6 as S  # noqa: E402  (import after path setup)

# ---- verbatim env pins from _curriculum_syn6.main() (lines 250-255) ----
os.environ["AARL_NET"] = S.VARIANT
os.environ["AARL_NPZ"] = S.NPZ
print(f"[syn6-reseed] pinned AARL_NET={os.environ['AARL_NET']} "
      f"AARL_NPZ={os.environ['AARL_NPZ']} DB={S.DB} study={STUDY} "
      f"n_trials={N_TRIALS} timeout_s={TIMEOUT_S}", flush=True)

# ---- verbatim from _curriculum_syn6.main() (lines 258-263) ----
prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
S.BASE_MUL = dict(prev["multipliers"])
S.BASE_MUL["renshaw"] = 0.0
S.BASE_MUL["syn6"] = 1.0
S.BASE_MUL["syn6_brainstem"] = 0.0

import optuna  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)
study = optuna.create_study(direction="maximize", storage=S.DB,
                            study_name=STUDY, load_if_exists=True,
                            sampler=optuna.samplers.TPESampler(
                                seed=21 + STAGE, n_startup_trials=8))
if len(study.trials) == 0:
    # ---- verbatim seed construction from _curriculum_syn6.main()
    # (lines 273-296) ----
    sk = S.STAGE_KEYS[STAGE]
    seed = dict(S.SEED_DEFAULTS[STAGE])
    seed["syn6"] = 1.0
    try:
        prevw = json.loads(open(f"curriculum_{S.VARIANT}_"
                                f"stage{STAGE}.json",
                                encoding="utf-8").read())["params"]
    except FileNotFoundError:
        try:
            prevw = json.loads(open(
                f"curriculum_{S.VARIANT}_stage{STAGE-1}.json",
                encoding="utf-8").read())["params"]
        except FileNotFoundError:
            prevw = {}
    for k in sk:
        if k in prevw:
            seed[k] = float(prevw[k])
    study.enqueue_trial(seed)
    print("seeded", json.dumps(seed), flush=True)

study.optimize(S.objective(STAGE), n_trials=N_TRIALS,
               gc_after_trial=True, timeout=TIMEOUT_S)
best = study.best_trial
print(f"== stage {STAGE} [{STUDY}] best {best.value:.3f} "
      f"(trial {best.number})", flush=True)
print(json.dumps(best.params, indent=1), flush=True)
with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump({"stage": STAGE, "score": best.value,
               "params": best.params, "trial": best.number,
               "study": STUDY, "db": S.DB, "npz": S.NPZ,
               "aarl_net": S.VARIANT, "g_syn6": 1.0,
               "g_syn6_brainstem": 0.0,
               "n_trials_in_study": len(study.trials),
               "provenance": "fresh-seed rescan of syn6 stage 5 after the "
                             "2026-10-03 wiring audit (synergy PASS; "
                             "build unchanged); seeded from the recorded "
                             "curr_syn6_s5_walk_contact winner"}, f,
              indent=2)
print(f"saved {OUT_JSON}", flush=True)
