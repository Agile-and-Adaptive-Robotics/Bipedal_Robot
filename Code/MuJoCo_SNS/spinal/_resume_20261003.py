"""RESUME DRIVER for the 2026-10-03 overnight curriculum campaign.

Resumes the two post-mutual-inhibition-fix studies WITHOUT modifying the
canonical curriculum scripts: the study-name redirect happens via
monkey-patch (s3k: optuna.create_study wrapper; syn6: module-level
STUDY_NAMES override), then the module's own main() runs unmodified.

Usage:
    python _resume_20261003.py s3k   <n_trials>   # stage 4 mutfix study
    python _resume_20261003.py syn6  <n_trials>   # stage 5 reseed study

Each process pins a unique AARL_NPZ so concurrent runs never clobber
spinal_run.npz or each other. The canonical stage json the script writes
is additionally copied to the dated campaign json name.
"""
import json
import os
import shutil
import sys

MODEL = sys.argv[1]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 30

if MODEL == "s3k":
    STUDY = "curr_s4_nocross_mutfix_20261003"
    os.environ["AARL_NPZ"] = f"spinal_run_{STUDY}_resume.npz"
    import optuna

    _orig_create = optuna.create_study

    def _create_override(**kw):
        kw["study_name"] = STUDY
        return _orig_create(**kw)

    optuna.create_study = _create_override
    sys.argv = ["_curriculum.py", "4", str(N)]
    import _curriculum as C
    C.main()
    shutil.copyfile("curriculum_stage4.json",
                    "curriculum_stage4_mutfix_20261003.json")
    print(f"[resume] refreshed curriculum_stage4_mutfix_20261003.json")

elif MODEL == "syn6":
    STUDY = "curr_syn6_s5_reseed_20261003"
    import _curriculum_syn6 as CS

    CS.STUDY_NAMES[5] = STUDY
    # main() pins AARL_NET/AARL_NPZ itself; only the npz name needs a
    # unique suffix for this resume process
    CS.NPZ = f"spinal_run_{STUDY}_resume.npz"
    sys.argv = ["_curriculum_syn6.py", "5", str(N)]
    CS.main()
    shutil.copyfile("curriculum_syn6_stage5.json",
                    "curriculum_syn6_stage5_reseed_20261003.json")
    print(f"[resume] refreshed curriculum_syn6_stage5_reseed_20261003.json")

else:
    raise SystemExit(f"unknown model {MODEL}")
