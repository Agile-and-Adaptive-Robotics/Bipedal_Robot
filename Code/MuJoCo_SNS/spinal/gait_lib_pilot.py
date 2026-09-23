"""Goal-5 reference-tension pilot - ONE config per process (2026-09-23).

Usage: python gait_lib_pilot.py <ref_name> <drive>
  ref_name: subject01 | falisse
Scores the tuned walker (--fitted --best10) at the given --drive against
the chosen reference by swapping kine_ref.REF_CACHE. Fresh process per
config (runner.main's --best loaders mutate the params module; re-entry
would double-apply pf_gain). Appends one CSV row to
reports_20260923/goal5_pilot_results.csv. npz goes to %TEMP%.
No core edits, no optuna db writes, no seed-file changes.
"""
import io
import os
import sys
import time
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

REF_NAME = sys.argv[1] if len(sys.argv) > 1 else "subject01"
DRIVE = float(sys.argv[2]) if len(sys.argv) > 2 else 2.93

import kine_ref
if REF_NAME == "falisse":
    from gait_lib_loader import load_reference_general
    FAL = Path(r"D:\temp\gait_lib_staging\falisse\predictsim_mtp-master"
               r"\Results\Case_40")
    kine_ref.REF_CACHE = load_reference_general(FAL / "motion.mot",
                                                FAL / "GRF.mot")
    assert kine_ref.REF_CACHE is not None
elif REF_NAME == "subject01":
    kine_ref.REF_CACHE = kine_ref.load_reference()
else:
    raise SystemExit(f"unknown ref {REF_NAME}")

import runner  # after REF_CACHE is set (runner imports kine_ref lazily? be safe)

tmp = Path(os.environ.get("TEMP", "."))
os.environ["AARL_NPZ"] = str(tmp / f"g5_pilot_{REF_NAME}_{DRIVE:.2f}.npz")

t0 = time.perf_counter()
m = runner.main(["--eval", "--fitted", "--best10", "--drive", repr(DRIVE)])
wall = time.perf_counter() - t0

row = {
    "ref": REF_NAME, "drive": DRIVE,
    "kine": m.get("kine_score"), "duty_r": m.get("duty_r"),
    "duty_l": m.get("duty_l"), "T_r": m.get("T_r"), "T_l": m.get("T_l"),
    "knee_min_r": m.get("knee_min_r"), "wall": round(wall, 1),
}
keys = list(row)
out = HERE / "reports_20260923" / "goal5_pilot_results.csv"
out.parent.mkdir(exist_ok=True)
new = not out.exists()
with out.open("a", encoding="utf-8") as f:
    if new:
        f.write(",".join(keys) + "\n")
    f.write(",".join("" if row[k] is None else str(row[k]) for k in keys)
            + "\n")
print("RESULT " + ", ".join(f"{k}={row[k]}" for k in keys))
