"""Goal-5 reference-tension pilot v2 (2026-09-23) - the PRODUCTION config.

The v10-era --best10 vehicle no longer walks under current physics (pilot
v1: 100% double-support, no cycles at any drive - recorded as a finding).
This pilot runs the CURRENT production winner (curriculum_stage3.json =
s3k trial-34 seed) through the exact stage-3 machinery (_curriculum
.set_stage + runner --eval) with kine_ref.REF_CACHE swapped per reference:

  configs: seed drive x {0.85, 1.0, 1.15}
  refs:    subject01 (of record) | Falisse Case_40 (predicted gait)

Validity gate: drive x1.0 vs subject01 must reproduce ~ -159.45 (the
recorded s3k score). Read-only: no optuna db writes, no seed edits, no
core-file changes; npz outputs go to %TEMP%.
"""
import io
import json
import os
import sys
import time
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")
HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402
import kine_ref  # noqa: E402
from gait_lib_loader import load_reference_general  # noqa: E402

FAL = Path(r"D:\temp\gait_lib_staging\falisse\predictsim_mtp-master"
           r"\Results\Case_40")
REFS = {
    "subject01": kine_ref.load_reference(),
    "falisse": load_reference_general(FAL / "motion.mot", FAL / "GRF.mot"),
}

import _curriculum as C  # noqa: E402  (imports optuna, optuna_walk_v10, runner)
import runner as R  # noqa: E402
import params as P  # noqa: E402

# The s3k study's suggest space (the 32 db keys) never sampled the keys
# OW.set_params additionally requires; the s3k-era code must have carried
# defaults for them. Reconstruct NEUTRALLY: current TAU/BAL defaults and
# multiplicative factors of 1.0. The reproduction gate (drive x1.0 vs
# subject01 ~ -159.45) judges whether this reconstruction is faithful.
DEF_RG_ADAPT = P.TAU["rg_adapt"]
DEF_KX = P.BAL["kx"]

seed = json.loads((HERE / "curriculum_stage3.json").read_text())["params"]
seed.setdefault("rg_adapt", DEF_RG_ADAPT)
seed.setdefault("kx", DEF_KX)
for _k in ("e2_pf", "f1_df", "f1_kf", "e2_adapt", "post_kneext",
           "post_hipext"):
    seed.setdefault(_k, 1.0)
print(f"seed: s3k trial 34, recorded score -159.45, drive {seed['drive']:.3f}, "
      f"neutral rg_adapt={DEF_RG_ADAPT:.3f} kx={DEF_KX:.1f}")

tmp = Path(os.environ.get("TEMP", "."))
rows = []
for scale in (0.85, 1.0, 1.15):
    p = dict(seed)
    p["drive"] = seed["drive"] * scale
    for name, ref in REFS.items():
        kine_ref.REF_CACHE = ref
        C.set_stage(3, p)
        os.environ["AARL_NPZ"] = str(
            tmp / f"g5_p2_{name}_{p['drive']:.3f}.npz")
        t0 = time.perf_counter()
        m = R.main(["--eval", "--drive", repr(p["drive"])])
        wall = time.perf_counter() - t0
        row = dict(ref=name, drive=round(p["drive"], 3),
                   kine=m.get("kine_score"), duty_r=m.get("duty_r"),
                   duty_l=m.get("duty_l"), T_r=m.get("T_r"),
                   knee_min_r=m.get("knee_min_r"), wall=round(wall, 1))
        rows.append(row)
        print("RESULT " + ", ".join(f"{k}={v}" for k, v in row.items()),
              flush=True)

keys = list(rows[0])
out_csv = HERE / "reports_20260923" / "goal5_pilot2_results.csv"
out_csv.parent.mkdir(exist_ok=True)
with out_csv.open("w", encoding="utf-8") as f:
    f.write(",".join(keys) + "\n")
    for r in rows:
        f.write(",".join("" if r[k] is None else str(r[k]) for k in keys)
                + "\n")

by = {(r["ref"], r["drive"]): r["kine"] for r in rows}
lines = ["# Goal-5 pilot v2: production s3k config vs two references",
         "", "Seed = curriculum_stage3.json (s3k trial 34, recorded "
         "-159.45 vs subject01). kine lower = better.", "",
         "| drive | vs subject01 | vs Falisse |", "|---|---|---|"]
for scale in (0.85, 1.0, 1.15):
    d = round(seed["drive"] * scale, 3)
    lines.append(f"| {d} | {by.get(('subject01', d))} "
                 f"| {by.get(('falisse', d))} |")
rep = by.get(("subject01", round(seed["drive"], 3)))
lines += ["",
          f"Reproduction gate: drive x1.0 vs subject01 scored {rep} "
          f"(recorded -159.45; small drift expected - the curriculum "
          f"study's own runs pass through the same path).",
          "", "Read the drive row where each reference's score is best: "
          "if Falisse prefers a higher drive (its cycle is faster, "
          "T 1.113 s vs 1.233 s), the two references pull the cadence "
          "lever in different directions - the reference-tension this "
          "pilot measures."]
out_md = HERE / "reports_20260923" / "goal5_pilot2.md"
out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")
print("saved", out_csv, "and", out_md)
