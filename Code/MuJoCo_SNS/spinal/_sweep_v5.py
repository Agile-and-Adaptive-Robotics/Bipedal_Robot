"""v5 phase-reset hand sweep: 4x4 grid over (phase_reset_e, phase_reset_f),
parallel batches of runner evals, results to v5_sweep.csv + stdout table.

Usage: python _sweep_v5.py [batch]
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import params  # noqa: F401  (baseline for reference values only)

HERE = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
PY = r"C:\Users\Ben Bolen\.conda\envs\myo\python.exe"
CSV = HERE + r"\v5_sweep.csv"
GRID = (0.0, 0.25, 0.5, 1.0)
FIELDS = ["phase_reset_e", "phase_reset_f", "kine_score", "duty", "duty_ref",
          "knee_min", "knee_min_ref", "rmse_hip", "rmse_knee", "rmse_ankle",
          "cadence", "tilt_max", "stayed_up", "nan", "n_cycles"]


def run_one(ge: float, gf: float) -> dict:
    cmd = [PY, "_eval_once.py", "--fitted", "--best", "--eval",
           "--phase-reset", repr(ge), repr(gf)]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=1800)
    line = next((ln for ln in p.stdout.splitlines()
                 if ln.startswith("METRICS ")), None)
    if line is None:
        print(f"({ge},{gf}) FAILED rc={p.returncode}: "
              f"{p.stdout[-300:]} | {p.stderr[-300:]}", flush=True)
        return dict(phase_reset_e=ge, phase_reset_f=gf, kine_score="FAIL")
    d = json.loads(line[8:])
    m, k = d["metrics"], d["kine"] or {}
    row = dict(phase_reset_e=ge, phase_reset_f=gf,
               kine_score=f"{m['kine_score']:.3f}",
               duty=round(k.get("duty", float("nan")), 3), duty_ref=0.61,
               knee_min=round(k.get("knee_min", float("nan")), 2),
               knee_min_ref=-69.7,
               rmse_hip=round(k.get("rmse_hip", float("nan")), 2),
               rmse_knee=round(k.get("rmse_knee", float("nan")), 2),
               rmse_ankle=round(k.get("rmse_ankle", float("nan")), 2),
               cadence=round(k.get("n_cycles", 0) / 7.5, 3),
               tilt_max=round(m.get("tilt_max", float("nan")), 1),
               stayed_up=(not m.get("nan")) and m.get("kz", 0) > 0.62,
               nan=m.get("nan"), n_cycles=k.get("n_cycles"))
    row["_dt"] = round(time.time() - t0)
    print(f"({ge:.2f},{gf:.2f}) kine {row['kine_score']:>9} duty "
          f"{row['duty']} knee_min {row['knee_min']:>7} cad "
          f"{row['cadence']} tilt {row['tilt_max']} up={row['stayed_up']} "
          f"[{row['_dt']}s]", flush=True)
    return row


def main():
    batch = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    jobs = [(ge, gf) for ge in GRID for gf in GRID]
    rows = []
    with ThreadPoolExecutor(max_workers=batch) as ex:
        futures = [ex.submit(run_one, ge, gf) for ge, gf in jobs]
        rows = [f.result() for f in futures]
    with open(CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nwrote {CSV} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
