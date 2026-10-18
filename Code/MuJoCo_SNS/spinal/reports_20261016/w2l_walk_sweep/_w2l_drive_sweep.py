r"""W2L ROBOT drive-bootstrap sweep (overnight campaign 2026-10-16, egg).

Stage: "bootstrap driver (robot)". Sweeps strong symmetric + asymmetric
L/R GAMMA (te/tf per side via --te-l/-r, --tf-l/-r) and ALPHA
(--alpha-l/-r) drive on test_w2l_ground.py --phase=walk under the
drop-protocol harness/virtual-walker (hold 3 s -> lower 1 s -> retained
wsup support; rig S leash/tilt per combo).

Usage (egg, from anywhere):
    D:\Anaconda\envs\myo\python.exe _w2l_drive_sweep.py <combos.json>

combos.json: list of {"id": str, "args": {knob: value, ...}} — knob
names use underscores, converted to --dashes CLI args of
test_w2l_ground.py. Every run's full stdout goes to logs/<id>.log; one
CSV row per run appended to results.csv (created with header on first
write). Sequential runs; ~7.5 s per 20 s sim.
"""
from __future__ import annotations

import csv
import io
import json
import os
import re
import subprocess
import sys
import time

PY = r"D:\Anaconda\envs\myo\python.exe"
TEST = (r"D:\GitHub\Bipedal_Robot\Code\MuJoCo_SNS\spinal\w2l_mujoco"
        r"\test_w2l_ground.py")
HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, "logs")
os.makedirs(LOGS, exist_ok=True)
CSV_PATH = os.path.join(HERE, "results.csv")

FIELDS = ["id", "args", "verdict", "finite", "fall", "window_s",
          "strikes_L", "strikes_R", "per_L", "per_R", "duty_L", "duty_R",
          "fwd_m", "v_ms", "tilt_max", "pelz_min", "harness_pct",
          "rg_bursts_L", "rg_bursts_R", "rg_per", "rg_r",
          "heelSN", "toeSN", "heel_in",
          "hip_L", "hip_R", "knee_L", "knee_R", "ank_L", "ank_R",
          "wall_s"]

RE_HEAD = re.compile(
    r"finite=(\w+) \| fall=(t=[\d.]+ s|no) \| pelvis z min ([\d.]+) m "
    r"\(start [\d.]+\d*\) \| tilt max ([\d.\-+]+) deg \| harness carries "
    r"([+\-][\d.]+)%")
RE_CYC = re.compile(
    r"cycles \(heel strikes >20 N, window ([\d.]+) s\): L (\d+) "
    r"\(interval ([\d.na]+) s\) \| R (\d+) \(interval ([\d.na]+) s\)")
RE_DUTY = re.compile(
    r"heel duty L ([\d.]+) R ([\d.]+) \| forward disp ([+\-][\d.]+) m "
    r"\(([+\-][\d.]+) m/s\)")
RE_RNG = re.compile(
    r"flexion ranges \(deg\): hip L ([\d.]+) R ([\d.]+) \| knee L "
    r"([\d.]+) R ([\d.]+) \| ankle L ([\d.]+) R ([\d.]+)")
RE_NEU = re.compile(
    r"neural: L RG-E bursts (\d+) \(period ([\d.na]+) s\), R RG-E bursts "
    r"(\d+); L/R RG-E r ([+\-.\dna]+)")
RE_ENC = re.compile(
    r"contact encoder evidence: heel SN max ([\d.]+) mV, toe SN max "
    r"([\d.]+) mV, heel port current max ([\d.]+) nA")
RE_VER = re.compile(r"walk verdict: (.+)")


def fmt(v):
    if isinstance(v, float) and v != v:
        return "nan"
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def run_one(cid: str, args: dict) -> dict:
    cmd = [PY, TEST, "--phase=walk"]
    for k, v in args.items():
        cmd.append(f"--{k.replace('_', '-')}={fmt(v)}")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600)
    wall = time.time() - t0
    out = p.stdout
    open(os.path.join(LOGS, f"{cid}.log"), "w", encoding="utf-8").write(
        out + "\n--- stderr ---\n" + p.stderr)
    row = dict(id=cid, args=json.dumps(args, sort_keys=True), verdict="PARSE-FAIL",
               finite="", fall="", window_s="", strikes_L="", strikes_R="",
               per_L="", per_R="", duty_L="", duty_R="", fwd_m="", v_ms="",
               tilt_max="", pelz_min="", harness_pct="", rg_bursts_L="",
               rg_bursts_R="", rg_per="", rg_r="", heelSN="", toeSN="",
               heel_in="", hip_L="", hip_R="", knee_L="", knee_R="",
               ank_L="", ank_R="", wall_s=f"{wall:.1f}")
    m = RE_HEAD.search(out)
    if m:
        row["finite"], row["fall"], row["pelz_min"], row["tilt_max"], \
            row["harness_pct"] = (m.group(1), m.group(2), m.group(3),
                                  m.group(4), m.group(5))
    m = RE_CYC.search(out)
    if m:
        (row["window_s"], row["strikes_L"], row["per_L"], row["strikes_R"],
         row["per_R"]) = m.groups()
    m = RE_DUTY.search(out)
    if m:
        row["duty_L"], row["duty_R"], row["fwd_m"], row["v_ms"] = m.groups()
    m = RE_RNG.search(out)
    if m:
        (row["hip_L"], row["hip_R"], row["knee_L"], row["knee_R"],
         row["ank_L"], row["ank_R"]) = m.groups()
    m = RE_NEU.search(out)
    if m:
        (row["rg_bursts_L"], row["rg_per"], row["rg_bursts_R"],
         row["rg_r"]) = m.groups()
    m = RE_ENC.search(out)
    if m:
        row["heelSN"], row["toeSN"], row["heel_in"] = m.groups()
    m = RE_VER.search(out)
    if m:
        row["verdict"] = m.group(1).strip()
    return row


def main() -> int:
    combos = json.load(open(sys.argv[1], encoding="utf-8"))
    new = not os.path.exists(CSV_PATH)
    fh = open(CSV_PATH, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()
    for c in combos:
        row = run_one(c["id"], c["args"])
        w.writerow(row)
        fh.flush()
        print(f"[{row['id']}] {row['verdict'][:26]:26s} "
              f"strikes {row['strikes_L']}/{row['strikes_R']} "
              f"duty {row['duty_L']}/{row['duty_R']} "
              f"fwd {row['fwd_m']} bursts {row['rg_bursts_L']}/"
              f"{row['rg_bursts_R']} r {row['rg_r']} "
              f"tilt {row['tilt_max']} ({row['wall_s']}s)", flush=True)
    fh.close()
    print("DONE", len(combos), "combos ->", CSV_PATH, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
