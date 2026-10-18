"""Orchestrator for the 2026-10-16 HUMAN L/R alpha+gamma drive sweep.

Grid (per the campaign phase-4 ask):
  gamma multiplier g in {1.5, 2.5, 4.0}  (x the tuned drive level; the
    runner feeds u[GAMMA_c_s] = drive*G[gamma_s], so g=1.5 -> 3.12 nA)
  x alpha a in {0.0, 0.3, 1.0}           (u[ALPHA_c_s] = drive*a)
  x L/R imbalance d in {0, +0.2, -0.2}   (left knobs *= (1+d),
                                          right knobs *= (1-d))
plus a knobs-0 baseline row (id "base") which must reproduce the anchor A
kine_score -161.56754173676563 bit-exactly (in-harness gate).

27 combos + baseline = 28. Each combo runs in a FRESH subprocess
(params.G is module-global); results append to
reports_20261016/drive_sweep_human/results.jsonl and already-done ids are
skipped on re-run (safe to chunk across ssh calls).

Usage: python drive_sweep_human_20261016.py [--limit N] [--only ID ...]
"""
import json
import os
import subprocess
import sys

PY = r"D:\Anaconda\envs\myo\python.exe"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "drive_sweep_human",
                   "results.jsonl")
os.makedirs(os.path.dirname(OUT), exist_ok=True)

combos = [dict(id="base", gamma_l=0.0, gamma_r=0.0,
               alpha_l=0.0, alpha_r=0.0)]
for g in (1.5, 2.5, 4.0):
    for a in (0.0, 0.3, 1.0):
        for d in (0.0, 0.2, -0.2):
            tag = "d0" if d == 0 else ("dp20" if d > 0 else "dm20")
            combos.append(dict(
                id=f"g{g:g}_a{a:g}_{tag}",
                gamma_l=round(g * (1 + d), 6),
                gamma_r=round(g * (1 - d), 6),
                alpha_l=round(a * (1 + d), 6),
                alpha_r=round(a * (1 - d), 6)))

done = set()
if os.path.exists(OUT):
    for ln in open(OUT, encoding="utf-8"):
        try:
            done.add(json.loads(ln)["spec"]["id"])
        except Exception:
            pass

args = sys.argv[1:]
limit, only = None, []
while args:
    x = args.pop(0)
    if x == "--limit":
        limit = int(args.pop(0))
    elif x == "--only":
        only.append(args.pop(0))
todo = [c for c in combos if c["id"] not in done
        and (not only or c["id"] in only)]
if limit is not None:
    todo = todo[:limit]
print(f"{len(combos)} combos total, {len(done)} already done, "
      f"running {len(todo)} now", flush=True)
env = dict(os.environ, PYTHONIOENCODING="utf-8")
for c in todo:
    print(f"=== combo {c['id']}: {json.dumps(c)} ===", flush=True)
    r = subprocess.run(
        [PY, os.path.join(HERE, "drive_sweep_one_20261016.py"),
         json.dumps(c)], cwd=HERE, env=env)
    if r.returncode != 0:
        print(f"COMBO {c['id']} FAILED rc={r.returncode}", flush=True)
print("ORCH DONE", flush=True)
