"""FIGURES task (2026-09-25): replay the goal-4 VARIANT stage-5 winner.

Replicates the exact objective code path of _curriculum_{variant}.py
stage 5 (env pins, BASE_MUL merge from best_walk_params_v10.json,
set_stage, runner args, kine+kz+tilt scoring) and leaves the run's npz
at reports_20260925/figs/tmp/figs_npz_{variant}.npz for the figure
scripts. No campaign artifact is touched (AARL_NPZ is figs-private;
no optuna db is opened).

Usage:
    python replay_s5.py w2lvar|syn6
Exit 0 iff the recomputed score matches the json within 1%.
"""
import importlib
import json
import os
import sys
import time

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
SPINAL = r"D:\Github\Bipedal_Robot\Code\MuJoCo_SNS\spinal"
sys.path.insert(0, SPINAL)
os.chdir(SPINAL)

variant = sys.argv[1]
assert variant in ("w2lvar", "syn6")
FIGTMP = os.path.join("reports_20260925", "figs", "tmp")
os.makedirs(FIGTMP, exist_ok=True)
npz_rel = f"reports_20260925/figs/tmp/figs_npz_{variant}.npz"
if os.path.exists(npz_rel):
    os.remove(npz_rel)

# ---- replicate main()'s setup (env pins + BASE_MUL) --------------------
os.environ["AARL_NET"] = variant
os.environ["AARL_NPZ"] = npz_rel
prev = json.loads(open("best_walk_params_v10.json",
                       encoding="utf-8").read())
BASE_MUL = dict(prev["multipliers"])
BASE_MUL["renshaw"] = 0.5 if variant == "w2lvar" else 0.0
if variant == "syn6":
    BASE_MUL["syn6"] = 1.0
    BASE_MUL["syn6_brainstem"] = 0.0

cur = importlib.import_module(f"_curriculum_{variant}")
J = json.loads(open(f"curriculum_{variant}_stage5.json",
                    encoding="utf-8").read())
p = {**BASE_MUL, **J["params"]}
print(f"[figs:{variant}] json score={J['score']!r} trial={J['trial']} "
      f"study={J['study']}", flush=True)
print(f"[figs:{variant}] pinned AARL_NET={os.environ['AARL_NET']} "
      f"AARL_NPZ={os.environ['AARL_NPZ']}", flush=True)

cur.set_stage(5, p)

# ---- stage-5 runner invocation (exact stage-5 args) --------------------
args = ["--eval", "--drive", repr(p["drive"])]
print(f"[figs:{variant}] runner args: {args}", flush=True)
t0 = time.time()
m = cur.R.main(args)
wall = time.time() - t0
print(f"[figs:{variant}] runner wall time {wall:.1f} s", flush=True)
km = m.get("kine") or {}
print(f"[figs:{variant}] kine_score={m['kine_score']:.4f} kz={m['kz']:.4f} "
      f"tilt_max={m['tilt_max']:.2f} duty={m['duty']:.4f} "
      f"n_r={km.get('n_cycles_r')} n_l={km.get('n_cycles_l')} "
      f"T_r={km.get('T_r', float('nan')):.3f} "
      f"lag_rl={km.get('lag_rl', float('nan')):.3f} "
      f"bilateral={km.get('bilateral')}", flush=True)

# ---- exact stage-5 objective recompute ---------------------------------
score = max(float(m["kine_score"]), -315.0)
if m["kz"] < 0.62:
    score -= 20.0
if m["tilt_max"] > 40.0:
    score -= 10.0
jval = float(J["score"])
delta = score - jval
pct = abs(delta) / max(abs(jval), 1e-9) * 100.0
print(f"[figs:{variant}] RECOMPUTED {score!r} vs json {jval!r} "
      f"delta {delta:+.6f} ({pct:.4f}%)")
print(f"[figs:{variant}] REPLAY:", "PASS" if pct <= 1.0 else "FAIL")
print(f"[figs:{variant}] npz kept at {npz_rel}")
sys.exit(0 if pct <= 1.0 else 1)
