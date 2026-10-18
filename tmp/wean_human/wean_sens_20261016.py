"""Gain-sensitivity probe (Ben's anti-lottery requirement 2): perturb the
finalist's top-k tuned knobs one at a time — reset each to its SEED value
(the s3k/bootstrap incumbent it was tuned from) and, where the seed is
nonzero, also zero it — and confirm the score degrades meaningfully.
A gain whose reset changes nothing was lottery weight, not structure.

Usage: python wean_sens_20261016.py <winner.json> <baseline_score>
Appends variants to wean_human/sens.jsonl (fresh subprocess per variant,
3 concurrent). seed values = s3k winner params + bootstrap knobs.
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "reports_20261016", "wean_human")
PY = r"D:\Anaconda\envs\myo\python.exe"

win = json.loads(open(sys.argv[1], encoding="utf-8").read())
baseline = float(sys.argv[2])
TOPK = 8

# seed values: FB keys from the anchor base (s3k winner), drive knobs from
# the bootstrap sweep winner
import wean_tune_20261016 as WT
base = WT.load_base()
seed = {k: float(base.get(k, 0.0)) for k in WT.FB_KEYS}
seed.update(WT.BOOT_KNOBS)

finalist_scores = {
    **{k: float(win["params"].get(k, 0.0)) for k in WT.FB_KEYS},
    **{k: float(win["knobs"][k]) for k in WT.DRIVE_KEYS},
}

# rank by relative movement from seed (lottery suspects move ~0)
def rel(k):
    s0 = seed.get(k, 0.0)
    v = finalist_scores[k]
    if s0 == 0.0:
        return abs(v)          # seed-zero: any tuned value is movement
    return abs(v - s0) / max(abs(s0), 1e-9)

ranked = sorted(finalist_scores, key=rel, reverse=True)[:TOPK]
print("top-k by movement from seed: " +
      ", ".join(f"{k}={finalist_scores[k]:.3f} (seed {seed.get(k, 0):.3f})"
                for k in ranked))

variants = []
for k in ranked:
    variants.append((k, "seed", seed.get(k, 0.0)))
    if seed.get(k, 0.0) != 0.0:
        variants.append((k, "zero", 0.0))

jsonl = os.path.join(OUT, "sens.jsonl")
log = open(os.path.join(OUT, "sens.log"), "a", encoding="utf-8",
           buffering=1)
log.write(f"baseline {baseline!r} winner {sys.argv[1]}\n")

# each variant overrides ONE G key (drive knobs route through knobs/scale
# for gamma/alpha; FB keys override ON TOP of the winner's tuned fb, which
# must ride along in spec["fb"] — without it set_stage would load s3k
# values and the variant would not be a single-key perturbation)
win_fb = {k: float(win["params"].get(k, 0.0)) for k in WT.FB_KEYS}
tasks = []
for k, kind, val in variants:
    spec = dict(tag=f"sens_{k}_{kind}", knobs=dict(win["knobs"]),
                scale=float(win.get("scale", 1.0)), jsonl=jsonl,
                fb=dict(win_fb))
    if k in ("gamma_l", "gamma_r", "alpha_l", "alpha_r", "alpha_to_mn"):
        spec["knobs"][k] = val
    else:
        spec["overrides"] = {k: val}
    tasks.append((k, kind, val, spec))

running = []
pending = list(tasks)
while pending or running:
    while pending and len(running) < 3:
        k, kind, val, spec = pending.pop(0)
        sp = os.path.join(OUT, f"_sens_spec_{k}_{kind}.json")
        json.dump(spec, open(sp, "w", encoding="utf-8"))
        wl = open(os.path.join(OUT, f"_sens_{k}_{kind}.log"), "w",
                  encoding="utf-8", buffering=1)
        pr = subprocess.Popen([PY, "-u", os.path.join(
            HERE, "wean_tune_20261016.py"), "one", sp], stdout=wl,
            stderr=wl, cwd=HERE,
            env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        running.append((k, kind, val, pr, wl))
        log.write(f"launched {k} {kind}={val:g} pid {pr.pid}\n")
        time.sleep(5)
    time.sleep(10)
    for item in list(running):
        k, kind, val, pr, wl = item
        if pr.poll() is not None:
            wl.close()
            running.remove(item)
            log.write(f"{k} {kind} EXIT {pr.returncode}\n")
log.write("SENS_DONE\n")
log.close()
print("sens rows in " + jsonl)
