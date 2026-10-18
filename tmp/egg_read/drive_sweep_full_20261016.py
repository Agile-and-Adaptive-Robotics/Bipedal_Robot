"""Full 22 s (default schedule) run of a drive-sweep winner — the
"stays up 20 s" verification the 16 s --eval window cannot give.

Usage: D:\\Anaconda\\envs\\myo\\python.exe drive_sweep_full_20261016.py <knobs.json> [npzname]

knobs.json = {"id":..., "knobs": {gamma_l, gamma_r, alpha_l, alpha_r}}
Non-eval run: runner prints its honest summary (COM height / stayed up /
FELL at t, per-leg cycle metrics); output npz is redirected via AARL_NPZ
(default drive_sweep_winner_full.npz) so spinal_run.npz is untouched.
Post-hoc: reload the npz, print com_z min + stays-up over the FULL window
and the walk-window kine metrics via kine_ref.compare.
"""
import json
import os
import sys

knobs = json.loads(open(sys.argv[1], encoding="utf-8").read())
npz = sys.argv[2] if len(sys.argv) > 2 else "drive_sweep_winner_full.npz"

import _curriculum as C
import params as P
import runner as R

base = dict(json.loads(
    open("best_walk_params_v10.json", encoding="utf-8").read())
    ["multipliers"])
base["renshaw"] = 0.5
s3k = json.loads(open(
    "reports_20260923/s3k_trial34_full_params.json",
    encoding="utf-8").read())["params"]
base.update(s3k)
base["full_rules"] = 0.0          # anchor A recipe (sweep base)
C.set_stage(4, base)
for k in ("gamma_l", "gamma_r", "alpha_l", "alpha_r"):
    P.G[k] = float(knobs["knobs"][k])

os.environ["AARL_NPZ"] = npz
R.main(["--drive", repr(base["drive"])])

import numpy as np
d = np.load(npz, allow_pickle=True)
t, com, ct, q, neuro = d["t"], d["com"], d["contact"], d["q"], d["neuro"]
zmin = float(np.nanmin(com[:, 2]))
print(f"FULLRUN id={knobs.get('id')} npz={npz} t_end={float(t[-1]):.2f} "
      f"com_z_min={zmin:.3f} stays_up_20s={zmin > 0.55}", flush=True)
try:
    import kine_ref
    kk = kine_ref.compare(t, q, neuro, 5.0,
                          ref=kine_ref.ref_cached(), contact=ct)
    scal = {k2: v for k2, v in (kk or {}).items()
            if isinstance(v, (int, float, bool))}
    print("FULLRUN kine: " + json.dumps(scal), flush=True)
except Exception as e:
    print(f"FULLRUN kine-ref failed: {e}", flush=True)
