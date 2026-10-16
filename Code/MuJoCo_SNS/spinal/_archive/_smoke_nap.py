"""Smoke test: full runner, NaP RG architecture, deafferented air 14 s
(default-ish gains: rg_nap_h 0.35, weak G_W 0.4, no afferent-central).
Success = finite, rhythmic (>=1 rise), no crash in the class-swap
compile."""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace")

import numpy as np

import _curriculum as C
import params as P
import runner as R

import json
mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]
p = {**mul, "drive": 2.5, "rg_nap_h": 0.35, "desc_e": 1.7, "desc_f": 1.4,
     "rg_to_pf": 2.4}
C.set_stage(1, p)
R.main(["--no-ground", "--no-afferents", "--no-interleg", "--time", "14",
        "--drive", "2.5"])
z = np.load("spinal_run.npz", allow_pickle=True)
t, q, neuro = z["t"], z["q"], z["neuro"]
m = (t >= 5.0) & (t <= 17.0)
ok = bool(np.all(np.isfinite(q[m])) and np.all(np.isfinite(neuro[m])))
knee = q[m, 4]
rge = neuro[m, 2]
on = rge > 0.5 * max(rge.max(), 1e-9)
rises = int(np.sum(np.diff(on.astype(int)) == 1))
print(f"[smoke] finite={ok} rises={rises} knee_min={float(knee.min()):.1f} "
      f"RG_E range {float(rge.min()):.2f}..{float(rge.max()):.2f} mV")
assert ok and rises >= 1, "smoke FAILED"
print("SMOKE PASS")
