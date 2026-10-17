"""Full 22 s (default production schedule) verification run of a wean
round winner — the "stays up 20 s" check + production-window kine, the
same verification drive_sweep_full_20261016.py gave the bootstrap.

Usage: python wean_full_20261016.py <winner.json> <outnpz>
"""
import json
import os
import sys

winner = json.loads(open(sys.argv[1], encoding="utf-8").read())
npz = sys.argv[2] if len(sys.argv) > 2 else "wean_winner_full.npz"

import _curriculum as C
import params as P
import runner as R

base = dict(winner["params"])
C.set_stage(4, base)
scale = float(winner.get("scale", 1.0))
for k in ("gamma_l", "gamma_r", "alpha_l", "alpha_r"):
    P.G[k] = float(winner["knobs"][k]) * scale
P.G["alpha_to_mn"] = float(winner["knobs"].get("alpha_to_mn", 0.5))

os.environ["AARL_NPZ"] = npz
R.main(["--drive", repr(base["drive"])])

import numpy as np
d = np.load(npz, allow_pickle=True)
t, com, ct, q, neuro = d["t"], d["com"], d["contact"], d["q"], d["neuro"]
zmin = float(np.nanmin(com[:, 2]))
print(f"FULLRUN id=round{winner.get('round')} scale={scale:g} "
      f"npz={npz} t_end={float(t[-1]):.2f} com_z_min={zmin:.3f} "
      f"stays_up_20s={zmin > 0.55}", flush=True)
try:
    import kine_ref
    kk = kine_ref.compare(t, q, neuro, 5.0, ref=kine_ref.ref_cached(),
                          contact=ct)
    scal = {k2: v for k2, v in (kk or {}).items()
            if isinstance(v, (int, float, bool))}
    print("FULLRUN kine: " + json.dumps(scal), flush=True)
except Exception as e:
    print(f"FULLRUN kine-ref failed: {e}", flush=True)
