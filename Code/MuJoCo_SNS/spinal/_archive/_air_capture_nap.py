"""Dedicated AIR-WALKING capture on the NaP architecture (default-ish
gains, deafferented) - saved to nap_air_walk.npz so the curriculum's
per-trial spinal_run.npz overwrites don't clobber it. (No stdout wrap
here: _curriculum wraps it on import - double-wrapping closes it.)"""
import json
import shutil
import sys

import _curriculum as C
import params as P
import runner as R

mul = json.loads(open("best_walk_params_v10.json",
                      encoding="utf-8").read())["multipliers"]
p = {**mul, "drive": 2.5, "rg_nap_h": 0.35, "desc_e": 1.7, "desc_f": 1.4,
     "rg_to_pf": 2.4}
C.set_stage(1, p)
print(f"NaP air walk: drive 2.5, tau_h {P.TAU['rg_nap_h']:.2f} s, "
      f"G_W {P.G['rg_weak_exc']:.2f}", flush=True)
R.main(["--no-ground", "--no-afferents", "--no-interleg", "--time", "22",
        "--drive", "2.5"])
shutil.copy("spinal_run.npz", "nap_air_walk.npz")
print("saved nap_air_walk.npz")
