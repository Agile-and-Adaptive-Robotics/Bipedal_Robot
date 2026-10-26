# Full onset lists: base vs p1/p3, ext + flx, for the four-way verdict.
import numpy as np


def onsets(t, v, thr=-45.0):
    on, last = [], -1e9
    above = v > thr
    for i in range(1, len(v)):
        if above[i] and not above[i - 1] and t[i] - last > 0.2 and t[i] > 1.0:
            on.append(round(float(t[i]), 3)); last = t[i]
    return on


for cond in ("base", "p1_rg", "p2_pf", "p3_rg_long"):
    d = np.load(f"deng_mujoco_{cond}.npz")
    oE, oF = onsets(d["t"], d["vE"]), onsets(d["t"], d["vF"])
    per = np.diff(oE)
    print(f"{cond}: ext {oE}")
    print(f"{'':12s} flx {oF}")
    if len(per):
        print(f"{'':12s} ext periods {np.round(per, 4)}")
