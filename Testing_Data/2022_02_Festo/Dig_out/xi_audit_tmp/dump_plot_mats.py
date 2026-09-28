# dump_plot_mats.py - torque curves in the two optimized plot mats (dtype-guarded)
import scipy.io as sio
import numpy as np
import os, datetime

BASE = r"D:\GitHub\Bipedal_Robot\Testing_Data\2022_02_Festo"

def mt(p):
    return datetime.datetime.fromtimestamp(os.path.getmtime(p)).strftime("%Y-%m-%d %H:%M")

for f, lbl in [("KneeFlx_20mm_optimized.mat", "FLEXOR"), ("KneeExt_20mm_optimized.mat", "EXTENSOR")]:
    p = os.path.join(BASE, f)
    m = sio.loadmat(p, squeeze_me=True, struct_as_record=False)
    print(f"\n=== {lbl}: {f} (mtime {mt(p)}) ===")
    ang = np.asarray(m['Angle']).flatten()
    print(f" Angle: [{ang.min():.2f}, {ang.max():.2f}] N={ang.size}")
    print(f" {'var':<14} {'min':>9} {'max':>9}")
    for k in sorted(m.keys()):
        if k.startswith('__') or k == 'Angle' or k == 'None':
            continue
        v = m[k]
        if hasattr(v, '_fieldnames') or hasattr(v, 'dtype') and v.dtype == object:
            continue
        try:
            a = np.asarray(v, dtype=float).flatten()
        except Exception:
            continue
        if a.size == ang.size:
            print(f" {k:<14} {np.nanmin(a):9.4f} {np.nanmax(a):9.4f}")
        elif a.size <= 6:
            print(f" {k:<14} = {a}")
    # candidate ratio table: robot/human torque curves
    pairs = {
        "FLEXOR": [("Bifemsh_T", None), ("H", None)],
        "EXTENSOR": [("Torque2", "Torque0"), ("Torque3", "Torque0"), ("G2", None), ("G3", None)],
    }.get(lbl, [])
    names = [n for n, _ in pairs]
    cols = {}
    for n in names:
        if n in m:
            try:
                cols[n] = np.asarray(m[n], dtype=float).flatten()
            except Exception:
                pass
    if len(cols) >= 2:
        print(f" {'ang':>7} " + " ".join(f"{n:>10}" for n in cols))
        for i in range(ang.size):
            print(f" {ang[i]:7.1f} " + " ".join(f"{cols[n][i]:10.3f}" for n in cols))
